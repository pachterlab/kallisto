#ifndef GPU_READ_LOADER_CUH
#define GPU_READ_LOADER_CUH

// GPU-accelerated FASTQ file reader
// Supports BGZF (GPU decompression via nvcomp) and gzip (multithreaded rapidgzip CPU decompression)

#include <cuda_runtime.h>
#include <nvcomp.h>
#include <nvcomp/deflate.h>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/execution_policy.h>

#include "cufq/common/bgzf.cuh"
#include "cufq/common/fastq_parser.cuh"
#include "cufq/common/utils.cuh"
#include "RapidGzipReader.h"
#include "common.h"
#include "BenchmarkStats.h"

#include <memory>

// ============================================================================
// Macros
// ============================================================================

#define NVCOMP_CHECK(call)                                                     \
    do {                                                                       \
        nvcompStatus_t status = call;                                          \
        if (status != nvcompSuccess) {                                         \
            fprintf(stderr, "nvcomp error at %s:%d: %d\n", __FILE__, __LINE__, \
                    (int)status);                                              \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)

// ============================================================================
// BatchDecompressor: Pre-allocated GPU buffers for nvcomp batched deflate
// ============================================================================
class BatchDecompressor {
public:
    BatchDecompressor(size_t max_batch_size, cudaStream_t stream)
        : stream_(stream), max_batch_size_(max_batch_size) {

        // BGZF guarantees max 64KB decompressed per block
        max_blocks_ = max_batch_size / (64 * 1024) + 1;

        // Pre-allocate device arrays for metadata
        CUDA_CHECK(cudaMalloc(&d_comp_ptrs_, max_blocks_ * sizeof(void*)));
        CUDA_CHECK(cudaMalloc(&d_decomp_ptrs_, max_blocks_ * sizeof(void*)));
        CUDA_CHECK(cudaMalloc(&d_comp_sizes_, max_blocks_ * sizeof(size_t)));
        CUDA_CHECK(cudaMalloc(&d_decomp_sizes_, max_blocks_ * sizeof(size_t)));
        CUDA_CHECK(cudaMalloc(&d_statuses_, max_blocks_ * sizeof(nvcompStatus_t)));
        CUDA_CHECK(cudaMalloc(&d_actual_sizes_, max_blocks_ * sizeof(size_t)));

        // Pre-allocate host pinned memory for metadata
        CUDA_CHECK(cudaMallocHost(&h_comp_ptrs_, max_blocks_ * sizeof(void*)));
        CUDA_CHECK(cudaMallocHost(&h_decomp_ptrs_, max_blocks_ * sizeof(void*)));
        CUDA_CHECK(cudaMallocHost(&h_comp_sizes_, max_blocks_ * sizeof(size_t)));
        CUDA_CHECK(cudaMallocHost(&h_decomp_sizes_, max_blocks_ * sizeof(size_t)));

        // Pre-allocate compressed data buffer on GPU
        CUDA_CHECK(cudaMalloc(&d_compressed_, max_batch_size));

        // Pre-allocate temp buffer
        nvcompBatchedDeflateDecompressOpts_t opts = nvcompBatchedDeflateDecompressDefaultOpts;
        NVCOMP_CHECK(nvcompBatchedDeflateDecompressGetTempSizeAsync(
            max_blocks_, 65536, opts, &temp_size_, max_batch_size));
        CUDA_CHECK(cudaMalloc(&d_temp_, temp_size_));
    }

    ~BatchDecompressor() {
        cudaFree(d_comp_ptrs_);
        cudaFree(d_decomp_ptrs_);
        cudaFree(d_comp_sizes_);
        cudaFree(d_decomp_sizes_);
        cudaFree(d_statuses_);
        cudaFree(d_actual_sizes_);
        cudaFree(d_compressed_);
        cudaFree(d_temp_);
        cudaFreeHost(h_comp_ptrs_);
        cudaFreeHost(h_decomp_ptrs_);
        cudaFreeHost(h_comp_sizes_);
        cudaFreeHost(h_decomp_sizes_);
    }

    size_t max_blocks() const { return max_blocks_; }

    // Decompress a batch of BGZF blocks using single bulk transfer
    // Returns total decompressed size
    size_t decompress_batch(
        const char* h_file_data,
        const std::vector<BgzfBlock>& blocks,
        size_t start_idx,
        size_t num_blocks,
        char* d_decompressed
    ) {
        if (num_blocks == 0) return 0;
        if (num_blocks > max_blocks_) {
            fprintf(stderr, "Error: num_blocks %zu exceeds max_blocks %zu\n",
                    num_blocks, max_blocks_);
            return 0;
        }

        // Calculate contiguous file range for this batch
        size_t range_start = blocks[start_idx].file_offset;
        size_t range_end = blocks[start_idx + num_blocks - 1].file_offset +
                          blocks[start_idx + num_blocks - 1].compressed_size;
        size_t range_size = range_end - range_start;

        if (range_size > max_batch_size_) {
            fprintf(stderr, "Error: Batch range %zu exceeds buffer %zu\n",
                    range_size, max_batch_size_);
            return 0;
        }

        // Single bulk H2D transfer
        CUDA_CHECK(cudaMemcpyAsync(d_compressed_, h_file_data + range_start,
                                   range_size, cudaMemcpyHostToDevice, stream_));

        // Compute pointers within the GPU buffer
        size_t decomp_offset = 0;
        size_t batch_uncompressed = 0;

        for (size_t i = 0; i < num_blocks; i++) {
            const BgzfBlock& block = blocks[start_idx + i];

            uint16_t xlen;
            memcpy(&xlen, h_file_data + block.file_offset + 10, sizeof(uint16_t));
            size_t header_size = 10 + 2 + xlen;
            size_t trailer_size = 8;

            size_t offset_in_range = block.file_offset - range_start + header_size;
            size_t deflate_size = block.compressed_size - header_size - trailer_size;

            h_comp_ptrs_[i] = d_compressed_ + offset_in_range;
            h_decomp_ptrs_[i] = d_decompressed + decomp_offset;
            h_comp_sizes_[i] = deflate_size;
            h_decomp_sizes_[i] = block.uncompressed_size;

            decomp_offset += block.uncompressed_size;
            batch_uncompressed += block.uncompressed_size;
        }

        // Transfer metadata
        CUDA_CHECK(cudaMemcpyAsync(d_comp_ptrs_, h_comp_ptrs_,
                                   num_blocks * sizeof(void*),
                                   cudaMemcpyHostToDevice, stream_));
        CUDA_CHECK(cudaMemcpyAsync(d_decomp_ptrs_, h_decomp_ptrs_,
                                   num_blocks * sizeof(void*),
                                   cudaMemcpyHostToDevice, stream_));
        CUDA_CHECK(cudaMemcpyAsync(d_comp_sizes_, h_comp_sizes_,
                                   num_blocks * sizeof(size_t),
                                   cudaMemcpyHostToDevice, stream_));
        CUDA_CHECK(cudaMemcpyAsync(d_decomp_sizes_, h_decomp_sizes_,
                                   num_blocks * sizeof(size_t),
                                   cudaMemcpyHostToDevice, stream_));

        // GPU decompression
        nvcompBatchedDeflateDecompressOpts_t opts = nvcompBatchedDeflateDecompressDefaultOpts;
        NVCOMP_CHECK(nvcompBatchedDeflateDecompressAsync(
            (const void* const*)d_comp_ptrs_,
            d_comp_sizes_,
            d_decomp_sizes_,
            d_actual_sizes_,
            num_blocks,
            d_temp_,
            temp_size_,
            (void* const*)d_decomp_ptrs_,
            opts,
            d_statuses_,
            stream_));

        return batch_uncompressed;
    }

private:
    cudaStream_t stream_;
    size_t max_batch_size_;
    size_t max_blocks_;
    size_t temp_size_;

    char* d_compressed_;
    void** d_comp_ptrs_;
    void** d_decomp_ptrs_;
    size_t* d_comp_sizes_;
    size_t* d_decomp_sizes_;
    nvcompStatus_t* d_statuses_;
    size_t* d_actual_sizes_;
    void* d_temp_;

    void** h_comp_ptrs_;
    void** h_decomp_ptrs_;
    size_t* h_comp_sizes_;
    size_t* h_decomp_sizes_;
};

// ============================================================================
// GPU kernel: count newlines in data
// ============================================================================
inline __global__ void gpu_read_count_newlines_kernel(
    const char* __restrict__ data,
    size_t size,
    uint32_t* count
) {
    __shared__ uint32_t block_count;
    if (threadIdx.x == 0) block_count = 0;
    __syncthreads();

    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t local_count = 0;

    for (size_t i = idx; i < size; i += blockDim.x * gridDim.x) {
        if (data[i] == '\n') local_count++;
    }

    atomicAdd(&block_count, local_count);
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicAdd(count, block_count);
    }
}

// ============================================================================
// GPU kernel: find record boundary (last complete FASTQ record)
// ============================================================================
inline __global__ void gpu_read_find_record_boundary_kernel(
    const char* __restrict__ data,
    size_t size,
    uint32_t skip_lines,
    size_t* result_pos
) {
    if (blockIdx.x != 0 || threadIdx.x != 0) return;

    if (skip_lines == 0) {
        *result_pos = size;
        return;
    }

    uint32_t newlines_to_find = skip_lines + 1;
    uint32_t newlines_found = 0;
    size_t pos = size;

    while (pos > 0 && newlines_found < newlines_to_find) {
        pos--;
        if (data[pos] == '\n') {
            newlines_found++;
        }
    }

    if (newlines_found == newlines_to_find) {
        *result_pos = pos + 1;
    } else {
        *result_pos = 0;
    }
}

// ============================================================================
// GPU kernel: walk back from a SEQ-line offset to find the '@' header start
// of that record. Used to slice surplus parsed records into a leftover buffer.
// ============================================================================
inline __global__ void gpu_read_find_record_header_start_kernel(
    const char* __restrict__ data,
    uint64_t seq_offset,
    size_t* __restrict__ out_pos
) {
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    if (seq_offset < 2) {
        *out_pos = 0;
        return;
    }
    // data[seq_offset-1] is '\n' that ends the '@' header line.
    // Walk back from seq_offset-2 to find the previous '\n' (end of prev qual)
    // or beginning of buffer. The byte immediately after that newline is '@'.
    size_t pos = seq_offset - 2;
    while (pos > 0 && data[pos] != '\n') pos--;
    *out_pos = (data[pos] == '\n') ? pos + 1 : pos;
}

// ============================================================================
// GPU kernel: convert SeqDescriptor + base_offset to read metadata arrays
// Preserves ALL reads in order (short reads get kmer_count=0) so that
// R1[i] and R2[i] remain paired. valid_count tracks stats only.
// ============================================================================
inline __global__ void convert_descriptors_kernel(
    const SeqDescriptor* __restrict__ descriptors,
    uint32_t num_sequences,
    uint64_t base_offset,
    uint32_t k,
    uint32_t start_idx,
    uint64_t* __restrict__ out_offsets,
    uint32_t* __restrict__ out_lengths,
    uint32_t* __restrict__ out_kmer_counts,
    uint32_t* __restrict__ valid_count
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_sequences) return;

    uint32_t len = descriptors[idx].length;
    uint32_t out_idx = start_idx + idx;
    out_offsets[out_idx] = base_offset + descriptors[idx].offset;
    out_lengths[out_idx] = len;
    if (len >= k) {
        out_kmer_counts[out_idx] = len - k + 1;
        atomicAdd(valid_count, 1);
    } else {
        out_kmer_counts[out_idx] = 0;
    }
}

// Functor to cast uint32_t to uint64_t (used with thrust::transform_iterator)
struct CastU32ToU64 {
    __host__ __device__ uint64_t operator()(uint32_t c) const {
        return static_cast<uint64_t>(c);
    }
};

// ============================================================================
// GPUReadLoader: main class
// ============================================================================
struct GPUReadLoader {
    // === Output data (on device) - consumed by DeviceKmerLoader ===
    char* d_reads;                           // Raw FASTQ data on GPU (single allocation for R1+R2)
    size_t d_reads_size;                     // Total size of d_reads
    thrust::device_vector<uint64_t> read_id_to_offset;
    thrust::device_vector<uint32_t> read_id_to_length;
    thrust::device_vector<uint64_t> read_id_to_kmer_first;
    thrust::device_vector<uint32_t> read_id_to_kmer_count;
    uint64_t kmer_count;
    uint32_t read_count;
    uint32_t r1_count;                       // Number of valid R1 reads (R2 starts at this index)

    // === Constructor ===
    GPUReadLoader(const ProgramOptions& opt, size_t batch_size_bytes = 256 * 1024 * 1024)
        : d_reads(nullptr), d_reads_size(0), kmer_count(0), read_count(0), r1_count(0),
          k_(Kmer::k), batch_size_(batch_size_bytes),
          initialized_(false), finished_(false), is_bgzf_(false),
          num_files_(opt.files.size()),
          block_idx_r1_(0), block_idx_r2_(0),
          partial_size_r1_(0), partial_size_r2_(0),
          h_file_r1_(nullptr), h_file_r2_(nullptr),
          file_size_r1_(0), file_size_r2_(0),
          d_decomp_r1_(nullptr), d_decomp_r2_(nullptr),
          d_partial_r1_(nullptr), d_partial_r2_(nullptr),
          d_valid_count_(nullptr),
          d_rb_count_r1_(nullptr), d_rb_boundary_r1_(nullptr),
          d_rb_count_r2_(nullptr), d_rb_boundary_r2_(nullptr),
          decomp_r1_(nullptr), decomp_r2_(nullptr),
          h_gz_buf_r1_(nullptr), h_gz_buf_r2_(nullptr),
          d_leftover_r1_(nullptr), d_leftover_r2_(nullptr),
          leftover_size_r1_(0), leftover_size_r2_(0),
          leftover_capacity_r1_(0), leftover_capacity_r2_(0)
    {
        files_ = opt.files;

        CUDA_CHECK(cudaStreamCreate(&stream_r1_));
        CUDA_CHECK(cudaStreamCreate(&stream_r2_));

        // Pre-allocate temp buffers for record boundary detection (separate per file)
        CUDA_CHECK(cudaMalloc(&d_rb_count_r1_, sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_rb_boundary_r1_, sizeof(size_t)));
        CUDA_CHECK(cudaMalloc(&d_rb_count_r2_, sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_rb_boundary_r2_, sizeof(size_t)));
    }

    ~GPUReadLoader() {
        // Free GPU buffers
        if (d_reads) cudaFree(d_reads);
        if (d_decomp_r1_) cudaFree(d_decomp_r1_);
        if (d_decomp_r2_) cudaFree(d_decomp_r2_);
        if (d_partial_r1_) cudaFree(d_partial_r1_);
        if (d_partial_r2_) cudaFree(d_partial_r2_);
        if (d_valid_count_) cudaFree(d_valid_count_);

        // Free pinned host memory
        if (h_file_r1_) cudaFreeHost(h_file_r1_);
        if (h_file_r2_) cudaFreeHost(h_file_r2_);
        if (h_gz_buf_r1_) cudaFreeHost(h_gz_buf_r1_);
        if (h_gz_buf_r2_) cudaFreeHost(h_gz_buf_r2_);

        // Free decompressors
        delete decomp_r1_;
        delete decomp_r2_;

        rapidgzip_r1_.reset();
        rapidgzip_r2_.reset();

        if (d_rb_count_r1_) cudaFree(d_rb_count_r1_);
        if (d_rb_boundary_r1_) cudaFree(d_rb_boundary_r1_);
        if (d_rb_count_r2_) cudaFree(d_rb_count_r2_);
        if (d_rb_boundary_r2_) cudaFree(d_rb_boundary_r2_);

        if (d_leftover_r1_) cudaFree(d_leftover_r1_);
        if (d_leftover_r2_) cudaFree(d_leftover_r2_);

        cudaStreamDestroy(stream_r1_);
        cudaStreamDestroy(stream_r2_);
    }

    // === Initialize: read files, detect format, set up decompressors ===
    void initialize() {
        if (initialized_) return;
        initialized_ = true;

        if (num_files_ == 0) {
            std::cerr << "Error: No input files specified" << std::endl;
            exit(1);
        }

        // Detect format from file header (avoid reading entire file for gzip)
        is_bgzf_ = is_bgzf_file(files_[0]);

        if (!is_bgzf_) {
            std::cerr << "[warning] Input files are not BGZF format. "
                      << "Using multithreaded rapidgzip CPU decompression. "
                      << "For GPU decompression use: bgzip -c file.fastq > file.fastq.bgz"
                      << std::endl;
        }

        const size_t PARTIAL_MAX = 64 * 1024;

        if (is_bgzf_) {
            read_file_to_pinned(files_[0], &h_file_r1_, &file_size_r1_);

            if (num_files_ >= 2) {
                read_file_to_pinned(files_[1], &h_file_r2_, &file_size_r2_);
            }
        }

        if (is_bgzf_) {
            // Parse BGZF block headers
            blocks_r1_ = parse_bgzf_blocks_from_memory(h_file_r1_, file_size_r1_);

            if (num_files_ >= 2) {
                blocks_r2_ = parse_bgzf_blocks_from_memory(h_file_r2_, file_size_r2_);
            }

            // Create batch decompressors
            decomp_r1_ = new BatchDecompressor(batch_size_, stream_r1_);
            if (num_files_ >= 2) {
                decomp_r2_ = new BatchDecompressor(batch_size_, stream_r2_);
            }

            // Allocate GPU decompression buffers
            CUDA_CHECK(cudaMalloc(&d_decomp_r1_, batch_size_ + PARTIAL_MAX));
            CUDA_CHECK(cudaMalloc(&d_partial_r1_, PARTIAL_MAX));
            if (num_files_ >= 2) {
                CUDA_CHECK(cudaMalloc(&d_decomp_r2_, batch_size_ + PARTIAL_MAX));
                CUDA_CHECK(cudaMalloc(&d_partial_r2_, PARTIAL_MAX));
            }
        } else {
            // gzip: multithreaded rapidgzip decompression
            rapidgzip_r1_ = std::make_unique<RapidGzipReader>(files_[0], 0, 4ULL * 1024 * 1024);

            CUDA_CHECK(cudaMallocHost(&h_gz_buf_r1_, batch_size_));

            CUDA_CHECK(cudaMalloc(&d_decomp_r1_, batch_size_ + PARTIAL_MAX));
            CUDA_CHECK(cudaMalloc(&d_partial_r1_, PARTIAL_MAX));

            if (num_files_ >= 2) {
                rapidgzip_r2_ = std::make_unique<RapidGzipReader>(files_[1], 0, 4ULL * 1024 * 1024);

                CUDA_CHECK(cudaMallocHost(&h_gz_buf_r2_, batch_size_));
                CUDA_CHECK(cudaMalloc(&d_decomp_r2_, batch_size_ + PARTIAL_MAX));
                CUDA_CHECK(cudaMalloc(&d_partial_r2_, PARTIAL_MAX));
            }
        }

        // Allocate output buffer for combined R1+R2 reads
        // We'll allocate on first use based on actual decompressed size
        // For now, allocate a generous buffer
        size_t combined_size = (num_files_ >= 2) ? 2 * (batch_size_ + PARTIAL_MAX)
                                                  : (batch_size_ + PARTIAL_MAX);
        CUDA_CHECK(cudaMalloc(&d_reads, combined_size));
        d_reads_size = combined_size;

        // Allocate valid count device variable
        CUDA_CHECK(cudaMalloc(&d_valid_count_, sizeof(uint32_t)));

        // Create FASTQ parsers and pre-allocate work buffers
        parser_r1_ = std::make_unique<FastqParser>(batch_size_ + PARTIAL_MAX);
        parser_r1_->ensure_work_buffers(batch_size_ + PARTIAL_MAX);
        if (num_files_ >= 2) {
            parser_r2_ = std::make_unique<FastqParser>(batch_size_ + PARTIAL_MAX);
            parser_r2_->ensure_work_buffers(batch_size_ + PARTIAL_MAX);
        }
    }

    // === Load next batch ===
    // Returns true if data was loaded, false if all files are exhausted
    bool load() {
        auto start_time = std::chrono::high_resolution_clock::now();

        if (!initialized_) {
            initialize();
        }

        if (finished_) return false;

        read_count = 0;
        kmer_count = 0;

        bool got_data = false;

        if (is_bgzf_) {
            got_data = load_bgzf_batch();
        } else {
            got_data = load_gzip_batch();
        }

        if (!got_data) {
            finished_ = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            g_benchmark_stats.io_decompress_ms += duration.count() / 1000.0;
            return false;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        g_benchmark_stats.io_decompress_ms += duration.count() / 1000.0;

        return true;
    }

private:
    // === Helper: read file into pinned memory ===
    static void read_file_to_pinned(const std::string& path, char** out_data, size_t* out_size) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f.is_open()) {
            std::cerr << "Error: Cannot open file: " << path << std::endl;
            exit(1);
        }
        *out_size = static_cast<size_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        CUDA_CHECK(cudaMallocHost(out_data, *out_size));
        f.read(*out_data, *out_size);
        f.close();
    }

    // === Record boundary detection (split into launch/finish for R1/R2 parallelism) ===

    // Launch newline counting kernel (async, returns immediately)
    void launch_newline_count(const char* d_data, size_t size, cudaStream_t stream,
                              uint32_t* d_count) {
        CUDA_CHECK(cudaMemsetAsync(d_count, 0, sizeof(uint32_t), stream));
        int block_size = 256;
        int num_blocks = std::min((int)((size + block_size - 1) / block_size), 256);
        gpu_read_count_newlines_kernel<<<num_blocks, block_size, 0, stream>>>(d_data, size, d_count);
    }

    // Sync stream, read newline count, launch boundary kernel if needed, return boundary
    size_t finish_record_boundary(const char* d_data, size_t size, cudaStream_t stream,
                                  uint32_t* d_count, size_t* d_boundary) {
        uint32_t total_newlines;
        CUDA_CHECK(cudaMemcpyAsync(&total_newlines, d_count, sizeof(uint32_t),
                                    cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));

        uint32_t extra_lines = total_newlines % 4;
        if (extra_lines == 0) return size;

        gpu_read_find_record_boundary_kernel<<<1, 1, 0, stream>>>(
            d_data, size, extra_lines, d_boundary);

        size_t boundary;
        CUDA_CHECK(cudaMemcpyAsync(&boundary, d_boundary, sizeof(size_t),
                                    cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));

        return boundary;
    }

    // === Helper: convert parsed sequences to read metadata arrays (GPU-based) ===
    // Uses convert_descriptors_kernel + thrust::exclusive_scan instead of D2H round-trip
    void build_read_metadata_gpu(
        FastqParser& parser,
        uint32_t num_seqs,
        uint64_t base_offset,
        uint32_t start_read_idx,
        uint64_t& running_kmer_count,
        cudaStream_t stream
    ) {
        if (num_seqs == 0) return;

        SeqDescriptor* d_descriptors = parser.device_descriptors();

        // Reset valid count
        CUDA_CHECK(cudaMemsetAsync(d_valid_count_, 0, sizeof(uint32_t), stream));

        // Filter reads < k and compute offset/length/kmer_count on GPU
        int threads = 256;
        int blocks = (num_seqs + threads - 1) / threads;
        convert_descriptors_kernel<<<blocks, threads, 0, stream>>>(
            d_descriptors, num_seqs, base_offset, k_, start_read_idx,
            thrust::raw_pointer_cast(read_id_to_offset.data()),
            thrust::raw_pointer_cast(read_id_to_length.data()),
            thrust::raw_pointer_cast(read_id_to_kmer_count.data()),
            d_valid_count_
        );

        // Synchronize to ensure kernel is complete before scanning
        CUDA_CHECK(cudaStreamSynchronize(stream));

        // Compute kmer_first as exclusive prefix sum of kmer_count for ALL reads
        // (including short reads with kmer_count=0) to preserve index alignment for pairing
        auto kmer_count_begin = read_id_to_kmer_count.begin() + start_read_idx;
        auto kmer_first_begin = read_id_to_kmer_first.begin() + start_read_idx;
        auto count_as_u64 = thrust::make_transform_iterator(
            kmer_count_begin, CastU32ToU64()
        );
        thrust::exclusive_scan(
            thrust::cuda::par.on(stream),
            count_as_u64, count_as_u64 + num_seqs,
            kmer_first_begin,
            running_kmer_count
        );

        // Compute total kmer count with a reduce
        uint64_t total_new_kmers = thrust::reduce(
            thrust::cuda::par.on(stream),
            count_as_u64, count_as_u64 + num_seqs,
            0ULL
        );
        running_kmer_count += total_new_kmers;

        read_count = start_read_idx + num_seqs;
        kmer_count = running_kmer_count;
    }

    // === BGZF batch loading ===
    bool load_bgzf_batch() {
        const size_t PARTIAL_MAX = 64 * 1024;

        // Check if R1 is done
        bool r1_blocks_done = (block_idx_r1_ >= blocks_r1_.size());
        bool r2_blocks_done = (num_files_ < 2) || (block_idx_r2_ >= blocks_r2_.size());
        bool r1_done = r1_blocks_done && (leftover_size_r1_ == 0);
        bool r2_done = r2_blocks_done && (leftover_size_r2_ == 0);

        if (r1_done && r2_done) return false;

        size_t max_blocks_limit = decomp_r1_->max_blocks();

        // === Collect R1 blocks (cap so leftover + new_decomp <= batch_size_) ===
        size_t r1_budget = (leftover_size_r1_ < batch_size_)
                               ? (batch_size_ - leftover_size_r1_)
                               : 0;
        size_t batch_start_r1 = block_idx_r1_;
        size_t batch_uncomp_r1 = 0;
        size_t batch_block_count_r1 = 0;
        while (block_idx_r1_ < blocks_r1_.size() &&
               batch_uncomp_r1 + blocks_r1_[block_idx_r1_].uncompressed_size <= r1_budget &&
               batch_block_count_r1 < max_blocks_limit) {
            batch_uncomp_r1 += blocks_r1_[block_idx_r1_].uncompressed_size;
            block_idx_r1_++;
            batch_block_count_r1++;
        }
        size_t num_blocks_r1 = block_idx_r1_ - batch_start_r1;

        // === Collect R2 blocks (cap so leftover + new_decomp <= batch_size_) ===
        size_t batch_start_r2 = block_idx_r2_;
        size_t num_blocks_r2 = 0;
        if (num_files_ >= 2 && !r2_blocks_done) {
            size_t r2_budget = (leftover_size_r2_ < batch_size_)
                                   ? (batch_size_ - leftover_size_r2_)
                                   : 0;
            size_t batch_uncomp_r2 = 0;
            size_t batch_block_count_r2 = 0;
            while (block_idx_r2_ < blocks_r2_.size() &&
                   batch_uncomp_r2 + blocks_r2_[block_idx_r2_].uncompressed_size <= r2_budget &&
                   batch_block_count_r2 < max_blocks_limit) {
                batch_uncomp_r2 += blocks_r2_[block_idx_r2_].uncompressed_size;
                block_idx_r2_++;
                batch_block_count_r2++;
            }
            num_blocks_r2 = block_idx_r2_ - batch_start_r2;
        }

        if (num_blocks_r1 == 0 && num_blocks_r2 == 0 &&
            leftover_size_r1_ == 0 && leftover_size_r2_ == 0) {
            return false;
        }

        bool is_last_batch_r1 = (block_idx_r1_ >= blocks_r1_.size());
        bool is_last_batch_r2 = (num_files_ < 2) || (block_idx_r2_ >= blocks_r2_.size());

        // === Copy leftover (surplus parsed records from previous batch) at the
        //     start of d_decomp, then partial fragment, then new decompressed.
        if (leftover_size_r1_ > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_decomp_r1_, d_leftover_r1_, leftover_size_r1_,
                                       cudaMemcpyDeviceToDevice, stream_r1_));
        }
        if (partial_size_r1_ > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_decomp_r1_ + leftover_size_r1_, d_partial_r1_,
                                       partial_size_r1_, cudaMemcpyDeviceToDevice,
                                       stream_r1_));
        }
        if (num_files_ >= 2 && leftover_size_r2_ > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_decomp_r2_, d_leftover_r2_, leftover_size_r2_,
                                       cudaMemcpyDeviceToDevice, stream_r2_));
        }
        if (num_files_ >= 2 && partial_size_r2_ > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_decomp_r2_ + leftover_size_r2_, d_partial_r2_,
                                       partial_size_r2_, cudaMemcpyDeviceToDevice,
                                       stream_r2_));
        }

        // === Decompress ===
        size_t new_decomp_size_r1 = 0;
        if (num_blocks_r1 > 0) {
            new_decomp_size_r1 = decomp_r1_->decompress_batch(
                h_file_r1_, blocks_r1_, batch_start_r1, num_blocks_r1,
                d_decomp_r1_ + leftover_size_r1_ + partial_size_r1_);
        }

        size_t new_decomp_size_r2 = 0;
        if (num_blocks_r2 > 0) {
            new_decomp_size_r2 = decomp_r2_->decompress_batch(
                h_file_r2_, blocks_r2_, batch_start_r2, num_blocks_r2,
                d_decomp_r2_ + leftover_size_r2_ + partial_size_r2_);
        }

        size_t total_size_r1 = leftover_size_r1_ + partial_size_r1_ + new_decomp_size_r1;
        size_t total_size_r2 = leftover_size_r2_ + partial_size_r2_ + new_decomp_size_r2;
        // The leftover and partial buffers have been consumed for this batch.
        leftover_size_r1_ = 0;
        leftover_size_r2_ = 0;

        // Sync both streams
        CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
        if (num_files_ >= 2) {
            CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
        }

        // === Find record boundaries (R1/R2 newline counting overlaps on GPU) ===
        size_t parse_size_r1 = total_size_r1;
        size_t parse_size_r2 = total_size_r2;
        bool need_boundary_r1 = !is_last_batch_r1 && total_size_r1 > 0;
        bool need_boundary_r2 = num_files_ >= 2 && !is_last_batch_r2 && total_size_r2 > 0;

        // Launch newline counting kernels on both streams (overlap on GPU)
        if (need_boundary_r1) {
            launch_newline_count(d_decomp_r1_, total_size_r1, stream_r1_, d_rb_count_r1_);
        }
        if (need_boundary_r2) {
            launch_newline_count(d_decomp_r2_, total_size_r2, stream_r2_, d_rb_count_r2_);
        }

        // Finish R1 boundary detection (sync + possible second kernel)
        if (need_boundary_r1) {
            parse_size_r1 = finish_record_boundary(d_decomp_r1_, total_size_r1, stream_r1_,
                                                   d_rb_count_r1_, d_rb_boundary_r1_);
            size_t new_partial = total_size_r1 - parse_size_r1;
            if (new_partial > 0 && new_partial <= PARTIAL_MAX) {
                CUDA_CHECK(cudaMemcpyAsync(d_partial_r1_, d_decomp_r1_ + parse_size_r1,
                                           new_partial, cudaMemcpyDeviceToDevice, stream_r1_));
            } else if (new_partial > PARTIAL_MAX) {
                new_partial = 0;
            }
            partial_size_r1_ = new_partial;
        } else {
            partial_size_r1_ = 0;
        }

        // Finish R2 boundary detection
        if (need_boundary_r2) {
            parse_size_r2 = finish_record_boundary(d_decomp_r2_, total_size_r2, stream_r2_,
                                                   d_rb_count_r2_, d_rb_boundary_r2_);
            size_t new_partial = total_size_r2 - parse_size_r2;
            if (new_partial > 0 && new_partial <= PARTIAL_MAX) {
                CUDA_CHECK(cudaMemcpyAsync(d_partial_r2_, d_decomp_r2_ + parse_size_r2,
                                           new_partial, cudaMemcpyDeviceToDevice, stream_r2_));
            } else if (new_partial > PARTIAL_MAX) {
                new_partial = 0;
            }
            partial_size_r2_ = new_partial;
        } else {
            partial_size_r2_ = 0;
        }

        // === Copy decompressed data into combined reads buffer ===
        // d_reads = [R1 data][R2 data]
        size_t needed = parse_size_r1 + parse_size_r2;
        if (needed > d_reads_size) {
            cudaFree(d_reads);
            CUDA_CHECK(cudaMalloc(&d_reads, needed));
            d_reads_size = needed;
        }
        if (parse_size_r1 > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_reads, d_decomp_r1_, parse_size_r1,
                                       cudaMemcpyDeviceToDevice, stream_r1_));
        }
        if (parse_size_r2 > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_reads + parse_size_r1, d_decomp_r2_, parse_size_r2,
                                       cudaMemcpyDeviceToDevice, stream_r2_));
        }
        CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
        if (num_files_ >= 2) {
            CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
        }

        // === Parse FASTQ on GPU ===
        uint32_t seqs_r1 = 0, seqs_r2 = 0;
        if (parse_size_r1 > 0) {
            seqs_r1 = parser_r1_->parse_device_data(d_reads, parse_size_r1, 0, stream_r1_);
        }
        if (parse_size_r2 > 0) {
            seqs_r2 = parser_r2_->parse_device_data(d_reads + parse_size_r1, parse_size_r2, 0, stream_r2_);
        }

        // === Build read metadata on GPU ===
        uint32_t total_seqs = seqs_r1 + seqs_r2;
        if (total_seqs == 0) return false;

        // Pre-allocate metadata arrays (upper bound; may have fewer valid reads)
        read_id_to_offset.resize(total_seqs);
        read_id_to_length.resize(total_seqs);
        read_id_to_kmer_first.resize(total_seqs);
        read_id_to_kmer_count.resize(total_seqs);

        uint64_t running_kmer_count = 0;
        read_count = 0;

        // Build R1 metadata on GPU
        build_read_metadata_gpu(*parser_r1_, seqs_r1, 0, 0, running_kmer_count, stream_r1_);
        uint32_t r1_valid = read_count;

        // Build R2 metadata on GPU (offsets adjusted by parse_size_r1)
        if (seqs_r2 > 0) {
            build_read_metadata_gpu(*parser_r2_, seqs_r2, parse_size_r1, r1_valid,
                                   running_kmer_count, stream_r2_);
        }

        uint32_t r2_valid = read_count - r1_valid;

        // Paired-end equalization: truncate to min(r1, r2) so pairs match.
        // Surplus parsed records on the larger side are saved as a leftover
        // buffer (their bytes from '@' header onwards) so that the next batch
        // can re-parse them. This preserves all reads instead of dropping them.
        if (num_files_ >= 2 && r1_valid != r2_valid) {
            uint32_t pair_count = std::min(r1_valid, r2_valid);

            // Save R1 surplus bytes for next batch
            if (pair_count > 0 && r1_valid > pair_count) {
                size_t boundary = find_record_start_offset(
                    d_decomp_r1_, *parser_r1_, pair_count,
                    d_rb_boundary_r1_, stream_r1_);
                size_t surplus = (total_size_r1 > boundary)
                                     ? (total_size_r1 - boundary)
                                     : 0;
                grow_leftover(&d_leftover_r1_, &leftover_capacity_r1_, surplus);
                if (surplus > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(d_leftover_r1_,
                                               d_decomp_r1_ + boundary,
                                               surplus,
                                               cudaMemcpyDeviceToDevice,
                                               stream_r1_));
                    CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
                }
                leftover_size_r1_ = surplus;
                partial_size_r1_ = 0;
            }
            // Save R2 surplus bytes for next batch
            if (pair_count > 0 && r2_valid > pair_count) {
                size_t boundary = find_record_start_offset(
                    d_decomp_r2_, *parser_r2_, pair_count,
                    d_rb_boundary_r2_, stream_r2_);
                size_t surplus = (total_size_r2 > boundary)
                                     ? (total_size_r2 - boundary)
                                     : 0;
                grow_leftover(&d_leftover_r2_, &leftover_capacity_r2_, surplus);
                if (surplus > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(d_leftover_r2_,
                                               d_decomp_r2_ + boundary,
                                               surplus,
                                               cudaMemcpyDeviceToDevice,
                                               stream_r2_));
                    CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
                }
                leftover_size_r2_ = surplus;
                partial_size_r2_ = 0;
            }

            if (pair_count < r1_valid) {
                thrust::copy(
                    read_id_to_offset.begin() + r1_valid,
                    read_id_to_offset.begin() + r1_valid + pair_count,
                    read_id_to_offset.begin() + pair_count);
                thrust::copy(
                    read_id_to_length.begin() + r1_valid,
                    read_id_to_length.begin() + r1_valid + pair_count,
                    read_id_to_length.begin() + pair_count);
                thrust::copy(
                    read_id_to_kmer_first.begin() + r1_valid,
                    read_id_to_kmer_first.begin() + r1_valid + pair_count,
                    read_id_to_kmer_first.begin() + pair_count);
                thrust::copy(
                    read_id_to_kmer_count.begin() + r1_valid,
                    read_id_to_kmer_count.begin() + r1_valid + pair_count,
                    read_id_to_kmer_count.begin() + pair_count);
            }
            read_count = 2 * pair_count;
            r1_valid = pair_count;
            running_kmer_count = 0;
            if (read_count > 0) {
                uint64_t last_kmer_first = 0;
                uint32_t last_kmer_cnt = 0;
                CUDA_CHECK(cudaMemcpy(&last_kmer_first,
                    thrust::raw_pointer_cast(read_id_to_kmer_first.data()) + read_count - 1,
                    sizeof(uint64_t), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&last_kmer_cnt,
                    thrust::raw_pointer_cast(read_id_to_kmer_count.data()) + read_count - 1,
                    sizeof(uint32_t), cudaMemcpyDeviceToHost));
                running_kmer_count = last_kmer_first + last_kmer_cnt;
            }
        }

        r1_count = r1_valid;
        kmer_count = running_kmer_count;

        // Trim arrays to actual size
        read_id_to_offset.resize(read_count);
        read_id_to_length.resize(read_count);
        read_id_to_kmer_first.resize(read_count);
        read_id_to_kmer_count.resize(read_count);

        return read_count > 0;
    }

    // === gzip batch loading (CPU decompression + GPU parsing) ===
    bool load_gzip_batch() {
        const size_t PARTIAL_MAX = 64 * 1024;

        bool r1_eof = !rapidgzip_r1_ || rapidgzip_r1_->eof();
        bool r2_eof = (num_files_ < 2) || !rapidgzip_r2_ || rapidgzip_r2_->eof();

        bool r1_done = r1_eof && partial_size_r1_ == 0 && leftover_size_r1_ == 0;
        bool r2_done = r2_eof && partial_size_r2_ == 0 && leftover_size_r2_ == 0;
        if (r1_done && r2_done) return false;

        // === CPU decompress R1 ===
        size_t parse_size_r1 = 0;
        if (!r1_eof || partial_size_r1_ > 0 || leftover_size_r1_ > 0) {
            // Copy leftover (surplus parsed records from previous batch) + partial
            if (leftover_size_r1_ > 0) {
                CUDA_CHECK(cudaMemcpyAsync(d_decomp_r1_, d_leftover_r1_, leftover_size_r1_,
                                           cudaMemcpyDeviceToDevice, stream_r1_));
            }
            if (partial_size_r1_ > 0) {
                CUDA_CHECK(cudaMemcpyAsync(d_decomp_r1_ + leftover_size_r1_, d_partial_r1_,
                                           partial_size_r1_, cudaMemcpyDeviceToDevice,
                                           stream_r1_));
            }

            // Cap read budget so leftover + partial + new <= batch_size_
            size_t r1_budget = 0;
            if (leftover_size_r1_ + partial_size_r1_ < batch_size_) {
                r1_budget = batch_size_ - leftover_size_r1_ - partial_size_r1_;
            }
            int bytes_r1 = 0;
            if (!r1_eof && r1_budget > 0) {
                ssize_t n = rapidgzip_r1_->read(h_gz_buf_r1_, r1_budget);
                bytes_r1 = (n > 0) ? static_cast<int>(n) : 0;
            }

            if (bytes_r1 > 0) {
                CUDA_CHECK(cudaMemcpyAsync(
                    d_decomp_r1_ + leftover_size_r1_ + partial_size_r1_,
                    h_gz_buf_r1_, bytes_r1, cudaMemcpyHostToDevice, stream_r1_));
            }

            size_t total_r1 = leftover_size_r1_ + partial_size_r1_ + bytes_r1;
            CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
            leftover_size_r1_ = 0;

            bool is_last_r1 = rapidgzip_r1_->eof();
            if (!is_last_r1 && total_r1 > 0) {
                launch_newline_count(d_decomp_r1_, total_r1, stream_r1_, d_rb_count_r1_);
                parse_size_r1 = finish_record_boundary(d_decomp_r1_, total_r1, stream_r1_,
                                                       d_rb_count_r1_, d_rb_boundary_r1_);
                size_t new_partial = total_r1 - parse_size_r1;
                if (new_partial > 0 && new_partial <= PARTIAL_MAX) {
                    CUDA_CHECK(cudaMemcpyAsync(d_partial_r1_, d_decomp_r1_ + parse_size_r1,
                                               new_partial, cudaMemcpyDeviceToDevice, stream_r1_));
                } else if (new_partial > PARTIAL_MAX) {
                    new_partial = 0;
                }
                partial_size_r1_ = new_partial;
            } else {
                parse_size_r1 = total_r1;
                partial_size_r1_ = 0;
            }
        }

        // === CPU decompress R2 ===
        size_t parse_size_r2 = 0;
        if (num_files_ >= 2 &&
            (!r2_eof || partial_size_r2_ > 0 || leftover_size_r2_ > 0)) {
            if (leftover_size_r2_ > 0) {
                CUDA_CHECK(cudaMemcpyAsync(d_decomp_r2_, d_leftover_r2_, leftover_size_r2_,
                                           cudaMemcpyDeviceToDevice, stream_r2_));
            }
            if (partial_size_r2_ > 0) {
                CUDA_CHECK(cudaMemcpyAsync(d_decomp_r2_ + leftover_size_r2_, d_partial_r2_,
                                           partial_size_r2_, cudaMemcpyDeviceToDevice,
                                           stream_r2_));
            }

            size_t r2_budget = 0;
            if (leftover_size_r2_ + partial_size_r2_ < batch_size_) {
                r2_budget = batch_size_ - leftover_size_r2_ - partial_size_r2_;
            }
            int bytes_r2 = 0;
            if (!r2_eof && r2_budget > 0) {
                ssize_t n = rapidgzip_r2_->read(h_gz_buf_r2_, r2_budget);
                bytes_r2 = (n > 0) ? static_cast<int>(n) : 0;
            }

            if (bytes_r2 > 0) {
                CUDA_CHECK(cudaMemcpyAsync(
                    d_decomp_r2_ + leftover_size_r2_ + partial_size_r2_,
                    h_gz_buf_r2_, bytes_r2, cudaMemcpyHostToDevice, stream_r2_));
            }

            size_t total_r2 = leftover_size_r2_ + partial_size_r2_ + bytes_r2;
            CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
            leftover_size_r2_ = 0;

            bool is_last_r2 = rapidgzip_r2_->eof();
            if (!is_last_r2 && total_r2 > 0) {
                launch_newline_count(d_decomp_r2_, total_r2, stream_r2_, d_rb_count_r2_);
                parse_size_r2 = finish_record_boundary(d_decomp_r2_, total_r2, stream_r2_,
                                                       d_rb_count_r2_, d_rb_boundary_r2_);
                size_t new_partial = total_r2 - parse_size_r2;
                if (new_partial > 0 && new_partial <= PARTIAL_MAX) {
                    CUDA_CHECK(cudaMemcpyAsync(d_partial_r2_, d_decomp_r2_ + parse_size_r2,
                                               new_partial, cudaMemcpyDeviceToDevice, stream_r2_));
                } else if (new_partial > PARTIAL_MAX) {
                    new_partial = 0;
                }
                partial_size_r2_ = new_partial;
            } else {
                parse_size_r2 = total_r2;
                partial_size_r2_ = 0;
            }
        }

        if (parse_size_r1 == 0 && parse_size_r2 == 0) return false;

        // === Copy into combined reads buffer ===
        size_t needed = parse_size_r1 + parse_size_r2;
        if (needed > d_reads_size) {
            cudaFree(d_reads);
            CUDA_CHECK(cudaMalloc(&d_reads, needed));
            d_reads_size = needed;
        }
        if (parse_size_r1 > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_reads, d_decomp_r1_, parse_size_r1,
                                       cudaMemcpyDeviceToDevice, stream_r1_));
        }
        if (parse_size_r2 > 0) {
            CUDA_CHECK(cudaMemcpyAsync(d_reads + parse_size_r1, d_decomp_r2_, parse_size_r2,
                                       cudaMemcpyDeviceToDevice, stream_r2_));
        }
        CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
        if (num_files_ >= 2) {
            CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
        }

        // === Parse FASTQ on GPU ===
        uint32_t seqs_r1 = 0, seqs_r2 = 0;
        if (parse_size_r1 > 0) {
            seqs_r1 = parser_r1_->parse_device_data(d_reads, parse_size_r1, 0, stream_r1_);
        }
        if (parse_size_r2 > 0) {
            seqs_r2 = parser_r2_->parse_device_data(d_reads + parse_size_r1, parse_size_r2, 0, stream_r2_);
        }

        // === Build read metadata on GPU ===
        uint32_t total_seqs = seqs_r1 + seqs_r2;
        if (total_seqs == 0) return false;

        read_id_to_offset.resize(total_seqs);
        read_id_to_length.resize(total_seqs);
        read_id_to_kmer_first.resize(total_seqs);
        read_id_to_kmer_count.resize(total_seqs);

        uint64_t running_kmer_count = 0;
        read_count = 0;

        build_read_metadata_gpu(*parser_r1_, seqs_r1, 0, 0, running_kmer_count, stream_r1_);
        uint32_t r1_valid = read_count;

        if (seqs_r2 > 0) {
            build_read_metadata_gpu(*parser_r2_, seqs_r2, parse_size_r1, r1_valid,
                                   running_kmer_count, stream_r2_);
        }

        uint32_t r2_valid = read_count - r1_valid;

        // Paired-end equalization: truncate to min(r1, r2) so pairs match.
        // Surplus parsed records on the larger side are saved as a leftover
        // buffer so the next batch can re-parse them (preserving all reads).
        if (num_files_ >= 2 && r1_valid != r2_valid) {
            uint32_t pair_count = std::min(r1_valid, r2_valid);

            // Total bytes in d_decomp_r1_ for this batch (parsed + new partial)
            size_t total_size_r1 = parse_size_r1 + partial_size_r1_;
            size_t total_size_r2 = parse_size_r2 + partial_size_r2_;

            if (pair_count > 0 && r1_valid > pair_count) {
                size_t boundary = find_record_start_offset(
                    d_decomp_r1_, *parser_r1_, pair_count,
                    d_rb_boundary_r1_, stream_r1_);
                size_t surplus = (total_size_r1 > boundary)
                                     ? (total_size_r1 - boundary)
                                     : 0;
                grow_leftover(&d_leftover_r1_, &leftover_capacity_r1_, surplus);
                if (surplus > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(d_leftover_r1_,
                                               d_decomp_r1_ + boundary,
                                               surplus,
                                               cudaMemcpyDeviceToDevice,
                                               stream_r1_));
                    CUDA_CHECK(cudaStreamSynchronize(stream_r1_));
                }
                leftover_size_r1_ = surplus;
                partial_size_r1_ = 0;
            }
            if (pair_count > 0 && r2_valid > pair_count) {
                size_t boundary = find_record_start_offset(
                    d_decomp_r2_, *parser_r2_, pair_count,
                    d_rb_boundary_r2_, stream_r2_);
                size_t surplus = (total_size_r2 > boundary)
                                     ? (total_size_r2 - boundary)
                                     : 0;
                grow_leftover(&d_leftover_r2_, &leftover_capacity_r2_, surplus);
                if (surplus > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(d_leftover_r2_,
                                               d_decomp_r2_ + boundary,
                                               surplus,
                                               cudaMemcpyDeviceToDevice,
                                               stream_r2_));
                    CUDA_CHECK(cudaStreamSynchronize(stream_r2_));
                }
                leftover_size_r2_ = surplus;
                partial_size_r2_ = 0;
            }

            // Rebuild metadata arrays: keep first pair_count R1 reads and first pair_count R2 reads
            if (pair_count < r1_valid) {
                thrust::copy(
                    read_id_to_offset.begin() + r1_valid,
                    read_id_to_offset.begin() + r1_valid + pair_count,
                    read_id_to_offset.begin() + pair_count);
                thrust::copy(
                    read_id_to_length.begin() + r1_valid,
                    read_id_to_length.begin() + r1_valid + pair_count,
                    read_id_to_length.begin() + pair_count);
                thrust::copy(
                    read_id_to_kmer_first.begin() + r1_valid,
                    read_id_to_kmer_first.begin() + r1_valid + pair_count,
                    read_id_to_kmer_first.begin() + pair_count);
                thrust::copy(
                    read_id_to_kmer_count.begin() + r1_valid,
                    read_id_to_kmer_count.begin() + r1_valid + pair_count,
                    read_id_to_kmer_count.begin() + pair_count);
            }
            read_count = 2 * pair_count;
            r1_valid = pair_count;
            running_kmer_count = 0;
            if (read_count > 0) {
                uint64_t last_kmer_first = 0;
                uint32_t last_kmer_cnt = 0;
                CUDA_CHECK(cudaMemcpy(&last_kmer_first,
                    thrust::raw_pointer_cast(read_id_to_kmer_first.data()) + read_count - 1,
                    sizeof(uint64_t), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&last_kmer_cnt,
                    thrust::raw_pointer_cast(read_id_to_kmer_count.data()) + read_count - 1,
                    sizeof(uint32_t), cudaMemcpyDeviceToHost));
                running_kmer_count = last_kmer_first + last_kmer_cnt;
            }
        }

        r1_count = r1_valid;
        kmer_count = running_kmer_count;

        read_id_to_offset.resize(read_count);
        read_id_to_length.resize(read_count);
        read_id_to_kmer_first.resize(read_count);
        read_id_to_kmer_count.resize(read_count);

        return read_count > 0;
    }

    // === Private members ===
    uint32_t k_;
    size_t batch_size_;
    bool initialized_;
    bool finished_;
    bool is_bgzf_;
    size_t num_files_;
    std::vector<std::string> files_;

    // BGZF state
    std::vector<BgzfBlock> blocks_r1_, blocks_r2_;
    size_t block_idx_r1_, block_idx_r2_;
    size_t partial_size_r1_, partial_size_r2_;
    char* h_file_r1_;
    char* h_file_r2_;
    size_t file_size_r1_, file_size_r2_;
    char* d_decomp_r1_;
    char* d_decomp_r2_;
    char* d_partial_r1_;
    char* d_partial_r2_;
    BatchDecompressor* decomp_r1_;
    BatchDecompressor* decomp_r2_;

    // gzip state (multithreaded rapidgzip)
    std::unique_ptr<RapidGzipReader> rapidgzip_r1_;
    std::unique_ptr<RapidGzipReader> rapidgzip_r2_;
    char* h_gz_buf_r1_;
    char* h_gz_buf_r2_;

    // Parsers
    std::unique_ptr<FastqParser> parser_r1_;
    std::unique_ptr<FastqParser> parser_r2_;

    // CUDA streams
    cudaStream_t stream_r1_;
    cudaStream_t stream_r2_;

    // Temp device memory
    uint32_t* d_valid_count_;

    // Pre-allocated temp buffers for record boundary detection (per-file for parallelism)
    uint32_t* d_rb_count_r1_;
    size_t* d_rb_boundary_r1_;
    uint32_t* d_rb_count_r2_;
    size_t* d_rb_boundary_r2_;

    // Leftover surplus parsed bytes carried over to the next batch when one
    // side ends up with more parsed records than the other. The contents are
    // [start of '@' header of first surplus record .. end of decompressed buffer]
    // (i.e. they include any cross-block-boundary partial fragment too, so when
    // leftover is non-zero the partial buffer is logically empty for that side).
    char* d_leftover_r1_;
    char* d_leftover_r2_;
    size_t leftover_size_r1_;
    size_t leftover_size_r2_;
    size_t leftover_capacity_r1_;
    size_t leftover_capacity_r2_;

    void grow_leftover(char** dptr, size_t* capacity, size_t needed) {
        if (*capacity >= needed) return;
        if (*dptr) cudaFree(*dptr);
        size_t new_cap = needed;
        // Round up to 1MiB granularity to avoid frequent reallocations
        const size_t align = 1024 * 1024;
        new_cap = ((new_cap + align - 1) / align) * align;
        CUDA_CHECK(cudaMalloc(dptr, new_cap));
        *capacity = new_cap;
    }

    // Find byte offset of the '@' header start for record `record_idx` in d_buf,
    // given the device-side seq-offset array (uint64_t per descriptor).
    // Issues a tiny kernel; sync on `stream` afterwards.
    size_t find_record_start_offset(
        const char* d_buf,
        const FastqParser& parser,
        uint32_t record_idx,
        size_t* d_out_scratch,
        cudaStream_t stream
    ) {
        // Copy descriptors[record_idx].offset to host
        SeqDescriptor desc;
        CUDA_CHECK(cudaMemcpyAsync(
            &desc,
            const_cast<FastqParser&>(parser).device_descriptors() + record_idx,
            sizeof(SeqDescriptor), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
        if (desc.offset == 0) return 0;
        // Tiny single-thread kernel walks back to '@' header start.
        gpu_read_find_record_header_start_kernel<<<1, 1, 0, stream>>>(
            d_buf, static_cast<uint64_t>(desc.offset), d_out_scratch);
        size_t out_pos = 0;
        CUDA_CHECK(cudaMemcpyAsync(&out_pos, d_out_scratch, sizeof(size_t),
                                   cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
        return out_pos;
    }
};

#endif // GPU_READ_LOADER_CUH
