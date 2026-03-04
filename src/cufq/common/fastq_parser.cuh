#pragma once

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <vector>
#include <cstring>
#include "utils.cuh"

// Kernel to mark newline positions (1 for '\n', 0 otherwise)
__global__ void mark_newlines_kernel(const char* __restrict__ data,
                                      uint32_t* __restrict__ marks,
                                      size_t size) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        marks[idx] = (data[idx] == '\n') ? 1 : 0;
    }
}

// Kernel to find sequence line starts and ends
// After inclusive scan, scan[i] = number of newlines up to and including position i
// Line type = (scan[i] - 1) % 4 for the character AFTER the newline
// We want line type 1 (sequence line)
__global__ void find_sequence_positions_kernel(
    const char* __restrict__ data,
    const uint32_t* __restrict__ scan,
    uint32_t* __restrict__ seq_starts,
    uint32_t* __restrict__ seq_count,
    size_t size,
    uint32_t newline_offset  // offset for batch continuation
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < size && data[idx] == '\n') {
        // This is a newline character
        // The line that ENDS here has type: (scan[idx] - 1 + newline_offset) % 4
        // Line 0 = read name, Line 1 = sequence, Line 2 = +, Line 3 = quality
        uint32_t line_type = (scan[idx] - 1 + newline_offset) % 4;
        
        // If this newline ends line 0 (read name), the next line is the sequence
        if (line_type == 0 && idx + 1 < size) {
            // The position after this newline is the start of a sequence
            uint32_t seq_idx = atomicAdd(seq_count, 1);
            seq_starts[seq_idx] = idx + 1;
        }
    }
}

// Kernel to compute sequence lengths (distance to next newline)
__global__ void compute_sequence_lengths_kernel(
    const char* __restrict__ data,
    const uint32_t* __restrict__ seq_starts,
    SeqDescriptor* __restrict__ descriptors,
    uint32_t num_sequences,
    size_t data_size
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_sequences) {
        uint32_t start = seq_starts[idx];
        uint32_t length = 0;
        
        // Find the next newline
        for (uint32_t pos = start; pos < data_size; pos++) {
            if (data[pos] == '\n') {
                length = pos - start;
                break;
            }
        }
        
        descriptors[idx].offset = start;
        descriptors[idx].length = length;
    }
}

// Alternative: vectorized length computation using scan results
__global__ void compute_sequence_lengths_fast_kernel(
    const uint32_t* __restrict__ scan,
    const uint32_t* __restrict__ seq_starts,
    SeqDescriptor* __restrict__ descriptors,
    uint32_t num_sequences,
    size_t data_size,
    uint32_t newline_offset
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_sequences) {
        uint32_t start = seq_starts[idx];
        descriptors[idx].offset = start;
        
        // Binary search for the next newline using the scan array
        // scan[pos] increases by 1 at each newline
        // We want the smallest pos > start where scan[pos] > scan[start-1]
        uint32_t target = (start > 0) ? scan[start - 1] + 1 : 1;
        
        uint32_t left = start;
        uint32_t right = min((uint32_t)data_size - 1, start + 500);  // Assume max seq len 500
        
        while (left < right) {
            uint32_t mid = (left + right) / 2;
            if (scan[mid] < target) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        descriptors[idx].length = left - start;
    }
}

// Main FASTQ parser class
class FastqParser {
public:
    FastqParser(size_t max_batch_size = DEFAULT_BATCH_SIZE)
        : max_batch_size_(max_batch_size),
          d_data_(nullptr),
          d_marks_(nullptr),
          d_scan_(nullptr),
          d_seq_starts_(nullptr),
          d_seq_count_(nullptr),
          d_descriptors_(nullptr),
          d_temp_storage_(nullptr),
          temp_storage_bytes_(0),
          allocated_size_(0) {
        CUDA_CHECK(cudaStreamCreate(&stream_));
    }
    
    ~FastqParser() {
        free_buffers();
        cudaStreamDestroy(stream_);
    }
    
    // Allocate GPU buffers for a given data size
    void allocate_buffers(size_t size) {
        if (size <= allocated_size_) return;
        
        free_buffers();
        
        CUDA_CHECK(cudaMalloc(&d_data_, size));
        CUDA_CHECK(cudaMalloc(&d_marks_, size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_scan_, size * sizeof(uint32_t)));
        
        // Estimate max sequences (one per 4 newlines, ~100 bytes per read minimum)
        size_t max_sequences = size / 50;
        CUDA_CHECK(cudaMalloc(&d_seq_starts_, max_sequences * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_seq_count_, sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_descriptors_, max_sequences * sizeof(SeqDescriptor)));
        
        // Determine temp storage for CUB
        CUDA_CHECK(cub::DeviceScan::InclusiveSum(
            nullptr, temp_storage_bytes_,
            d_marks_, d_scan_, size, stream_));
        CUDA_CHECK(cudaMalloc(&d_temp_storage_, temp_storage_bytes_));
        
        allocated_size_ = size;
    }
    
    // Pre-allocate work buffers (marks, scan, etc.) without allocating d_data_
    // Call this during initialization to avoid cudaMalloc during hot path
    void ensure_work_buffers(size_t size) {
        if (size <= allocated_size_) return;
        
        if (d_marks_) cudaFree(d_marks_);
        if (d_scan_) cudaFree(d_scan_);
        if (d_seq_starts_) cudaFree(d_seq_starts_);
        if (d_seq_count_) cudaFree(d_seq_count_);
        if (d_descriptors_) cudaFree(d_descriptors_);
        if (d_temp_storage_) cudaFree(d_temp_storage_);
        
        CUDA_CHECK(cudaMalloc(&d_marks_, size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_scan_, size * sizeof(uint32_t)));
        
        size_t max_sequences = size / 50;
        CUDA_CHECK(cudaMalloc(&d_seq_starts_, max_sequences * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_seq_count_, sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_descriptors_, max_sequences * sizeof(SeqDescriptor)));
        
        CUDA_CHECK(cub::DeviceScan::InclusiveSum(
            nullptr, temp_storage_bytes_,
            d_marks_, d_scan_, size, stream_));
        CUDA_CHECK(cudaMalloc(&d_temp_storage_, temp_storage_bytes_));
        
        allocated_size_ = size;
    }
    
    void free_buffers() {
        if (d_data_) cudaFree(d_data_);
        if (d_marks_) cudaFree(d_marks_);
        if (d_scan_) cudaFree(d_scan_);
        if (d_seq_starts_) cudaFree(d_seq_starts_);
        if (d_seq_count_) cudaFree(d_seq_count_);
        if (d_descriptors_) cudaFree(d_descriptors_);
        if (d_temp_storage_) cudaFree(d_temp_storage_);
        
        d_data_ = nullptr;
        d_marks_ = nullptr;
        d_scan_ = nullptr;
        d_seq_starts_ = nullptr;
        d_seq_count_ = nullptr;
        d_descriptors_ = nullptr;
        d_temp_storage_ = nullptr;
        allocated_size_ = 0;
    }
    
    // Parse a batch of FASTQ data already on the GPU
    // Returns the number of sequences found
    // newline_offset: cumulative newline count from previous batches mod 4
    uint32_t parse_batch_gpu(size_t size, uint32_t newline_offset = 0) {
        if (size == 0) return 0;
        
        const int block_size = 256;
        const int num_blocks = (size + block_size - 1) / block_size;
        
        // Step 1: Mark newlines
        mark_newlines_kernel<<<num_blocks, block_size, 0, stream_>>>(
            d_data_, d_marks_, size);
        
        // Step 2: Inclusive scan
        CUDA_CHECK(cub::DeviceScan::InclusiveSum(
            d_temp_storage_, temp_storage_bytes_,
            d_marks_, d_scan_, size, stream_));
        
        // Step 3: Reset sequence count
        CUDA_CHECK(cudaMemsetAsync(d_seq_count_, 0, sizeof(uint32_t), stream_));
        
        // Step 4: Find sequence positions
        find_sequence_positions_kernel<<<num_blocks, block_size, 0, stream_>>>(
            d_data_, d_scan_, d_seq_starts_, d_seq_count_, size, newline_offset);
        
        // Get sequence count
        uint32_t num_sequences;
        CUDA_CHECK(cudaMemcpyAsync(&num_sequences, d_seq_count_, sizeof(uint32_t),
                                    cudaMemcpyDeviceToHost, stream_));
        CUDA_CHECK(cudaStreamSynchronize(stream_));
        
        if (num_sequences == 0) return 0;

        // Step 4b: Sort seq_starts so descriptors are in file order
        // (atomicAdd in find_sequence_positions_kernel produces non-deterministic order)
        thrust::device_ptr<uint32_t> starts_ptr(d_seq_starts_);
        thrust::sort(thrust::cuda::par.on(stream_), starts_ptr, starts_ptr + num_sequences);
        
        // Step 5: Compute sequence lengths
        const int seq_blocks = (num_sequences + block_size - 1) / block_size;
        compute_sequence_lengths_kernel<<<seq_blocks, block_size, 0, stream_>>>(
            d_data_, d_seq_starts_, d_descriptors_, num_sequences, size);
        
        CUDA_CHECK(cudaStreamSynchronize(stream_));
        
        return num_sequences;
    }
    
    // Parse from host data
    uint32_t parse_batch(const char* h_data, size_t size, uint32_t newline_offset = 0) {
        allocate_buffers(size);
        
        // Copy data to GPU
        CUDA_CHECK(cudaMemcpyAsync(d_data_, h_data, size,
                                    cudaMemcpyHostToDevice, stream_));
        
        return parse_batch_gpu(size, newline_offset);
    }
    
    // Parse from device data (no copy needed - data already on GPU)
    // The external_d_data pointer must remain valid during parsing
    uint32_t parse_device_data(const char* external_d_data, size_t size, 
                               uint32_t newline_offset = 0, cudaStream_t ext_stream = 0) {
        if (size == 0) return 0;
        
        cudaStream_t use_stream = ext_stream ? ext_stream : stream_;
        
        // Allocate working buffers only (not d_data_ since we use external)
        if (size > allocated_size_) {
            // Free old buffers
            if (d_marks_) cudaFree(d_marks_);
            if (d_scan_) cudaFree(d_scan_);
            if (d_seq_starts_) cudaFree(d_seq_starts_);
            if (d_seq_count_) cudaFree(d_seq_count_);
            if (d_descriptors_) cudaFree(d_descriptors_);
            if (d_temp_storage_) cudaFree(d_temp_storage_);
            
            CUDA_CHECK(cudaMalloc(&d_marks_, size * sizeof(uint32_t)));
            CUDA_CHECK(cudaMalloc(&d_scan_, size * sizeof(uint32_t)));
            
            size_t max_sequences = size / 50;
            CUDA_CHECK(cudaMalloc(&d_seq_starts_, max_sequences * sizeof(uint32_t)));
            CUDA_CHECK(cudaMalloc(&d_seq_count_, sizeof(uint32_t)));
            CUDA_CHECK(cudaMalloc(&d_descriptors_, max_sequences * sizeof(SeqDescriptor)));
            
            CUDA_CHECK(cub::DeviceScan::InclusiveSum(
                nullptr, temp_storage_bytes_,
                d_marks_, d_scan_, size, use_stream));
            CUDA_CHECK(cudaMalloc(&d_temp_storage_, temp_storage_bytes_));
            
            allocated_size_ = size;
        }
        
        const int block_size = 256;
        const int num_blocks = (size + block_size - 1) / block_size;
        
        // Step 1: Mark newlines (use external data)
        mark_newlines_kernel<<<num_blocks, block_size, 0, use_stream>>>(
            external_d_data, d_marks_, size);
        
        // Step 2: Inclusive scan
        CUDA_CHECK(cub::DeviceScan::InclusiveSum(
            d_temp_storage_, temp_storage_bytes_,
            d_marks_, d_scan_, size, use_stream));
        
        // Step 3: Reset sequence count
        CUDA_CHECK(cudaMemsetAsync(d_seq_count_, 0, sizeof(uint32_t), use_stream));
        
        // Step 4: Find sequence positions
        find_sequence_positions_kernel<<<num_blocks, block_size, 0, use_stream>>>(
            external_d_data, d_scan_, d_seq_starts_, d_seq_count_, size, newline_offset);
        
        // Get sequence count
        uint32_t num_sequences;
        CUDA_CHECK(cudaMemcpyAsync(&num_sequences, d_seq_count_, sizeof(uint32_t),
                                    cudaMemcpyDeviceToHost, use_stream));
        CUDA_CHECK(cudaStreamSynchronize(use_stream));
        
        if (num_sequences == 0) return 0;

        // Step 4b: Sort seq_starts so descriptors are in file order
        thrust::device_ptr<uint32_t> starts_ptr(d_seq_starts_);
        thrust::sort(thrust::cuda::par.on(use_stream), starts_ptr, starts_ptr + num_sequences);
        
        // Step 5: Compute sequence lengths
        const int seq_blocks = (num_sequences + block_size - 1) / block_size;
        compute_sequence_lengths_kernel<<<seq_blocks, block_size, 0, use_stream>>>(
            external_d_data, d_seq_starts_, d_descriptors_, num_sequences, size);
        
        CUDA_CHECK(cudaStreamSynchronize(use_stream));
        
        // Store reference to external data for later use
        d_data_ = const_cast<char*>(external_d_data);
        
        return num_sequences;
    }
    
    // Get the total newline count from the last parse (for batch continuation)
    uint32_t get_total_newlines(size_t size) {
        uint32_t total;
        CUDA_CHECK(cudaMemcpy(&total, d_scan_ + size - 1, sizeof(uint32_t),
                              cudaMemcpyDeviceToHost));
        return total;
    }
    
    // Get total newlines using external scan buffer
    uint32_t get_total_newlines_from_scan(size_t size) {
        uint32_t total;
        CUDA_CHECK(cudaMemcpy(&total, d_scan_ + size - 1, sizeof(uint32_t),
                              cudaMemcpyDeviceToHost));
        return total;
    }
    
    // Copy sequence descriptors to host
    void get_descriptors(SeqDescriptor* h_descriptors, uint32_t num_sequences) {
        CUDA_CHECK(cudaMemcpy(h_descriptors, d_descriptors_,
                              num_sequences * sizeof(SeqDescriptor),
                              cudaMemcpyDeviceToHost));
    }
    
    // Copy specific sequences to host buffer
    void extract_sequences(char* h_output, const SeqDescriptor* h_descriptors,
                          uint32_t num_sequences, size_t max_output_size) {
        // This copies the raw data for extracting sequences on CPU
        // For large-scale use, you'd want to do this on GPU
        std::vector<char> h_data(allocated_size_);
        CUDA_CHECK(cudaMemcpy(h_data.data(), d_data_, allocated_size_,
                              cudaMemcpyDeviceToHost));
        
        size_t offset = 0;
        for (uint32_t i = 0; i < num_sequences && offset < max_output_size; i++) {
            size_t len = std::min((size_t)h_descriptors[i].length,
                                  max_output_size - offset - 1);
            memcpy(h_output + offset, h_data.data() + h_descriptors[i].offset, len);
            h_output[offset + len] = '\n';
            offset += len + 1;
        }
    }
    
    // Accessors
    char* device_data() { return d_data_; }
    const char* device_data() const { return d_data_; }
    SeqDescriptor* device_descriptors() { return d_descriptors_; }
    cudaStream_t stream() { return stream_; }
    size_t allocated_size() const { return allocated_size_; }
    
private:
    size_t max_batch_size_;
    char* d_data_;
    uint32_t* d_marks_;
    uint32_t* d_scan_;
    uint32_t* d_seq_starts_;
    uint32_t* d_seq_count_;
    SeqDescriptor* d_descriptors_;
    void* d_temp_storage_;
    size_t temp_storage_bytes_;
    size_t allocated_size_;
    cudaStream_t stream_;
};

