#ifndef GPU_PIPELINE_CUH
#define GPU_PIPELINE_CUH

#include "GPUKernels.cuh"
#include "GPUIndex.cuh"
// #include "Kmer.hpp"
#include "KmerIndex.h"
#include "ProcessReads.h"
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/fill.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/gather.h>
#include <thrust/reduce.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cuda_runtime.h>
// GPU pipeline uses the same FASTQ reading mechanism as the CPU path
// via FastqSequenceReader (defined in ProcessReads.h).
#include <chrono>
#include "BenchmarkStats.h"

// Forward declaration (full definition in GPUReadLoader.cuh)
struct GPUReadLoader;

// Read loading using the same FASTQ machinery as the CPU path.
// For paired-end: reads are stored as [R1_0, R1_1, ..., R1_n, R2_0, R2_1, ..., R2_n]
// matching the layout used by GPUReadLoader. r1_count gives the boundary.
// pair_limit specifies how many read PAIRS to load (not individual reads).
struct HostReadLoader {
  FastqSequenceReader reader;
  std::vector<char> buffer;
  std::vector<uint64_t> read_id_to_offset;
  std::vector<uint32_t> read_id_to_length;
  std::vector<uint64_t> read_id_to_kmer_first;
  std::vector<uint32_t> read_id_to_kmer_count;
  std::vector<char> reads;
  uint64_t kmer_count;
  uint32_t read_count;
  uint32_t r1_count;
  uint32_t pair_limit;
  bool is_paired;

  HostReadLoader(const ProgramOptions& opt,
                 uint32_t pair_limit,
                 uint32_t read_length_hint = 200)
    : reader(opt),
      kmer_count(0),
      read_count(0),
      r1_count(0),
      pair_limit(pair_limit),
      is_paired(!opt.single_end && opt.files.size() >= 2)
  {
    uint32_t max_reads = is_paired ? (pair_limit * 2) : pair_limit;
    read_id_to_offset.resize(max_reads);
    read_id_to_length.resize(max_reads);
    read_id_to_kmer_first.resize(max_reads);
    read_id_to_kmer_count.resize(max_reads);

    buffer.resize(static_cast<size_t>(max_reads) *
                  static_cast<size_t>(read_length_hint) * 2);
    reads.reserve(static_cast<size_t>(max_reads) *
                  static_cast<size_t>(read_length_hint));
  }

  bool load() {
    auto start_time = std::chrono::high_resolution_clock::now();

    reads.clear();
    read_count = 0;
    r1_count = 0;
    kmer_count = 0;

    if (reader.empty()) {
      return false;
    }

    // fetchSequences returns interleaved [R1_0, R2_0, R1_1, R2_1, ...] for paired-end.
    // We collect into separate temporary vectors, then lay out [R1s..., R2s...].

    struct ReadEntry {
      uint32_t length;
      uint64_t temp_offset;  // offset into reads vector at time of insertion
    };
    std::vector<ReadEntry> r1_entries, r2_entries;
    r1_entries.reserve(pair_limit);
    if (is_paired) r2_entries.reserve(pair_limit);

    std::vector<std::pair<const char*, int>> seqs;
    std::vector<std::pair<const char*, int>> names;
    std::vector<std::pair<const char*, int>> quals;
    std::vector<uint32_t> flags;
    std::vector<std::string> umis;
    int readbatch_id = -1;

    uint32_t pairs_loaded = 0;

    while (pairs_loaded < pair_limit) {
      if (reader.empty()) break;

      bool has_more = reader.fetchSequences(
        buffer.data(),
        static_cast<int>(buffer.size()),
        seqs, names, quals, flags, umis, readbatch_id,
        /*full=*/false, /*comments=*/false
      );

      if (seqs.empty()) {
        if (!has_more) break;
        continue;
      }

      if (is_paired) {
        // seqs comes interleaved: [R1, R2, R1, R2, ...]
        for (size_t j = 0; j + 1 < seqs.size(); j += 2) {
          if (pairs_loaded >= pair_limit) break;
          const char* s1 = seqs[j].first;
          int len1 = seqs[j].second;
          const char* s2 = seqs[j + 1].first;
          int len2 = seqs[j + 1].second;

          uint64_t off1 = static_cast<uint64_t>(reads.size());
          reads.insert(reads.end(), s1, s1 + len1);
          r1_entries.push_back({static_cast<uint32_t>(len1), off1});

          uint64_t off2 = static_cast<uint64_t>(reads.size());
          reads.insert(reads.end(), s2, s2 + len2);
          r2_entries.push_back({static_cast<uint32_t>(len2), off2});

          pairs_loaded++;
        }
      } else {
        for (size_t j = 0; j < seqs.size(); j++) {
          if (pairs_loaded >= pair_limit) break;
          const char* s1 = seqs[j].first;
          int len1 = seqs[j].second;

          uint64_t off1 = static_cast<uint64_t>(reads.size());
          reads.insert(reads.end(), s1, s1 + len1);
          r1_entries.push_back({static_cast<uint32_t>(len1), off1});

          pairs_loaded++;
        }
      }

      if (!has_more) break;
    }

    // Now build the metadata arrays in [R1_0..R1_n, R2_0..R2_n] order.
    uint32_t n_r1 = static_cast<uint32_t>(r1_entries.size());
    uint32_t n_r2 = static_cast<uint32_t>(r2_entries.size());
    uint32_t total = n_r1 + n_r2;
    if (total == 0) return false;

    if (total > static_cast<uint32_t>(read_id_to_offset.size())) {
      read_id_to_offset.resize(total);
      read_id_to_length.resize(total);
      read_id_to_kmer_first.resize(total);
      read_id_to_kmer_count.resize(total);
    }

    uint64_t running_kmer = 0;
    uint32_t k = Kmer::k;

    // R1 entries: indices [0, n_r1)
    for (uint32_t i = 0; i < n_r1; i++) {
      read_id_to_offset[i] = r1_entries[i].temp_offset;
      read_id_to_length[i] = r1_entries[i].length;
      uint32_t nk = (r1_entries[i].length >= k) ? (r1_entries[i].length - k + 1) : 0;
      read_id_to_kmer_first[i] = running_kmer;
      read_id_to_kmer_count[i] = nk;
      running_kmer += nk;
    }

    // R2 entries: indices [n_r1, n_r1 + n_r2)
    for (uint32_t i = 0; i < n_r2; i++) {
      uint32_t idx = n_r1 + i;
      read_id_to_offset[idx] = r2_entries[i].temp_offset;
      read_id_to_length[idx] = r2_entries[i].length;
      uint32_t nk = (r2_entries[i].length >= k) ? (r2_entries[i].length - k + 1) : 0;
      read_id_to_kmer_first[idx] = running_kmer;
      read_id_to_kmer_count[idx] = nk;
      running_kmer += nk;
    }

    read_count = total;
    r1_count = n_r1;
    kmer_count = running_kmer;

    // Trim vectors to actual size so downstream code doesn't see stale data
    read_id_to_offset.resize(total);
    read_id_to_length.resize(total);
    read_id_to_kmer_first.resize(total);
    read_id_to_kmer_count.resize(total);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    g_benchmark_stats.io_decompress_ms += duration.count() / 1000.0;

    return !reads.empty();
  }
};

struct DeviceKmerLoader {
  thrust::device_vector<char> reads;
  thrust::device_vector<uint64_t> read_id_to_offset;
  thrust::device_vector<uint32_t> read_id_to_length;
  thrust::device_vector<uint64_t> read_id_to_kmer_first;
  thrust::device_vector<uint32_t> read_id_to_kmer_count;
  thrust::device_vector<uint64_t> kmers;
  uint64_t read_count;
  uint32_t r1_count;  // Number of R1 reads (R2 starts at this index); 0 for single-end

  DeviceKmerLoader() : read_count(0), r1_count(0) {}

  void clear() {
    read_count = 0;
    r1_count = 0;
  }

  void copy_from_host_to_device(const HostReadLoader& h) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    reads.resize(h.reads.size());
    thrust::copy(h.reads.begin(), h.reads.end(), reads.begin());
    
    read_id_to_offset.resize(h.read_id_to_offset.size());
    thrust::copy(h.read_id_to_offset.begin(), h.read_id_to_offset.end(), read_id_to_offset.begin());
    
    read_id_to_length.resize(h.read_id_to_length.size());
    thrust::copy(h.read_id_to_length.begin(), h.read_id_to_length.end(), read_id_to_length.begin());

    read_id_to_kmer_first.resize(h.read_id_to_kmer_first.size());
    thrust::copy(h.read_id_to_kmer_first.begin(), h.read_id_to_kmer_first.end(), read_id_to_kmer_first.begin());
    
    read_id_to_kmer_count.resize(h.read_id_to_kmer_count.size());
    thrust::copy(h.read_id_to_kmer_count.begin(), h.read_id_to_kmer_count.end(), read_id_to_kmer_count.begin());
    
    kmers.resize(h.kmer_count);

    read_count = h.read_count;
    r1_count = h.r1_count;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    g_benchmark_stats.host_h2d_copy_ms += duration.count() / 1000.0;
  }

  void load(const HostReadLoader& h) {
    copy_from_host_to_device(h);
  }

  // Load from GPUReadLoader - data is already on GPU
  // Implementation in GPUProcessReads.cu to avoid CUDA kernel multiple definitions
  void load_from_gpu(GPUReadLoader& gpu);

  void run() {
    Kmer empty;
    empty.set_empty();
    uint64_t empty_kmer_value = to_ullong(empty);
    
    int threads = 1024;
    int blocks = (read_count + threads - 1) / threads;
    
    kmer_kernel<<<blocks, threads>>>(
      reads.data().get(),
      read_id_to_offset.data().get(),
      read_id_to_length.data().get(),
      read_id_to_kmer_count.data().get(),
      read_id_to_kmer_first.data().get(),
      kmers.data().get(),
      empty_kmer_value,
      read_count,
      Kmer::k
    );
  }
};

// Pipeline stages
struct ReadECCollapser {
  thrust::device_vector<int> read_ecs;
  thrust::device_vector<uint64_t> read_ec_offsets;
  thrust::device_vector<int> temp_ecs_;  // Cache for unique ECs between kernels
  uint64_t read_count;
  
  ReadECCollapser() : read_count(0) {}
  
  void collapse(const DeviceKmerLoader& d_loader, 
                const thrust::device_vector<int>& d_ecs) {
    read_count = d_loader.read_count;
    
    if (read_count == 0) {
      read_ec_offsets.resize(1, 0);
      read_ecs.clear();
      return;
    }
    
    thrust::device_vector<uint64_t> output_sizes(read_count);
    
    // Allocate temp buffer for caching unique ECs (reuse across batches)
    size_t temp_size = read_count * MAX_ECS_PER_READ;
    if (temp_ecs_.size() < temp_size) {
      temp_ecs_.resize(temp_size);
    }
    
    int threads = 256;
    int blocks = (read_count + threads - 1) / threads;
    
    cudaEvent_t collapse_start, collapse_stop;
    cudaEventCreate(&collapse_start);
    cudaEventCreate(&collapse_stop);
    cudaEventRecord(collapse_start);
    
    // First kernel: dedup, sort, and cache results in temp_ecs_
    compute_ec_collapse_sizes_kernel<<<blocks, threads>>>(
      d_loader.read_id_to_kmer_first.data().get(),
      d_loader.read_id_to_kmer_count.data().get(),
      d_ecs.data().get(),
      output_sizes.data().get(),
      temp_ecs_.data().get(),
      read_count
    );
    
    read_ec_offsets.resize(read_count + 1);
    auto host_scan_start = std::chrono::high_resolution_clock::now();
    thrust::exclusive_scan(
      output_sizes.begin(),
      output_sizes.end(),
      read_ec_offsets.begin(),
      0ULL
    );
    auto host_scan_end = std::chrono::high_resolution_clock::now();
    auto host_scan_duration = std::chrono::duration_cast<std::chrono::microseconds>(host_scan_end - host_scan_start);
    g_benchmark_stats.host_thrust_ops_ms += host_scan_duration.count() / 1000.0;
    
    uint64_t total_size = 0;
    if (read_count > 0) {
      total_size = read_ec_offsets[read_count - 1] + output_sizes[read_count - 1];
      read_ec_offsets[read_count] = total_size;
    }
    
    if (total_size == 0) {
      read_ecs.resize(1);
    } else {
      read_ecs.resize(total_size);
    }
    
    // Second kernel: just copy from temp_ecs_ to output (no dedup needed)
    collapse_ecs_per_read_kernel<<<blocks, threads>>>(
      temp_ecs_.data().get(),
      read_ec_offsets.data().get(),
      read_ecs.data().get(),
      read_count
    );
    
    cudaEventRecord(collapse_stop);
    cudaEventSynchronize(collapse_stop);
    float collapse_ms = 0;
    cudaEventElapsedTime(&collapse_ms, collapse_start, collapse_stop);
    g_benchmark_stats.gpu_ec_collapse_ms += collapse_ms;
    cudaEventDestroy(collapse_start);
    cudaEventDestroy(collapse_stop);
  }
};

struct ReadTranscriptIntersector {
  thrust::device_vector<int> read_transcripts;
  thrust::device_vector<uint64_t> read_transcript_offsets;
  thrust::device_vector<uint64_t> read_transcript_sizes;
  thrust::device_vector<uint8_t> read_had_mapped_kmers;
  uint64_t read_count;
  
  ReadTranscriptIntersector() : read_count(0) {}

  void intersect_gpu(const ReadECCollapser& collapser, const GPUECMap& gpu_ecmap) {
    read_count = collapser.read_count;
    
    if (read_count == 0) {
      read_transcript_offsets.resize(1, 0);
      read_transcript_sizes.clear();
      read_transcripts.clear();
      read_had_mapped_kmers.clear();
      return;
    }
    
    thrust::device_vector<uint64_t> output_sizes(read_count);
    thrust::device_vector<uint64_t> smallest_ec_indices(read_count);
    
    // Optimization 6: Calculate optimal launch configuration once and reuse for both kernels
    // This reduces kernel launch overhead and ensures consistent grid configuration
    // Using 1024 threads per block is a good default for maximizing occupancy
    const int threads_per_block = 256;
    const int blocks = (read_count + threads_per_block - 1) / threads_per_block;
    
    cudaEvent_t intersect_start, intersect_stop;
    cudaEventCreate(&intersect_start);
    cudaEventCreate(&intersect_stop);
    cudaEventRecord(intersect_start);
    
    compute_intersection_sizes_kernel<<<blocks, threads_per_block>>>(
      collapser.read_ecs.data().get(),
      collapser.read_ec_offsets.data().get(),
      gpu_ecmap.offsets.data().get(),
      output_sizes.data().get(),
      smallest_ec_indices.data().get(),
      read_count,
      static_cast<uint64_t>(gpu_ecmap.num_ecs),
      static_cast<uint64_t>(gpu_ecmap.transcripts.size()),
      static_cast<uint64_t>(collapser.read_ecs.size())
    );
    
    read_transcript_offsets.resize(read_count + 1);
    auto host_scan_start = std::chrono::high_resolution_clock::now();
    thrust::exclusive_scan(
      output_sizes.begin(),
      output_sizes.end(),
      read_transcript_offsets.begin(),
      0ULL
    );
    auto host_scan_end = std::chrono::high_resolution_clock::now();
    auto host_scan_duration = std::chrono::duration_cast<std::chrono::microseconds>(host_scan_end - host_scan_start);
    g_benchmark_stats.host_thrust_ops_ms += host_scan_duration.count() / 1000.0;
    
    uint64_t total_size = 0;
    if (read_count > 0) {
      total_size = read_transcript_offsets[read_count - 1] + output_sizes[read_count - 1];
      read_transcript_offsets[read_count] = total_size;
    }
    
    if (total_size == 0) {
      read_transcripts.resize(1);
    } else {
      read_transcripts.resize(total_size);
    }
    
    // Optimization 6: Reuse the same launch configuration for consistency
    intersect_transcripts_kernel<<<blocks, threads_per_block>>>(
      collapser.read_ecs.data().get(),
      collapser.read_ec_offsets.data().get(),
      gpu_ecmap.transcripts.data().get(),
      gpu_ecmap.offsets.data().get(),
      read_transcript_offsets.data().get(),
      read_transcripts.data().get(),
      output_sizes.data().get(),
      smallest_ec_indices.data().get(),
      read_count,
      static_cast<uint64_t>(gpu_ecmap.num_ecs),
      static_cast<uint64_t>(gpu_ecmap.transcripts.size()),
      static_cast<uint64_t>(collapser.read_ecs.size()),
      static_cast<uint64_t>(read_transcripts.size())
    );
    
    cudaEventRecord(intersect_stop);
    cudaEventSynchronize(intersect_stop);
    float intersect_ms = 0;
    cudaEventElapsedTime(&intersect_ms, intersect_start, intersect_stop);
    g_benchmark_stats.gpu_transcript_intersection_ms += intersect_ms;
    cudaEventDestroy(intersect_start);
    cudaEventDestroy(intersect_stop);
    
    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess) {
      std::cerr << "CUDA launch error for intersect_transcripts_kernel: " 
                << cudaGetErrorString(launch_err) << std::endl;
      exit(1);
    }
    
    read_transcript_sizes = output_sizes;

    read_had_mapped_kmers.resize(read_count);
    const uint64_t* ec_off = collapser.read_ec_offsets.data().get();
    thrust::transform(
      thrust::counting_iterator<uint64_t>(0),
      thrust::counting_iterator<uint64_t>(read_count),
      read_had_mapped_kmers.begin(),
      [ec_off] __device__ (uint64_t i) -> uint8_t {
        return (ec_off[i + 1] > ec_off[i]) ? 1 : 0;
      }
    );
  }
};

// Paired-end intersector: intersects R1[i] and R2[i] transcript sets
struct PairIntersector {
  thrust::device_vector<int> pair_transcripts;
  thrust::device_vector<uint64_t> pair_transcript_offsets;
  thrust::device_vector<uint64_t> pair_transcript_sizes;
  uint32_t pair_count;

  PairIntersector() : pair_count(0) {}

  void intersect_pairs(const ReadTranscriptIntersector& intersector, uint32_t r1_count) {
    pair_count = r1_count;

    if (pair_count == 0) {
      pair_transcript_offsets.resize(1, 0);
      pair_transcript_sizes.clear();
      pair_transcripts.clear();
      return;
    }

    cudaEvent_t pair_start, pair_stop;
    cudaEventCreate(&pair_start);
    cudaEventCreate(&pair_stop);
    cudaEventRecord(pair_start);

    const int threads = 256;
    const int blocks = (pair_count + threads - 1) / threads;

    // Pass 1: compute intersection sizes
    thrust::device_vector<uint64_t> sizes(pair_count);
    compute_pair_intersection_sizes_kernel<<<blocks, threads>>>(
      intersector.read_transcripts.data().get(),
      intersector.read_transcript_offsets.data().get(),
      intersector.read_transcript_sizes.data().get(),
      intersector.read_had_mapped_kmers.data().get(),
      sizes.data().get(),
      r1_count
    );

    // Prefix sum to get offsets
    pair_transcript_offsets.resize(pair_count + 1);
    thrust::exclusive_scan(sizes.begin(), sizes.end(),
                           pair_transcript_offsets.begin(), 0ULL);

    uint64_t total_size = 0;
    if (pair_count > 0) {
      total_size = pair_transcript_offsets[pair_count - 1] + sizes[pair_count - 1];
      pair_transcript_offsets[pair_count] = total_size;
    }

    if (total_size == 0) {
      pair_transcripts.resize(1);
    } else {
      pair_transcripts.resize(total_size);
    }

    // Pass 2: write intersections
    pair_transcript_sizes.resize(pair_count);
    intersect_pairs_kernel<<<blocks, threads>>>(
      intersector.read_transcripts.data().get(),
      intersector.read_transcript_offsets.data().get(),
      intersector.read_transcript_sizes.data().get(),
      intersector.read_had_mapped_kmers.data().get(),
      pair_transcript_offsets.data().get(),
      pair_transcripts.data().get(),
      pair_transcript_sizes.data().get(),
      r1_count
    );

    cudaEventRecord(pair_stop);
    cudaEventSynchronize(pair_stop);
    float pair_ms = 0;
    cudaEventElapsedTime(&pair_ms, pair_start, pair_stop);
    cudaEventDestroy(pair_start);
    cudaEventDestroy(pair_stop);
  }
};

struct ReadECLookup {
  thrust::device_vector<int> read_ecs_final;
  thrust::device_vector<uint64_t> read_hashes;
  uint64_t read_count;
  
  ReadECLookup() : read_count(0) {}

  // Helper: compute hashes, find via dynamic_map, override special cases, verify against GPUECMap
  void lookup_impl(const int* d_transcripts, const uint64_t* d_offsets,
                   const uint64_t* d_sizes, uint64_t max_transcripts,
                   uint64_t count, GPUECMapInv& gpu_ecmapinv, const GPUECMap& gpu_ecmap) {
    read_count = count;

    if (read_count == 0) {
      read_ecs_final.clear();
      read_hashes.clear();
      return;
    }

    read_ecs_final.resize(read_count);
    read_hashes.resize(read_count);

    cudaEvent_t lookup_start, lookup_stop;
    cudaEventCreate(&lookup_start);
    cudaEventCreate(&lookup_stop);
    cudaEventRecord(lookup_start);

    // Step 1: Compute hashes for each transcript set
    const int* d_tx = d_transcripts;
    const uint64_t* d_off = d_offsets;
    const uint64_t* d_sz = d_sizes;
    uint64_t d_max = max_transcripts;

    thrust::counting_iterator<uint64_t> indices(0);
    thrust::transform(
      thrust::device,
      indices, indices + read_count,
      read_hashes.begin(),
      [=] __device__ (uint64_t id) -> uint64_t {
        uint64_t tx_start = d_off[id];
        uint64_t tx_size = d_sz[id];
        if (tx_size == 0 || tx_start + tx_size > d_max) {
          return UINT64_MAX;
        }
        return hash_sorted_vector_well_defined_device(d_tx + tx_start, tx_size);
      }
    );

    // Step 2: Bulk find via dynamic_map -> candidate ec_ids
    gpu_ecmapinv.hash_map.find(read_hashes.begin(), read_hashes.end(),
                               read_ecs_final.begin());

    // Step 3: Override empty transcript sets -> ec = -1
    thrust::transform(
      thrust::device,
      indices, indices + read_count,
      read_ecs_final.begin(),
      read_ecs_final.begin(),
      [=] __device__ (uint64_t id, int candidate_ec) -> int {
        uint64_t tx_start = d_off[id];
        uint64_t tx_size = d_sz[id];
        if (tx_size == 0 || tx_start + tx_size > d_max) {
          return -1;
        }
        return candidate_ec;
      }
    );

    // Step 4: Verify against GPUECMap to catch hash collisions
    // All non-empty transcript sets (including singletons) go through verification
    const int* d_ecmap_tx = gpu_ecmap.transcripts.data().get();
    const uint64_t* d_ecmap_off = gpu_ecmap.offsets.data().get();
    size_t d_ecmap_num_ecs = gpu_ecmap.num_ecs;
    size_t d_ecmap_max_tx = gpu_ecmap.transcripts.size();

    thrust::transform(
      thrust::device,
      indices, indices + read_count,
      read_ecs_final.begin(),
      read_ecs_final.begin(),
      [=] __device__ (uint64_t id, int ec_id) -> int {
        if (ec_id < 0) return ec_id;

        uint64_t tx_start = d_off[id];
        uint64_t tx_size = d_sz[id];

        if (static_cast<size_t>(ec_id) >= d_ecmap_num_ecs) return -1;

        uint64_t stored_start = d_ecmap_off[ec_id];
        uint64_t stored_end = d_ecmap_off[ec_id + 1];
        uint64_t stored_size = stored_end - stored_start;

        if (stored_start + stored_size > d_ecmap_max_tx) return -1;

        bool match = verify_transcript_lists_equal(
          d_tx + tx_start, tx_size,
          d_ecmap_tx + stored_start, stored_size
        );
        return match ? ec_id : -1;
      }
    );

    cudaEventRecord(lookup_stop);
    cudaEventSynchronize(lookup_stop);
    float lookup_ms = 0;
    cudaEventElapsedTime(&lookup_ms, lookup_start, lookup_stop);
    g_benchmark_stats.gpu_ec_lookup_ms += lookup_ms;
    cudaEventDestroy(lookup_start);
    cudaEventDestroy(lookup_stop);
  }

  void lookup(const ReadTranscriptIntersector& intersector,
              GPUECMapInv& gpu_ecmapinv, const GPUECMap& gpu_ecmap) {
    lookup_impl(
      intersector.read_transcripts.data().get(),
      intersector.read_transcript_offsets.data().get(),
      intersector.read_transcript_sizes.data().get(),
      intersector.read_transcripts.size(),
      intersector.read_count,
      gpu_ecmapinv, gpu_ecmap
    );
  }

  void lookup_pairs(const PairIntersector& pairs,
                    GPUECMapInv& gpu_ecmapinv, const GPUECMap& gpu_ecmap) {
    lookup_impl(
      pairs.pair_transcripts.data().get(),
      pairs.pair_transcript_offsets.data().get(),
      pairs.pair_transcript_sizes.data().get(),
      pairs.pair_transcripts.size(),
      pairs.pair_count,
      gpu_ecmapinv, gpu_ecmap
    );
  }
};

struct ECCounter {
  thrust::device_vector<int> ec_counts;
  size_t num_ecs;
  
  ECCounter(size_t num_ecs) : num_ecs(num_ecs) {
    ec_counts.resize(num_ecs, 0);
  }

  void ensure_capacity(size_t new_num_ecs) {
    if (new_num_ecs > num_ecs) {
      ec_counts.resize(new_num_ecs, 0);
      num_ecs = new_num_ecs;
    }
  }
  
  void count_batch(const ReadECLookup& ec_lookup) {
    if (ec_lookup.read_count == 0) {
      return;
    }
    
    cudaEvent_t count_start, count_stop;
    cudaEventCreate(&count_start);
    cudaEventCreate(&count_stop);
    cudaEventRecord(count_start);
    
    thrust::device_vector<int> valid_ecs(ec_lookup.read_ecs_final.size());
    auto result_end = thrust::remove_copy_if(
      ec_lookup.read_ecs_final.begin(),
      ec_lookup.read_ecs_final.end(),
      valid_ecs.begin(),
      [] __device__ (int ec) {
        return ec == -1;
      }
    );
    
    size_t num_valid = result_end - valid_ecs.begin();
    if (num_valid == 0) {
      cudaEventRecord(count_stop);
      cudaEventSynchronize(count_stop);
      float count_ms = 0;
      cudaEventElapsedTime(&count_ms, count_start, count_stop);
      g_benchmark_stats.gpu_ec_counting_ms += count_ms;
      cudaEventDestroy(count_start);
      cudaEventDestroy(count_stop);
      return;
    }
    valid_ecs.resize(num_valid);
    
    thrust::sort(valid_ecs.begin(), valid_ecs.end());
    
    thrust::device_vector<int> unique_ecs(num_valid);
    thrust::device_vector<int> counts(num_valid);
    
    auto reduce_end = thrust::reduce_by_key(
      valid_ecs.begin(),
      valid_ecs.end(),
      thrust::make_constant_iterator(1),
      unique_ecs.begin(),
      counts.begin()
    );
    
    size_t num_unique = reduce_end.first - unique_ecs.begin();
    unique_ecs.resize(num_unique);
    counts.resize(num_unique);
    
    int* d_ec_counts = ec_counts.data().get();
    size_t d_num_ecs = num_ecs;
    
    thrust::for_each(
      thrust::device,
      thrust::make_zip_iterator(thrust::make_tuple(unique_ecs.begin(), counts.begin())),
      thrust::make_zip_iterator(thrust::make_tuple(unique_ecs.end(), counts.end())),
      [=] __device__ (thrust::tuple<int, int> t) {
        int ec_id = thrust::get<0>(t);
        int count = thrust::get<1>(t);
        if (ec_id >= 0 && static_cast<size_t>(ec_id) < d_num_ecs) {
          atomicAdd(&d_ec_counts[ec_id], count);
        }
      }
    );
    
    cudaEventRecord(count_stop);
    cudaEventSynchronize(count_stop);
    float count_ms = 0;
    cudaEventElapsedTime(&count_ms, count_start, count_stop);
    g_benchmark_stats.gpu_ec_counting_ms += count_ms;
    cudaEventDestroy(count_start);
    cudaEventDestroy(count_stop);
  }
  
  void write_counts(const std::string& filename) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<int> h_counts(num_ecs);
    thrust::copy(ec_counts.begin(), ec_counts.end(), h_counts.begin());
    
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
      std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
      exit(1);
    }
    
    for (size_t ec_id = 0; ec_id < num_ecs; ++ec_id) {
      outfile << ec_id << "\t" << h_counts[ec_id] << "\n";
    }
    
    outfile.close();
    std::cout << "  Wrote EC counts to: " << filename << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    g_benchmark_stats.host_file_write_ms += duration.count() / 1000.0;
  }
};

// GPU-accelerated handler for novel equivalence classes from paired-end intersections.
// Uses cuco::dynamic_map for O(1) GPU-side dedup and lookup.
struct NewECHandler {
  int next_ec_id;
  int total_new_ec_count;
  int total_new_read_count;

  NewECHandler(int initial_num_ecs)
    : next_ec_id(initial_num_ecs), total_new_ec_count(0), total_new_read_count(0) {}

  // Generic implementation: works on either PairIntersector or ReadTranscriptIntersector.
  // Caller passes the device pointers/sizes for the per-read transcript set.
  void handle_batch_impl(uint64_t source_count,
                         const int* tx_data, uint64_t tx_data_size,
                         const uint64_t* tx_offsets,
                         const uint64_t* tx_sizes,
                         ReadECLookup& ec_lookup,
                         GPUECMapInv& gpu_ecmapinv, GPUECMap& gpu_ecmap,
                         ECCounter& ec_counter) {
    if (ec_lookup.read_count == 0 || source_count == 0) return;

    uint64_t n = ec_lookup.read_count;

    // Step 1: Filter novel indices (ec == -1 and tx_size > 0)
    const int* d_ecs = ec_lookup.read_ecs_final.data().get();
    const uint64_t* d_sizes = tx_sizes;

    thrust::device_vector<uint64_t> novel_indices(n);
    auto novel_end = thrust::copy_if(
      thrust::device,
      thrust::counting_iterator<uint64_t>(0),
      thrust::counting_iterator<uint64_t>(n),
      novel_indices.begin(),
      [=] __device__ (uint64_t i) {
        return d_ecs[i] == -1 && d_sizes[i] > 0;
      }
    );
    size_t num_novel = novel_end - novel_indices.begin();
    if (num_novel == 0) return;
    novel_indices.resize(num_novel);
    total_new_read_count += static_cast<int>(num_novel);

    // Step 2: Gather hashes at novel indices
    const uint64_t* d_read_hashes = ec_lookup.read_hashes.data().get();
    thrust::device_vector<uint64_t> novel_hashes(num_novel);
    thrust::gather(thrust::device,
                   novel_indices.begin(), novel_indices.end(),
                   ec_lookup.read_hashes.begin(),
                   novel_hashes.begin());

    // Step 3: Sort + unique to get deduplicated hashes
    thrust::device_vector<uint64_t> unique_hashes(novel_hashes);
    thrust::sort(unique_hashes.begin(), unique_hashes.end());
    auto unique_end = thrust::unique(unique_hashes.begin(), unique_hashes.end());
    size_t num_unique = unique_end - unique_hashes.begin();
    unique_hashes.resize(num_unique);

    // Step 4: Check which unique hashes are already in the map (collision check)
    thrust::device_vector<int> found_ids(num_unique);
    gpu_ecmapinv.hash_map.find(unique_hashes.begin(), unique_hashes.end(), found_ids.begin());

    // Truly new hashes: those not found in the map (found_id == -1)
    thrust::device_vector<uint64_t> truly_new_hashes(num_unique);
    auto truly_new_end = thrust::copy_if(
      thrust::device,
      unique_hashes.begin(), unique_hashes.end(),
      found_ids.begin(),
      truly_new_hashes.begin(),
      [] __device__ (int found_id) {
        return found_id == -1;
      }
    );
    size_t num_truly_new = truly_new_end - truly_new_hashes.begin();
    truly_new_hashes.resize(num_truly_new);

    if (num_truly_new > 0) {
      // Step 5: For each truly new hash, find a representative pair to extract transcript data
      // We need to copy transcript data to host for GPUECMap::append_ecs
      // For each unique hash, find first novel pair with that hash

      // Copy truly_new_hashes to host
      std::vector<uint64_t> h_new_hashes(num_truly_new);
      thrust::copy(truly_new_hashes.begin(), truly_new_hashes.end(), h_new_hashes.begin());

      // Copy novel_indices and novel_hashes to host for matching
      std::vector<uint64_t> h_novel_indices(num_novel);
      thrust::copy(novel_indices.begin(), novel_indices.end(), h_novel_indices.begin());
      std::vector<uint64_t> h_novel_hashes(num_novel);
      thrust::copy(novel_hashes.begin(), novel_hashes.end(), h_novel_hashes.begin());

      // Copy per-read transcript data to host
      std::vector<uint64_t> h_offsets(source_count + 1);
      thrust::copy(thrust::device_pointer_cast(tx_offsets),
                   thrust::device_pointer_cast(tx_offsets) + source_count + 1,
                   h_offsets.begin());
      std::vector<uint64_t> h_sizes(source_count);
      thrust::copy(thrust::device_pointer_cast(tx_sizes),
                   thrust::device_pointer_cast(tx_sizes) + source_count,
                   h_sizes.begin());

      uint64_t total_tx = h_offsets[source_count];
      std::vector<int> h_transcripts(total_tx);
      if (total_tx > 0) {
        thrust::copy(thrust::device_pointer_cast(tx_data),
                     thrust::device_pointer_cast(tx_data) + total_tx,
                     h_transcripts.begin());
      }

      // Build hash -> first novel pair index map on host
      std::unordered_map<uint64_t, uint64_t> hash_to_first_pair;
      for (size_t i = 0; i < num_novel; ++i) {
        uint64_t h = h_novel_hashes[i];
        if (hash_to_first_pair.find(h) == hash_to_first_pair.end()) {
          hash_to_first_pair[h] = h_novel_indices[i];
        }
      }

      // Build new EC transcript data and insert pairs
      std::vector<cuco::pair<uint64_t, int>> h_insert_pairs;
      h_insert_pairs.reserve(num_truly_new);

      std::vector<int> all_new_tx;
      std::vector<uint64_t> new_offsets;
      new_offsets.reserve(num_truly_new + 1);
      uint64_t tx_offset = 0;

      for (size_t i = 0; i < num_truly_new; ++i) {
        uint64_t hash = h_new_hashes[i];
        int new_id = next_ec_id++;
        h_insert_pairs.push_back({hash, new_id});

        auto it = hash_to_first_pair.find(hash);
        if (it != hash_to_first_pair.end()) {
          uint64_t pair_idx = it->second;
          uint64_t start = h_offsets[pair_idx];
          uint64_t sz = h_sizes[pair_idx];
          new_offsets.push_back(tx_offset);
          all_new_tx.insert(all_new_tx.end(),
                            h_transcripts.begin() + start,
                            h_transcripts.begin() + start + sz);
          tx_offset += sz;
        } else {
          new_offsets.push_back(tx_offset);
        }
      }
      new_offsets.push_back(tx_offset);

      // Step 6: Insert into dynamic_map
      thrust::device_vector<cuco::pair<uint64_t, int>> d_insert_pairs(
        h_insert_pairs.begin(), h_insert_pairs.end());
      gpu_ecmapinv.insert_new_ecs(d_insert_pairs);

      // Step 7: Grow GPUECMap
      gpu_ecmap.append_ecs(all_new_tx, new_offsets, static_cast<int>(num_truly_new));

      // Step 8: Grow ECCounter
      ec_counter.ensure_capacity(next_ec_id);

      total_new_ec_count += static_cast<int>(num_truly_new);
    }

    // Step 9: Re-find novel pair hashes to get assigned ec_ids
    thrust::device_vector<int> novel_ec_ids(num_novel);
    gpu_ecmapinv.hash_map.find(novel_hashes.begin(), novel_hashes.end(), novel_ec_ids.begin());

    // Step 10: Verify re-found results against updated GPUECMap
    const int* d_pair_tx = tx_data;
    const uint64_t* d_pair_off = tx_offsets;
    const uint64_t* d_pair_sz = tx_sizes;
    uint64_t d_pair_max_tx = tx_data_size;
    const int* d_ecmap_tx = gpu_ecmap.transcripts.data().get();
    const uint64_t* d_ecmap_off = gpu_ecmap.offsets.data().get();
    size_t d_ecmap_num_ecs = gpu_ecmap.num_ecs;
    size_t d_ecmap_max_tx = gpu_ecmap.transcripts.size();
    const uint64_t* d_novel_idx = novel_indices.data().get();

    thrust::transform(
      thrust::device,
      thrust::counting_iterator<uint64_t>(0),
      thrust::counting_iterator<uint64_t>(num_novel),
      novel_ec_ids.begin(),
      novel_ec_ids.begin(),
      [=] __device__ (uint64_t i, int ec_id) -> int {
        if (ec_id < 0) return -1;

        uint64_t pair_idx = d_novel_idx[i];
        uint64_t tx_start = d_pair_off[pair_idx];
        uint64_t tx_size = d_pair_sz[pair_idx];

        if (tx_size <= 1 || tx_start + tx_size > d_pair_max_tx) return ec_id;
        if (static_cast<size_t>(ec_id) >= d_ecmap_num_ecs) return -1;

        uint64_t stored_start = d_ecmap_off[ec_id];
        uint64_t stored_end = d_ecmap_off[ec_id + 1];
        uint64_t stored_size = stored_end - stored_start;

        if (stored_start + stored_size > d_ecmap_max_tx) return -1;

        bool match = verify_transcript_lists_equal(
          d_pair_tx + tx_start, tx_size,
          d_ecmap_tx + stored_start, stored_size
        );
        return match ? ec_id : -1;
      }
    );

    // Step 11: Scatter verified ec_ids back into read_ecs_final
    int* d_read_ecs = ec_lookup.read_ecs_final.data().get();
    const int* d_novel_ecs = novel_ec_ids.data().get();

    thrust::for_each(
      thrust::device,
      thrust::counting_iterator<uint64_t>(0),
      thrust::counting_iterator<uint64_t>(num_novel),
      [=] __device__ (uint64_t i) {
        uint64_t pair_idx = d_novel_idx[i];
        d_read_ecs[pair_idx] = d_novel_ecs[i];
      }
    );
    cudaStreamSynchronize(0);
  }

  // Paired-end wrapper
  void handle_batch(const PairIntersector& pairs, ReadECLookup& ec_lookup,
                    GPUECMapInv& gpu_ecmapinv, GPUECMap& gpu_ecmap,
                    ECCounter& ec_counter) {
    handle_batch_impl(
      pairs.pair_count,
      pairs.pair_transcripts.data().get(),
      pairs.pair_transcripts.size(),
      pairs.pair_transcript_offsets.data().get(),
      pairs.pair_transcript_sizes.data().get(),
      ec_lookup, gpu_ecmapinv, gpu_ecmap, ec_counter);
  }

  // Single-end wrapper: novel ECs are minted from the raw per-read intersector output.
  void handle_batch(const ReadTranscriptIntersector& intersector, ReadECLookup& ec_lookup,
                    GPUECMapInv& gpu_ecmapinv, GPUECMap& gpu_ecmap,
                    ECCounter& ec_counter) {
    handle_batch_impl(
      intersector.read_count,
      intersector.read_transcripts.data().get(),
      intersector.read_transcripts.size(),
      intersector.read_transcript_offsets.data().get(),
      intersector.read_transcript_sizes.data().get(),
      ec_lookup, gpu_ecmapinv, gpu_ecmap, ec_counter);
  }

  int total_new_ecs() const { return total_new_ec_count; }
  int total_new_reads() const { return total_new_read_count; }
};

#endif // GPU_PIPELINE_CUH