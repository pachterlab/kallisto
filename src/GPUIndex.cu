#include "GPUIndex.cuh"
#include "GPUPipeline.cuh"
#include "GPUKernels.cuh"
// #include "Kmer.hpp"
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <cuda_runtime.h>
#include <chrono>

cuco::static_map<uint64_t, int> build_kmer_to_ec_map(const KmerIndex& index) {
  auto setup_start = std::chrono::high_resolution_clock::now();
  
  std::vector<cuco::pair<uint64_t, int>> host_kmers;

  size_t kmer_count = 0;
  // Count the number of kmers in the index
  for (const auto &it : index.dbg) {
    size_t contig_length = it.size - index.k + 1;
    kmer_count += contig_length;
  }

  host_kmers.reserve(kmer_count);

  std::vector<SparseVector<uint32_t>> vals;

  for (const const_UnitigMap<Node>& contig : index.dbg) {
    auto n = contig.getData();
    std::string seq = contig.referenceUnitigToString();

    n->ec.get_vals(vals);
    int j = 0;
    size_t contigpos = 0;
    while (contigpos < contig.len) {
      auto mc = n->ec.get_block_at(contigpos);
      const auto& val = vals[j];
      const Roaring& trs = val.getIndices();
      auto blocklen = mc.second - mc.first + index.k - 1;
      auto blockseq = seq.substr(mc.first, blocklen);
      auto ec = index.ecmapinv.find(trs);
      if (ec != index.ecmapinv.end()) {
        KmerIterator kit(blockseq.c_str()), kit_end;
        for (auto it = kit; it != kit_end; ++it) {
          host_kmers.push_back({to_ullong(it->first.rep()), ec->second});
        }
      }
      
      contigpos = mc.second;
      ++j;
    }


    

  }
  // for (const auto& it : index.kmap) {
  //   int ec = index.dbGraph.ecs[it.second.contig];
  //   host_kmers.push_back({to_ullong(it.first), ec});
  //   host_kmers.push_back({to_ullong(it.first.twin()), ec});
  // }
  
  cudaEvent_t setup_start_gpu, setup_stop_gpu;
  cudaEventCreate(&setup_start_gpu);
  cudaEventCreate(&setup_stop_gpu);
  cudaEventRecord(setup_start_gpu);
  
  thrust::device_vector<cuco::pair<uint64_t, int>> device_kmers(host_kmers.begin(),
                                                                host_kmers.end());
  cudaStreamSynchronize(0);
  
  size_t capacity = device_kmers.size() * 2.0;

  Kmer empty;
  empty.set_deleted();
  uint64_t empty_key_sentinel = to_ullong(empty);
  int empty_value_sentinel = -1;

  auto cuco_map = cuco::static_map<uint64_t, int>{
    capacity,
    cuco::empty_key{empty_key_sentinel},
    cuco::empty_value{empty_value_sentinel}
  };

  cuco_map.insert(device_kmers.begin(), device_kmers.end());
  cudaEventRecord(setup_stop_gpu);
  cudaEventSynchronize(setup_stop_gpu);
  
  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_kmer_to_ec_map_ms += setup_duration.count() / 1000.0;
  
  cudaEventDestroy(setup_start_gpu);
  cudaEventDestroy(setup_stop_gpu);
  
  return cuco_map;
}

uint64_t hash_sorted_vector_well_defined_cpu(const Roaring& r) {
  uint64_t hash = 0xcbf29ce484222325ULL;  // FNV-1a offset basis
  for (uint32_t x : r) {
    hash ^= static_cast<uint64_t>(x);
    hash *= 0x100000001b3ULL;  // FNV-1a prime
  }
  // Avoid sentinel value collision
  if (hash == UINT64_MAX) {
    hash = UINT64_MAX - 2;
  }
  return hash;
}

GPUECMap::GPUECMap(const KmerIndex& index) {
  auto setup_start = std::chrono::high_resolution_clock::now();
  
  // Reconstruct forward map (EC id -> transcript list) from ecmapinv
  // ecmapinv maps Roaring bitmap of transcript IDs -> EC id
  int32_t max_ec = -1;
  for (const auto& entry : index.ecmapinv) {
    if (entry.second > max_ec) {
      max_ec = entry.second;
    }
  }
  
  size_t num_ec = (max_ec >= 0) ? static_cast<size_t>(max_ec + 1) : 0;
  
  // Build temporary map: ec_id -> pointer to Roaring transcript set
  std::vector<const Roaring*> ec_to_transcripts(num_ec, nullptr);
  for (const auto& entry : index.ecmapinv) {
    ec_to_transcripts[entry.second] = &entry.first;
  }
  
  // Flatten into offsets + transcripts arrays (same layout as old ecmap)
  std::vector<int> h_transcripts;
  std::vector<uint64_t> h_offsets(num_ec + 1, 0);
  
  uint64_t offset = 0;
  for (size_t ec = 0; ec < num_ec; ++ec) {
    h_offsets[ec] = offset;
    if (ec_to_transcripts[ec] != nullptr) {
      for (uint32_t tr : *ec_to_transcripts[ec]) {
        h_transcripts.push_back(static_cast<int>(tr));
      }
      offset += ec_to_transcripts[ec]->cardinality();
    }
  }
  h_offsets[num_ec] = offset;
  num_ecs = num_ec;
  
  transcripts.assign(h_transcripts.begin(), h_transcripts.end());
  offsets.assign(h_offsets.begin(), h_offsets.end());
  cudaStreamSynchronize(0);
  
  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_gpu_ecmap_ms += setup_duration.count() / 1000.0;
}

GPUECMapInv::GPUECMapInv(const KmerIndex& index)
  : hash_map(std::max<size_t>(index.ecmapinv.size() * 2, 128),
             cuco::empty_key<uint64_t>{UINT64_MAX},
             cuco::empty_value<int>{-1}),
    num_ecs(0) {
  auto setup_start = std::chrono::high_resolution_clock::now();

  std::vector<cuco::pair<uint64_t, int>> host_pairs;
  host_pairs.reserve(index.ecmapinv.size());

  uint64_t duplicate_hashes = 0;
  std::unordered_set<uint64_t> seen_hashes;

  for (const auto& entry : index.ecmapinv) {
    uint64_t hash = hash_sorted_vector_well_defined_cpu(entry.first);
    if (!seen_hashes.insert(hash).second) {
      duplicate_hashes++;
    }
    host_pairs.push_back({hash, entry.second});
  }
  num_ecs = index.ecmapinv.size();

  if (duplicate_hashes > 0) {
    std::cerr << "  WARNING: " << duplicate_hashes << " hash collisions detected in index!" << std::endl;
    std::cerr << "  NOTE: Verification against GPUECMap will handle collisions correctly." << std::endl;
  }

  thrust::device_vector<cuco::pair<uint64_t, int>> device_pairs(host_pairs.begin(), host_pairs.end());
  hash_map.insert(device_pairs.begin(), device_pairs.end());
  cudaStreamSynchronize(0);

  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_gpu_ecmapinv_ms += setup_duration.count() / 1000.0;
}

cuco::static_map<uint64_t, int> build_kmer_to_ec_map_from_contigs(const GPUIndex& index) {
  auto setup_start = std::chrono::high_resolution_clock::now();

  // Build flat sequence and block arrays
  std::vector<char> flat_seq;
  std::vector<uint64_t> block_seq_offsets;
  std::vector<uint32_t> block_num_kmers;
  std::vector<int> block_ec_ids;

  uint64_t seq_offset = 0;
  uint64_t output_offset = 0;

  for (size_t c = 0; c < index.num_contigs; ++c) {
    const std::string& seq = index.contig_sequences[c];
    const auto& blocks = index.contig_ec_blocks[c];
    uint64_t contig_start = seq_offset;
    flat_seq.insert(flat_seq.end(), seq.begin(), seq.end());
    seq_offset += seq.size();

    for (const auto& blk : blocks) {
      uint32_t num_kmers = blk.end - blk.start;
      if (num_kmers == 0) continue;
      uint64_t block_seq_start = contig_start + blk.start;
      block_seq_offsets.push_back(block_seq_start);
      block_num_kmers.push_back(num_kmers);
      block_ec_ids.push_back(blk.ec_id);
      output_offset += num_kmers;
    }
  }

  uint64_t num_blocks = block_seq_offsets.size();
  if (num_blocks == 0 || output_offset == 0) {
    Kmer empty;
    empty.set_deleted();
    return cuco::static_map<uint64_t, int>{
      1, cuco::empty_key{to_ullong(empty)}, cuco::empty_value{-1}};
  }

  thrust::device_vector<char> d_seq(flat_seq.begin(), flat_seq.end());
  thrust::device_vector<uint64_t> d_block_seq_offsets(block_seq_offsets.begin(), block_seq_offsets.end());
  thrust::device_vector<uint32_t> d_block_num_kmers(block_num_kmers.begin(), block_num_kmers.end());
  thrust::device_vector<int> d_block_ec_ids(block_ec_ids.begin(), block_ec_ids.end());

  std::vector<uint64_t> block_output_offsets(num_blocks + 1);
  uint64_t off = 0;
  for (size_t i = 0; i < num_blocks; ++i) {
    block_output_offsets[i] = off;
    off += block_num_kmers[i];
  }
  block_output_offsets[num_blocks] = off;
  thrust::device_vector<uint64_t> d_block_output_offsets(block_output_offsets.begin(), block_output_offsets.end());

  thrust::device_vector<uint64_t> d_out_kmers(output_offset);
  thrust::device_vector<int> d_out_ecs(output_offset);

  Kmer empty;
  empty.set_deleted();
  uint64_t empty_kmer = to_ullong(empty);

  contig_kmer_kernel<<<(num_blocks + 255) / 256, 256>>>(
      thrust::raw_pointer_cast(d_seq.data()),
      thrust::raw_pointer_cast(d_block_seq_offsets.data()),
      thrust::raw_pointer_cast(d_block_num_kmers.data()),
      thrust::raw_pointer_cast(d_block_ec_ids.data()),
      thrust::raw_pointer_cast(d_block_output_offsets.data()),
      thrust::raw_pointer_cast(d_out_kmers.data()),
      thrust::raw_pointer_cast(d_out_ecs.data()),
      empty_kmer, -1,
      static_cast<uint32_t>(index.k),
      num_blocks);
  cudaStreamSynchronize(0);

  thrust::device_vector<cuco::pair<uint64_t, int>> device_pairs(output_offset);
  auto zip_in = thrust::make_zip_iterator(
      thrust::make_tuple(d_out_kmers.begin(), d_out_ecs.begin()));
  thrust::transform(
      zip_in, zip_in + output_offset,
      device_pairs.begin(),
      [] __device__ (thrust::tuple<uint64_t, int> t) {
        return cuco::pair<uint64_t, int>{thrust::get<0>(t), thrust::get<1>(t)};
      });

  auto remove_end = thrust::remove_if(
      device_pairs.begin(),
      device_pairs.end(),
      [=] __device__ (cuco::pair<uint64_t, int> p) {
        return p.second == -1 || p.first == empty_kmer;
      });
  size_t valid_count = thrust::distance(device_pairs.begin(), remove_end);
  device_pairs.resize(valid_count);

  size_t capacity = device_pairs.size() * 2;
  auto cuco_map = cuco::static_map<uint64_t, int>{
    capacity,
    cuco::empty_key{empty_kmer},
    cuco::empty_value{-1}
  };
  cuco_map.insert(device_pairs.begin(), device_pairs.end());

  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_kmer_to_ec_map_ms += setup_duration.count() / 1000.0;

  return cuco_map;
}

GPUECMap::GPUECMap(const GPUIndex& index) {
  auto setup_start = std::chrono::high_resolution_clock::now();

  int32_t max_ec = -1;
  for (const auto& entry : index.ecmapinv) {
    if (entry.second > max_ec) max_ec = entry.second;
  }
  size_t num_ec = (max_ec >= 0) ? static_cast<size_t>(max_ec + 1) : 0;

  std::vector<const Roaring*> ec_to_transcripts(num_ec, nullptr);
  for (const auto& entry : index.ecmapinv) {
    ec_to_transcripts[entry.second] = &entry.first;
  }

  std::vector<int> h_transcripts;
  std::vector<uint64_t> h_offsets(num_ec + 1, 0);
  uint64_t offset = 0;
  for (size_t ec = 0; ec < num_ec; ++ec) {
    h_offsets[ec] = offset;
    if (ec_to_transcripts[ec] != nullptr) {
      for (uint32_t tr : *ec_to_transcripts[ec]) {
        h_transcripts.push_back(static_cast<int>(tr));
      }
      offset += ec_to_transcripts[ec]->cardinality();
    }
  }
  h_offsets[num_ec] = offset;
  num_ecs = num_ec;

  transcripts.assign(h_transcripts.begin(), h_transcripts.end());
  offsets.assign(h_offsets.begin(), h_offsets.end());
  cudaStreamSynchronize(0);

  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_gpu_ecmap_ms += setup_duration.count() / 1000.0;
}

GPUECMapInv::GPUECMapInv(const GPUIndex& index)
  : hash_map(std::max<size_t>(index.ecmapinv.size() * 2, 128),
             cuco::empty_key<uint64_t>{UINT64_MAX},
             cuco::empty_value<int>{-1}),
    num_ecs(0) {
  auto setup_start = std::chrono::high_resolution_clock::now();

  std::vector<cuco::pair<uint64_t, int>> host_pairs;
  host_pairs.reserve(index.ecmapinv.size());

  for (const auto& entry : index.ecmapinv) {
    uint64_t hash = hash_sorted_vector_well_defined_cpu(entry.first);
    host_pairs.push_back({hash, entry.second});
  }
  num_ecs = index.ecmapinv.size();

  thrust::device_vector<cuco::pair<uint64_t, int>> device_pairs(host_pairs.begin(), host_pairs.end());
  hash_map.insert(device_pairs.begin(), device_pairs.end());
  cudaStreamSynchronize(0);

  auto setup_end = std::chrono::high_resolution_clock::now();
  auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
  g_benchmark_stats.setup_gpu_ecmapinv_ms += setup_duration.count() / 1000.0;
}

void GPUECMapInv::insert_new_ecs(const thrust::device_vector<cuco::pair<uint64_t, int>>& pairs) {
  if (pairs.empty()) return;
  hash_map.insert(pairs.begin(), pairs.end());
  num_ecs += pairs.size();
}

void GPUECMap::append_ecs(const std::vector<int>& new_transcripts,
                          const std::vector<uint64_t>& new_offsets,
                          int num_new_ecs) {
  if (num_new_ecs == 0) return;

  size_t old_tx_size = transcripts.size();
  size_t old_num_ecs = num_ecs;

  // Append transcript data
  size_t new_tx_size = old_tx_size + new_transcripts.size();
  transcripts.resize(new_tx_size);
  thrust::copy(new_transcripts.begin(), new_transcripts.end(),
               transcripts.begin() + old_tx_size);

  // Append offsets: new_offsets are relative, shift by old_tx_size
  // new_offsets has num_new_ecs + 1 entries (like a standard offset array)
  size_t new_num_ecs = old_num_ecs + num_new_ecs;
  offsets.resize(new_num_ecs + 1);
  std::vector<uint64_t> shifted_offsets(num_new_ecs + 1);
  for (int i = 0; i <= num_new_ecs; ++i) {
    shifted_offsets[i] = new_offsets[i] + old_tx_size;
  }
  thrust::copy(shifted_offsets.begin(), shifted_offsets.end(),
               offsets.begin() + old_num_ecs);

  num_ecs = new_num_ecs;
  cudaStreamSynchronize(0);
}