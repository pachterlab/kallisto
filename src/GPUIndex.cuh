#ifndef GPU_INDEX_CUH
#define GPU_INDEX_CUH

#include "KmerIndex.h"
#include "GPUIndexFormat.h"
#include <cuco/static_map.cuh>
#include <cuco/dynamic_map.cuh>
#include <thrust/device_vector.h>
#include <vector>
#include <cstdint>
#include <cstddef>

// Forward declarations
struct GPUIndex;

// Quick-and-dirty helper: interpret a Kmer as a 64-bit value by taking
// the first 64 bits of its storage. This is equivalent to accessing
// longs[0] in the original Kmer implementation and is only used for
// providing a sentinel "empty" k-mer value to GPU code.
inline uint64_t to_ullong(const Kmer& km) {
  static_assert(sizeof(uint64_t) <= sizeof(Kmer),
                "Kmer must be 64 bits wide");
  return *reinterpret_cast<const uint64_t*>(&km);
}

// Forward declarations
struct GPUECMap;
struct GPUECMapInv;

// Hash table building
cuco::static_map<uint64_t, int> build_kmer_to_ec_map(const KmerIndex& index);
cuco::static_map<uint64_t, int> build_kmer_to_ec_map_from_contigs(const GPUIndex& index);

// CPU hash function: compute a stable hash from a (sorted) Roaring bitmap
// of transcript IDs. Iteration over Roaring yields transcript IDs in
// sorted order, so this is well-defined.
uint64_t hash_sorted_vector_well_defined_cpu(const Roaring& r);

// EC Map structures (forward map: ec_id -> transcript list)
struct GPUECMap {
  thrust::device_vector<int> transcripts;
  thrust::device_vector<uint64_t> offsets;
  size_t num_ecs;
  
  GPUECMap(const KmerIndex& index);
  GPUECMap(const GPUIndex& index);

  // Append new ECs (transcript data + offsets) for novel equivalence classes
  void append_ecs(const std::vector<int>& new_transcripts,
                  const std::vector<uint64_t>& new_offsets,
                  int num_new_ecs);
};

// Inverse EC map: transcript-set-hash -> ec_id, using cuco::dynamic_map
// for GPU-side incremental insertion of novel ECs.
struct GPUECMapInv {
  cuco::dynamic_map<uint64_t, int> hash_map;
  size_t num_ecs;

  GPUECMapInv(const KmerIndex& index);
  GPUECMapInv(const GPUIndex& index);

  void insert_new_ecs(const thrust::device_vector<cuco::pair<uint64_t, int>>& pairs);
};

#endif // GPU_INDEX_CUH

