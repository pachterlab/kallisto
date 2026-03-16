#pragma once

#include "MetalUtils.h"
#include "KmerIndex.h"
#include "GPUIndexFormat.h"

#include <vector>
#include <unordered_map>
#include <cstdint>

// Quick-and-dirty: extract 64-bit value from a Kmer (same as CUDA version)
inline uint64_t to_ullong_metal(const Kmer& km) {
    static_assert(sizeof(uint64_t) <= sizeof(Kmer), "Kmer must be >= 64 bits");
    return *reinterpret_cast<const uint64_t*>(&km);
}

// Hash function used on both CPU and GPU (FNV-1a over sorted int list)
uint64_t hash_sorted_vector_cpu(const int* vec, size_t n);
uint64_t hash_sorted_roaring_cpu(const Roaring& r);

// ============================================================================
// KmerSlot — flat open-addressing hash table entry
// 16-byte aligned, identical layout to the MSL struct in kallisto.metal
//
// fwd_run: # consecutive k-mers (including this one) that share the same EC
//          when traversing the de Bruijn graph in the FORWARD direction.
//          Use for SENSE reads (canonical == forward strand).  Min 1.
// bwd_run: same count going BACKWARD (i.e., skip for ANTISENSE reads).
//          Min 1.  Both capped at 255.
// This allows the process_reads_kernel to skip entire EC blocks with a single
// hash lookup, reducing typical lookups per read from ~120 to 2–5.
// ============================================================================

struct KmerSlot {
    uint64_t key;      // canonical k-mer value; empty_sentinel if empty
    int32_t  value;    // ec_id (or -1)
    uint8_t  fwd_run;  // forward run length  (1 = no skip possible)
    uint8_t  bwd_run;  // backward run length (1 = no skip possible)
    int16_t  _pad;
};
static_assert(sizeof(KmerSlot) == 16, "KmerSlot must be 16 bytes");

// ============================================================================
// MetalHashTable — bucket hash table (4 slots per 64-byte cache line)
//
// Capacity: 2^27 = 134M slots in 2^25 = 33.5M buckets × 64B = 2.14 GB
// Load factor: 116M / 134M ≈ 86.5%; avg 3.46 kmers/bucket
// Expected bucket probes: ~1.27 (vs 2.3 for per-slot double hashing)
//
// Each KmerSlot carries fwd_run and bwd_run so the per-read GPU kernel can
// skip entire EC blocks in one lookup (sense and antisense respectively),
// reducing lookups per read from ~120 → 2–5.  Working set per batch:
//   100K reads × 4 lookups × 64B ≈ 25 MB  (fits in M4 Max SLC)
//   hot transcripts subset ≈ 1–8 MB        (fits in M1 8MB SLC)
// ============================================================================

struct MetalHashTable {
    MetalBuffer<KmerSlot> slots;      // flat hash table
    size_t                capacity = 0;   // always a power of 2
    size_t                num_kmers = 0;
    uint64_t              empty_sentinel = 0;

    MetalHashTable() = default;

    /// Build from a KmerIndex (collect pairs, insert with double hashing)
    void build(const KmerIndex& index);

    /// Build from a GPUIndex
    void build(const GPUIndex& index);

    /// CPU double-hash lookup — returns EC ID or -1
    int lookup(uint64_t kmer) const;

    /// Save/load binary sidecar. Returns false on failure.
    bool save_cache(const std::string& path) const;
    bool load_cache(const std::string& path);
};

// Alias so old code that refers to MetalKmerTable still compiles
using MetalKmerTable = MetalHashTable;
using SortedKmerTable = MetalHashTable;

// ============================================================================
// MetalECMap — forward map: ec_id -> transcript list (CSR layout)
// Both CPU and GPU access via shared MTLBuffers
// ============================================================================

struct MetalECMap {
    MetalBuffer<int>     transcripts;   // flat transcript IDs
    MetalBuffer<uint64_t> offsets;      // offsets[ec_id] .. offsets[ec_id+1]
    size_t               num_ecs = 0;

    MetalECMap() = default;
    void build(const KmerIndex& index);
    void build(const GPUIndex& index);

    /// Append novel ECs (relative offsets will be shifted)
    void append_ecs(const std::vector<int>& new_transcripts,
                    const std::vector<uint64_t>& new_offsets,
                    int num_new_ecs);
};

// ============================================================================
// MetalECMapInv — inverse map: FNV-1a hash of transcript set -> ec_id
// CPU-side only (std::unordered_map); GPU queries hash-then-CPU-verify
// ============================================================================

struct MetalECMapInv {
    std::unordered_map<uint64_t, int> map;
    size_t num_ecs = 0;

    MetalECMapInv() = default;
    void build(const KmerIndex& index);
    void build(const GPUIndex& index);

    /// Look up; returns -1 if not found
    int find(uint64_t hash) const;

    /// Insert new EC (hash -> ec_id)
    void insert(uint64_t hash, int ec_id);
};
