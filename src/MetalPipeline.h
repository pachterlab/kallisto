#pragma once

#include "MetalUtils.h"
#include "MetalIndex.h"
#include "BenchmarkStats.h"
#include "ProcessReads.h"

#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <iostream>

// ============================================================================
// HostReadLoader — CPU-side FASTQ reader (identical logic to CUDA version)
//
// Note: load() takes an external FastqSequenceReader so two HostReadLoader
// instances can share a single reader for correct double-buffering.
// ============================================================================

struct HostReadLoader {
    std::vector<char>     reads;
    std::vector<uint64_t> read_id_to_offset;
    std::vector<uint32_t> read_id_to_length;
    std::vector<uint64_t> read_id_to_kmer_first;
    std::vector<uint32_t> read_id_to_kmer_count;
    uint64_t kmer_count  = 0;
    uint32_t read_count  = 0;
    uint32_t r1_count    = 0;
    uint32_t pair_limit  = 0;
    bool     is_paired   = false;

    // Persistent fetch buffer — sequences in pending_seqs point into this buffer.
    // We only overwrite fetch_buf (via fetchSequences) when pending_seqs is exhausted.
    std::vector<char>                        fetch_buf;
    std::vector<std::pair<const char*, int>> pending_seqs;
    size_t                                   pending_offset   = 0;
    bool                                     pending_has_more = true;
    int                                      readbatch_id     = -1;

    HostReadLoader(const ProgramOptions& opt,
                   uint32_t pair_limit_,
                   uint32_t read_length_hint = 200);

    bool load(FastqSequenceReader& reader);
};

// ============================================================================
// DeviceKmerLoader — uploads reads to shared Metal buffers and runs kmer_kernel
// ============================================================================

struct DeviceKmerLoader {
    MetalBuffer<char>     reads;
    MetalBuffer<uint64_t> read_id_to_offset;
    MetalBuffer<uint32_t> read_id_to_length;
    MetalBuffer<uint64_t> read_id_to_kmer_first;
    MetalBuffer<uint32_t> read_id_to_kmer_count;
    MetalBuffer<uint64_t> kmers;
    uint64_t read_count = 0;
    uint32_t r1_count   = 0;

    DeviceKmerLoader() = default;

    void load(const HostReadLoader& h);
    void run(uint64_t empty_kmer_value, uint32_t k);
};

// ============================================================================
// ReadECCollapser — dedup + sort ECs per read, produce ragged array
// ============================================================================

struct ReadECCollapser {
    MetalBuffer<int>      read_ecs;
    MetalBuffer<uint64_t> read_ec_offsets;
    MetalBuffer<int>      temp_ecs_;
    uint64_t read_count = 0;

    ReadECCollapser() = default;

    // Original path: separate kmer_kernel + lookup kernel have already run.
    void collapse(const DeviceKmerLoader& loader,
                  const MetalBuffer<int>& d_ecs);

    // New optimised path: single process_reads_kernel replaces kmer_kernel,
    // kmer_lookup_kernel, and compute_ec_collapse_sizes_kernel in one pass.
    // Uses run-length skipping to reduce lookups from ~120/read to ~2–5/read.
    void process_reads(const DeviceKmerLoader& loader,
                       const MetalHashTable& kmer_table,
                       uint32_t k);
};

// ============================================================================
// ReadTranscriptIntersector — intersect per-read ECs with EC->transcript map
// ============================================================================

struct ReadTranscriptIntersector {
    MetalBuffer<int>      read_transcripts;
    MetalBuffer<uint64_t> read_transcript_offsets;
    MetalBuffer<uint64_t> read_transcript_sizes;
    MetalBuffer<uint8_t>  read_had_mapped_kmers;
    uint64_t read_count = 0;

    ReadTranscriptIntersector() = default;

    void intersect(const ReadECCollapser& collapser, const MetalECMap& ecmap);
};

// ============================================================================
// PairIntersector — intersect R1 and R2 transcript sets
// ============================================================================

struct PairIntersector {
    MetalBuffer<int>      pair_transcripts;
    MetalBuffer<uint64_t> pair_transcript_offsets;
    MetalBuffer<uint64_t> pair_transcript_sizes;
    uint32_t pair_count = 0;

    PairIntersector() = default;

    void intersect_pairs(const ReadTranscriptIntersector& intersector,
                         uint32_t r1_count);
};

// ============================================================================
// ReadECLookup — hash transcript sets, look up in inverse EC map (CPU-side)
// ============================================================================

struct ReadECLookup {
    MetalBuffer<int>      read_ecs_final;
    MetalBuffer<uint64_t> read_hashes;
    uint64_t read_count = 0;

    ReadECLookup() = default;

    void lookup(const ReadTranscriptIntersector& intersector,
                MetalECMapInv& ecmapinv,
                const MetalECMap& ecmap);

    void lookup_pairs(const PairIntersector& pairs,
                      MetalECMapInv& ecmapinv,
                      const MetalECMap& ecmap);

private:
    void lookup_impl(const MetalBuffer<int>&      transcripts,
                     const MetalBuffer<uint64_t>& offsets,
                     const MetalBuffer<uint64_t>& sizes,
                     uint64_t count,
                     MetalECMapInv& ecmapinv,
                     const MetalECMap& ecmap);
};

// ============================================================================
// ECCounter — accumulate EC counts across batches, write final output
// ============================================================================

struct ECCounter {
    MetalBuffer<int> ec_counts;
    size_t num_ecs = 0;

    explicit ECCounter(size_t n);

    void ensure_capacity(size_t new_num_ecs);
    void count_batch(const ReadECLookup& lookup);
    void write_counts(const std::string& filename);
};

// ============================================================================
// NewECHandler — CPU-side handler for novel paired-end ECs
// ============================================================================

struct NewECHandler {
    int next_ec_id          = 0;
    int total_new_ec_count  = 0;
    int total_new_read_count = 0;

    explicit NewECHandler(int initial_num_ecs) : next_ec_id(initial_num_ecs) {}

    void handle_batch(const PairIntersector& pairs,
                      ReadECLookup&   ec_lookup,
                      MetalECMapInv&  ecmapinv,
                      MetalECMap&     ecmap,
                      ECCounter&      ec_counter);

    int total_new_ecs()   const { return total_new_ec_count; }
    int total_new_reads() const { return total_new_read_count; }
};
