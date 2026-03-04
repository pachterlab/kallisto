#ifndef GPU_KERNELS_CUH
#define GPU_KERNELS_CUH

#include <cstdint>
#include <cuda_runtime.h>


// Combined encode + validity check: returns 0-3 for valid bases, 4 for invalid
// This avoids separate is_valid() calls and reduces branches
__device__ __forceinline__ uint8_t encode_base_or_invalid(char b) {
  // Using bit manipulation: A=0x41, C=0x43, G=0x47, T=0x54
  // Check if it's one of ACGT first
  uint8_t is_acg = (b == 'A') | (b == 'C') | (b == 'G');
  uint8_t is_t = (b == 'T');
  if (is_acg | is_t) {
    // A(0x41)->0, C(0x43)->1, G(0x47)->2, T(0x54)->3
    // Formula: ((b >> 1) & 0x3) with fixup for G and T
    return (b == 'A') ? 0 : (b == 'C') ? 1 : (b == 'G') ? 2 : 3;
  }
  return 4;  // Invalid marker
}

__device__ uint64_t hash_sorted_vector_well_defined_device(const int* __restrict__ vec, uint64_t size);
__device__ bool verify_transcript_lists_equal(
    const int* __restrict__ list1, uint64_t size1,
    const int* __restrict__ list2, uint64_t size2);
__device__ int64_t binary_search_hash(const uint64_t* __restrict__ sorted_hashes, uint64_t num_hashes, uint64_t target_hash);

// Kernel declarations
__global__ void kmer_kernel(
    const char* __restrict__ reads,
    const uint64_t* __restrict__ offsets,
    const uint32_t* __restrict__ lengths,
    const uint32_t* __restrict__ kmers_counts,
    const uint64_t* __restrict__ kmers_offsets,
    uint64_t* __restrict__ out_kmers,
    uint64_t empty_kmer_value,
    uint32_t num_reads,
    uint32_t k);

// MAX_ECS_PER_READ must match the value in GPUKernels.cu
#ifndef MAX_ECS_PER_READ
#define MAX_ECS_PER_READ 256
#endif

// First kernel: dedups, sorts, caches results in temp_ecs buffer
__global__ void compute_ec_collapse_sizes_kernel(
    const uint64_t* __restrict__ read_id_to_kmer_first,
    const uint32_t* __restrict__ read_id_to_kmer_count,
    const int* __restrict__ ecs,
    uint64_t* __restrict__ output_sizes,
    int* __restrict__ temp_ecs,  // Output: cached unique ECs (MAX_ECS_PER_READ * num_reads)
    uint64_t num_reads);

// Second kernel: just copies from temp_ecs to output (no dedup needed)
__global__ void collapse_ecs_per_read_kernel(
    const int* __restrict__ temp_ecs,  // Cached unique ECs from first kernel
    const uint64_t* __restrict__ output_offsets,
    int* __restrict__ output_ecs,
    uint64_t num_reads);

__global__ void compute_intersection_sizes_kernel(
    const int* __restrict__ read_ecs,
    const uint64_t* __restrict__ read_ec_offsets,
    const uint64_t* __restrict__ ecmap_offsets,
    uint64_t* __restrict__ output_sizes,
    uint64_t* __restrict__ smallest_ec_indices,    // Output: index of smallest EC per read
    uint64_t num_reads,
    uint64_t num_ecs,
    uint64_t max_transcripts,
    uint64_t read_ecs_size);

__global__ void intersect_transcripts_kernel(
    const int* __restrict__ read_ecs,
    const uint64_t* __restrict__ read_ec_offsets,
    const int* __restrict__ ecmap_transcripts,
    const uint64_t* __restrict__ ecmap_offsets,
    const uint64_t* __restrict__ output_offsets,
    int* __restrict__ output_transcripts,
    uint64_t* __restrict__ output_sizes,
    const uint64_t* __restrict__ smallest_ec_indices,  // Input: pre-computed smallest EC index
    uint64_t num_reads,
    uint64_t num_ecs,
    uint64_t max_transcripts,
    uint64_t read_ecs_size,
    uint64_t max_output_transcripts);

// Paired-end intersection kernels
// Pass 1: compute the size of the intersection for each R1/R2 pair
__global__ void compute_pair_intersection_sizes_kernel(
    const int* __restrict__ read_transcripts,
    const uint64_t* __restrict__ read_transcript_offsets,
    const uint64_t* __restrict__ read_transcript_sizes,
    const uint8_t* __restrict__ read_had_mapped_kmers,
    uint64_t* __restrict__ pair_sizes,
    uint32_t r1_count);

// Pass 2: write the actual intersection for each pair
__global__ void intersect_pairs_kernel(
    const int* __restrict__ read_transcripts,
    const uint64_t* __restrict__ read_transcript_offsets,
    const uint64_t* __restrict__ read_transcript_sizes,
    const uint8_t* __restrict__ read_had_mapped_kmers,
    const uint64_t* __restrict__ pair_offsets,
    int* __restrict__ pair_transcripts,
    uint64_t* __restrict__ pair_sizes_out,
    uint32_t r1_count);

// Contig k-mer extraction: for each EC block, extract k-mers and output (canonical_kmer, ec_id)
// block_seq_offsets: start position of each block in flat sequence
// block_num_kmers: number of k-mers per block
// block_output_offsets: output index for each block (cumulative)
// Output: out_pairs - array of (kmer, ec_id) pairs
__global__ void contig_kmer_kernel(
    const char* __restrict__ seq,
    const uint64_t* __restrict__ block_seq_offsets,
    const uint32_t* __restrict__ block_num_kmers,
    const int* __restrict__ block_ec_ids,
    const uint64_t* __restrict__ block_output_offsets,
    uint64_t* __restrict__ out_kmers,
    int* __restrict__ out_ecs,
    uint64_t empty_kmer_value,
    int empty_ec_value,
    uint32_t k,
    uint64_t num_blocks);

#endif // GPU_KERNELS_CUH

