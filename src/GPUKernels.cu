#include "GPUKernels.cuh"

#ifndef MAX_ECS_PER_READ
#define MAX_ECS_PER_READ 512
#endif

// Device utility function implementations (non-inline)
__device__ uint64_t hash_sorted_vector_well_defined_device(const int* __restrict__ vec, uint64_t size) {
  uint64_t hash = 0xcbf29ce484222325ULL;  // FNV-1a offset basis
  for (uint64_t i = 0; i < size; ++i) {
    hash ^= static_cast<uint64_t>(vec[i]);
    hash *= 0x100000001b3ULL;  // FNV-1a prime
  }
  // Avoid sentinel value collision
  if (hash == UINT64_MAX) {
    hash = UINT64_MAX - 2;
  }
  return hash;
}

__device__ bool verify_transcript_lists_equal(
    const int* __restrict__ list1, uint64_t size1,
    const int* __restrict__ list2, uint64_t size2) {
  if (size1 != size2) {
    return false;
  }
  
  for (uint64_t i = 0; i < size1; ++i) {
    if (list1[i] != list2[i]) {
      return false;
    }
  }
  
  return true;
}

__device__ int64_t binary_search_hash(const uint64_t* __restrict__ sorted_hashes, uint64_t num_hashes, uint64_t target_hash) {
  int64_t left = 0;
  int64_t right = static_cast<int64_t>(num_hashes) - 1;
  
  while (left <= right) {
    int64_t mid = left + (right - left) / 2;
    uint64_t mid_hash = sorted_hashes[mid];
    
    if (mid_hash == target_hash) {
      return mid;
    } else if (mid_hash < target_hash) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
  
  return -1;
}

__global__ void kmer_kernel(const char* __restrict__ reads,
                            const uint64_t* __restrict__ offsets,
                            const uint32_t* __restrict__ lengths,
                            const uint32_t* __restrict__ kmers_counts,
                            const uint64_t* __restrict__ kmers_offsets,
                            uint64_t* __restrict__ out_kmers,
                            uint64_t empty_kmer_value,
                            uint32_t num_reads,
                            uint32_t k)
{
  uint32_t r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= num_reads) return;

  const char* read   = reads + offsets[r];
  uint32_t    L      = lengths[r];
  uint32_t    n_kmers = kmers_counts[r];
  uint64_t    out0   = kmers_offsets[r];

  if (L < k || n_kmers == 0) {
    return;
  }

  // Precompute bit positions and masks (constant for all k-mers in this read)
  const uint32_t fwd_new_bit_pos = 2 * (32 - k);  // Where new base goes in fwd
  const uint32_t rev_high_bit_pos = 62;           // Where new complement goes in rev (always position 0 = bits 62-63)
  // Mask to clear low bits outside k-mer region after right shift
  // K-mer uses bits [2*(32-k), 63], so we need to clear bits [0, 2*(32-k)-1]
  const uint64_t low_bits_mask = ~((1ULL << fwd_new_bit_pos) - 1);

  uint64_t fwd = 0;
  uint64_t rev = 0;
  int      bad = 0;

  // Build first k-mer (both forward and reverse complement)
  for (uint32_t i = 0; i < k; ++i) {
    char c = read[i];
    uint8_t code = encode_base_or_invalid(c);
    if (code > 3) {
      ++bad;
      code = 0;  // Use 'A' encoding for invalid bases
    }
    
    // Forward: bases go left to right in high bits
    uint32_t fwd_bit_pos = 2 * (31 - i);
    fwd |= ((uint64_t)code) << fwd_bit_pos;
    
    // Reverse complement: complement bases go right to left
    // Position i in forward becomes position (k-1-i) in reverse
    uint8_t comp_code = code ^ 0x3;  // A<->T (0<->3), C<->G (1<->2)
    uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
    rev |= ((uint64_t)comp_code) << rev_bit_pos;
  }

  uint64_t canon0 = (fwd < rev) ? fwd : rev;
  out_kmers[out0] = (bad > 0) ? empty_kmer_value : canon0;

  // Process remaining k-mers with incremental updates
  for (uint32_t i = 1; i < n_kmers; ++i) {
    // Check outgoing base (leaving the window)
    char c_out = read[i - 1];
    uint8_t code_out = encode_base_or_invalid(c_out);
    if (code_out > 3) {
      --bad;
    }

    // Check incoming base (entering the window)
    char c_in = read[i + k - 1];
    uint8_t code_in = encode_base_or_invalid(c_in);
    if (code_in > 3) {
      ++bad;
      code_in = 0;  // Use 'A' encoding for invalid bases
    }

    // Update forward k-mer: shift left, add new base at low end
    fwd = (fwd << 2) | (((uint64_t)code_in) << fwd_new_bit_pos);

    // Update reverse complement incrementally:
    uint8_t comp_code_in = code_in ^ 0x3;
    rev = ((rev >> 2) | (((uint64_t)comp_code_in) << rev_high_bit_pos)) & low_bits_mask;

    uint64_t canon = (fwd < rev) ? fwd : rev;
    out_kmers[out0 + i] = (bad > 0) ? empty_kmer_value : canon;
  }
}

// EC collapse with caching:
// - Writes unique ECs to temp_ecs buffer (avoids redo in second kernel)
// - Uses linear search for small counts, binary search for larger
// - Adjacent duplicate fast path
// - Array stays sorted as we insert
__global__ void compute_ec_collapse_sizes_kernel(
    const uint64_t* __restrict__ read_id_to_kmer_first,
    const uint32_t* __restrict__ read_id_to_kmer_count,
    const int* __restrict__ ecs,
    uint64_t* __restrict__ output_sizes,
    int* __restrict__ temp_ecs,  // Output: cached unique ECs per read (MAX_ECS_PER_READ * num_reads)
    uint64_t num_reads) {
  
  uint64_t read_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (read_id >= num_reads) return;
  
  uint64_t kmer_start = read_id_to_kmer_first[read_id];
  uint32_t kmer_count = read_id_to_kmer_count[read_id];
  
  if (kmer_count == 0) {
    output_sizes[read_id] = 0;
    return;
  }
  
  int unique_ecs[MAX_ECS_PER_READ];
  int count = 0;
  int prev_ec = -1;  // For adjacent duplicate detection
  
  for (uint32_t i = 0; i < kmer_count; ++i) {
    int ec = ecs[kmer_start + i];
    
    if (ec == -1) continue;
    
    // Fast path: skip if same as previous (common for adjacent k-mers)
    if (ec == prev_ec) continue;
    prev_ec = ec;
    
    // For small counts, linear search is faster than binary search
    if (count < 8) {
      bool found = false;
      int insert_pos = count;  // Default: append at end
      for (int j = 0; j < count; ++j) {
        if (unique_ecs[j] == ec) {
          found = true;
          break;
        }
        if (unique_ecs[j] > ec && insert_pos == count) {
          insert_pos = j;
        }
      }
      if (found) continue;
      
      // Insert at insert_pos
      if (count < MAX_ECS_PER_READ) {
        for (int j = count; j > insert_pos; --j) {
          unique_ecs[j] = unique_ecs[j-1];
        }
        unique_ecs[insert_pos] = ec;
        count++;
      }
    } else {
      // Binary search for insertion point
      int lo = 0, hi = count;
      while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (unique_ecs[mid] < ec) lo = mid + 1;
        else hi = mid;
      }
      
      // Check if duplicate
      if (lo < count && unique_ecs[lo] == ec) continue;
      
      // Insert at position lo
      if (count < MAX_ECS_PER_READ) {
        for (int j = count; j > lo; --j) {
          unique_ecs[j] = unique_ecs[j-1];
        }
        unique_ecs[lo] = ec;
        count++;
      }
    }
  }
  
  output_sizes[read_id] = count;
  
  // Cache the unique ECs for the second kernel (avoid redoing dedup)
  int* my_temp = temp_ecs + read_id * MAX_ECS_PER_READ;
  for (int i = 0; i < count; ++i) {
    my_temp[i] = unique_ecs[i];
  }
}

// Fast copy kernel: just copies cached unique ECs from temp buffer to output
// All dedup work was done in compute_ec_collapse_sizes_kernel
__global__ void collapse_ecs_per_read_kernel(
    const int* __restrict__ temp_ecs,  // Cached unique ECs from first kernel
    const uint64_t* __restrict__ output_offsets,
    int* __restrict__ output_ecs,
    uint64_t num_reads) {
  
  uint64_t read_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (read_id >= num_reads) return;
  
  uint64_t output_start = output_offsets[read_id];
  uint64_t output_size = output_offsets[read_id + 1] - output_start;
  
  if (output_size == 0) {
    return;
  }
  
  // Just copy from temp buffer (already sorted, already deduped)
  const int* my_temp = temp_ecs + read_id * MAX_ECS_PER_READ;
  for (uint64_t i = 0; i < output_size; ++i) {
    output_ecs[output_start + i] = my_temp[i];
  }
}

__global__ void compute_intersection_sizes_kernel(
    const int* __restrict__ read_ecs,
    const uint64_t* __restrict__ read_ec_offsets,
    const uint64_t* __restrict__ ecmap_offsets,
    uint64_t* __restrict__ output_sizes,
    uint64_t* __restrict__ smallest_ec_indices,
    uint64_t num_reads,
    uint64_t num_ecs,
    uint64_t max_transcripts,
    uint64_t read_ecs_size) {
  
  uint64_t read_id = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (read_id >= num_reads) {
    return;
  }
  
  output_sizes[read_id] = 0;
  smallest_ec_indices[read_id] = 0;
  
  uint64_t ec_start = read_ec_offsets[read_id];
  uint64_t ec_end = read_ec_offsets[read_id + 1];
  
  if (ec_end <= ec_start || ec_start >= read_ecs_size) {
    return;
  }
  
  uint64_t ec_count = ec_end - ec_start;
  
  // Check first EC for validity
  int first_ec = read_ecs[ec_start];
  
  uint64_t first_tx_start = ecmap_offsets[first_ec];
  uint64_t first_tx_end = ecmap_offsets[first_ec + 1];
  
  if (ec_count == 1) {
    output_sizes[read_id] = first_tx_end - first_tx_start;
    smallest_ec_indices[read_id] = ec_start;
    return;
  }
  
  // Find minimum EC size and its index across all ECs for this read
  uint64_t first_ec_size = first_tx_end - first_tx_start;
  uint64_t min_ec_size = first_ec_size;
  uint64_t min_ec_idx = ec_start;
  
  for (uint64_t i = ec_start + 1; i < ec_end; ++i) {
    int ec = read_ecs[i];
    
    uint64_t ec_uint = static_cast<uint64_t>(ec);
    uint64_t tx_start = ecmap_offsets[ec_uint];
    uint64_t tx_end = ecmap_offsets[ec_uint + 1];
    
    uint64_t ec_size = tx_end - tx_start;
    if (ec_size < min_ec_size) {
      min_ec_size = ec_size;
      min_ec_idx = i;
    }
  }
  
  output_sizes[read_id] = min_ec_size;
  smallest_ec_indices[read_id] = min_ec_idx;
}

__global__ void intersect_transcripts_kernel(
    const int* __restrict__ read_ecs,
    const uint64_t* __restrict__ read_ec_offsets,
    const int* __restrict__ ecmap_transcripts,
    const uint64_t* __restrict__ ecmap_offsets,
    const uint64_t* __restrict__ output_offsets,
    int* __restrict__ output_transcripts,
    uint64_t* __restrict__ output_sizes,
    const uint64_t* __restrict__ smallest_ec_indices,
    uint64_t num_reads,
    uint64_t num_ecs,
    uint64_t max_transcripts,
    uint64_t read_ecs_size,
    uint64_t max_output_transcripts) {
  
  uint64_t read_id = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (read_id >= num_reads) {
    return;
  }
  
  uint64_t ec_start = read_ec_offsets[read_id];
  uint64_t ec_end = read_ec_offsets[read_id + 1];
  
  if (ec_end <= ec_start || ec_start >= read_ecs_size) {
    output_sizes[read_id] = 0;
    return;
  }
  
  uint64_t ec_count = ec_end - ec_start;
  uint64_t output_start = output_offsets[read_id];
  uint64_t output_end = output_offsets[read_id + 1];
  
  if (output_end <= output_start) {
    output_sizes[read_id] = 0;
    return;
  }
  
  // Use pre-computed smallest EC index from compute_intersection_sizes_kernel
  uint64_t smallest_ec_idx = smallest_ec_indices[read_id];
  
  int smallest_ec = read_ecs[smallest_ec_idx];
  
  uint64_t smallest_tx_start = ecmap_offsets[smallest_ec];
  uint64_t smallest_tx_end = ecmap_offsets[smallest_ec + 1];
  
  uint64_t smallest_ec_size = smallest_tx_end - smallest_tx_start;
  
  if (ec_count == 1) {
    uint64_t idx = 0;
    uint64_t output_limit = (output_end < output_start + max_output_transcripts) ? output_end : output_start + max_output_transcripts;
    for (uint64_t i = smallest_tx_start; i < smallest_tx_end && output_start + idx < output_limit; ++i, ++idx) {
      output_transcripts[output_start + idx] = ecmap_transcripts[i];
    }
    output_sizes[read_id] = idx;
    return;
  }
  
  uint64_t current_size = smallest_ec_size;
  uint64_t output_limit = (output_end < output_start + max_output_transcripts) ? output_end : output_start + max_output_transcripts;
  
  // Copy smallest EC transcripts (this fits in min_ec_size allocation)
  uint64_t copy_count = (current_size < output_limit - output_start) ? current_size : output_limit - output_start;
  for (uint64_t i = 0; i < copy_count; ++i) {
    output_transcripts[output_start + i] = ecmap_transcripts[smallest_tx_start + i];
  }
  current_size = copy_count;
  
  // Early exit for empty intersections
  if (current_size == 0) {
    output_sizes[read_id] = 0;
    return;
  }
  
  // Intersect with remaining ECs (skip the one we started with)
  for (uint64_t ec_idx = ec_start; ec_idx < ec_end && current_size > 0; ++ec_idx) {
    // Skip the EC we started with
    if (ec_idx == smallest_ec_idx) {
      continue;
    }
    
    int ec = read_ecs[ec_idx];
    
    uint64_t tx_start = ecmap_offsets[ec];
    uint64_t tx_end = ecmap_offsets[ec + 1];
    
    uint64_t ec_size = tx_end - tx_start;
    uint64_t write_idx = 0;
    
    // Early termination: check if ranges overlap at all
    // Both arrays are sorted, so if max(A) < min(B) or max(B) < min(A), intersection is empty
    int first_current = output_transcripts[output_start];
    int last_current = output_transcripts[output_start + current_size - 1];
    int first_ec = ecmap_transcripts[tx_start];
    int last_ec = ecmap_transcripts[tx_end - 1];
    
    if (last_current < first_ec || last_ec < first_current) {
      current_size = 0;
      break;
    }
    
    // Use binary search when EC is significantly larger than current intersection
    // Ratio-based threshold: binary search when ec_size > current_size * 5
    // Binary search: O(current_size * log(ec_size)) vs Two-pointer: O(current_size + ec_size)
    if (ec_size > current_size * 5) {
      for (uint64_t i = 0; i < current_size && output_start + write_idx < output_limit; ++i) {
        int val = output_transcripts[output_start + i];
        
        // Binary search for val in EC's transcript list
        int64_t lo = 0, hi = static_cast<int64_t>(ec_size) - 1;
        bool found = false;
        while (lo <= hi) {
          int64_t mid = lo + (hi - lo) / 2;
          int b = ecmap_transcripts[tx_start + mid];
          if (b == val) { found = true; break; }
          else if (b < val) lo = mid + 1;
          else hi = mid - 1;
        }
        
        if (found) {
          output_transcripts[output_start + write_idx++] = val;
        }
      }
    } else {
      // Two-pointer merge for larger intersections or small EC lists
      uint64_t read_i = 0;
      uint64_t j = tx_start;
      
      while (read_i < current_size && j < tx_end && output_start + write_idx < output_limit) {
        int a = output_transcripts[output_start + read_i];
        int b = ecmap_transcripts[j];
        
        if (a < b) {
          ++read_i;  // Skip a, not in intersection
        } else if (b < a) {
          ++j;       // Skip b, not in intersection
        } else {
          // Match found - write directly to compacted position
          output_transcripts[output_start + write_idx] = a;
          ++write_idx;
          ++read_i;
          ++j;
        }
      }
    }
    
    // Update size after intersection
    current_size = write_idx;
    if (current_size == 0) {
      break;
    }
  }
  
  output_sizes[read_id] = current_size;
}

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
    uint64_t num_blocks) {
  uint64_t block_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (block_id >= num_blocks) return;

  const char* block_seq = seq + block_seq_offsets[block_id];
  uint32_t n_kmers = block_num_kmers[block_id];
  int ec_id = block_ec_ids[block_id];
  uint64_t out_base = block_output_offsets[block_id];

  if (n_kmers == 0) return;

  const uint32_t fwd_new_bit_pos = 2 * (32 - k);
  const uint32_t rev_high_bit_pos = 62;
  const uint64_t low_bits_mask = ~((1ULL << fwd_new_bit_pos) - 1);

  uint64_t fwd = 0;
  uint64_t rev = 0;
  int bad = 0;

  for (uint32_t i = 0; i < k; ++i) {
    char c = block_seq[i];
    uint8_t code = encode_base_or_invalid(c);
    if (code > 3) {
      ++bad;
      code = 0;
    }
    uint32_t fwd_bit_pos = 2 * (31 - i);
    fwd |= ((uint64_t)code) << fwd_bit_pos;
    uint8_t comp_code = code ^ 0x3;
    uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
    rev |= ((uint64_t)comp_code) << rev_bit_pos;
  }

  uint64_t canon = (bad == 0) ? ((fwd < rev) ? fwd : rev) : empty_kmer_value;
  out_kmers[out_base] = canon;
  out_ecs[out_base] = (bad == 0) ? ec_id : empty_ec_value;
  ++out_base;

  for (uint32_t i = 1; i < n_kmers; ++i) {
    char c_out = block_seq[i - 1];
    uint8_t code_out = encode_base_or_invalid(c_out);
    if (code_out > 3) --bad;

    char c_in = block_seq[i + k - 1];
    uint8_t code_in = encode_base_or_invalid(c_in);
    if (code_in > 3) {
      ++bad;
      code_in = 0;
    }

    fwd = (fwd << 2) | (((uint64_t)code_in) << fwd_new_bit_pos);
    uint8_t comp_code_in = code_in ^ 0x3;
    rev = ((rev >> 2) | (((uint64_t)comp_code_in) << rev_high_bit_pos)) & low_bits_mask;

    canon = (bad == 0) ? ((fwd < rev) ? fwd : rev) : empty_kmer_value;
    out_kmers[out_base] = canon;
    out_ecs[out_base] = (bad == 0) ? ec_id : empty_ec_value;
    ++out_base;
  }
}

// =============================================================================
// Paired-end intersection kernels
// =============================================================================

// Pass 1: compute the size of the intersection for each R1/R2 pair
// For pair i: R1 = read i, R2 = read (r1_count + i)
// Fallback logic (matches CPU kallisto):
//   both empty -> 0
//   only R1 has transcripts -> use R1 size
//   only R2 has transcripts -> use R2 size
//   both non-empty -> sorted set intersection size
__global__ void compute_pair_intersection_sizes_kernel(
    const int* __restrict__ read_transcripts,
    const uint64_t* __restrict__ read_transcript_offsets,
    const uint64_t* __restrict__ read_transcript_sizes,
    const uint8_t* __restrict__ read_had_mapped_kmers,
    uint64_t* __restrict__ pair_sizes,
    uint32_t r1_count)
{
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= r1_count) return;

    uint32_t r1_idx = i;
    uint32_t r2_idx = r1_count + i;

    uint64_t r1_size = read_transcript_sizes[r1_idx];
    uint64_t r2_size = read_transcript_sizes[r2_idx];

    if (r1_size == 0 && r2_size == 0) {
        pair_sizes[i] = 0;
        return;
    }
    if (r1_size == 0) {
        if (read_had_mapped_kmers[r1_idx]) {
            pair_sizes[i] = 0;
        } else {
            pair_sizes[i] = r2_size;
        }
        return;
    }
    if (r2_size == 0) {
        if (read_had_mapped_kmers[r2_idx]) {
            pair_sizes[i] = 0;
        } else {
            pair_sizes[i] = r1_size;
        }
        return;
    }

    // Both non-empty: compute sorted set intersection size
    const int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
    const int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];

    uint64_t count = 0;
    uint64_t a = 0, b = 0;
    while (a < r1_size && b < r2_size) {
        int va = r1_list[a];
        int vb = r2_list[b];
        if (va == vb) {
            ++count;
            ++a;
            ++b;
        } else if (va < vb) {
            ++a;
        } else {
            ++b;
        }
    }

    pair_sizes[i] = count;
}

// Pass 2: write the actual intersection (or fallback copy) for each pair
__global__ void intersect_pairs_kernel(
    const int* __restrict__ read_transcripts,
    const uint64_t* __restrict__ read_transcript_offsets,
    const uint64_t* __restrict__ read_transcript_sizes,
    const uint8_t* __restrict__ read_had_mapped_kmers,
    const uint64_t* __restrict__ pair_offsets,
    int* __restrict__ pair_transcripts,
    uint64_t* __restrict__ pair_sizes_out,
    uint32_t r1_count)
{
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= r1_count) return;

    uint32_t r1_idx = i;
    uint32_t r2_idx = r1_count + i;

    uint64_t r1_size = read_transcript_sizes[r1_idx];
    uint64_t r2_size = read_transcript_sizes[r2_idx];
    uint64_t out_offset = pair_offsets[i];

    if (r1_size == 0 && r2_size == 0) {
        pair_sizes_out[i] = 0;
        return;
    }

    if (r1_size == 0) {
        if (!read_had_mapped_kmers[r1_idx]) {
            const int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];
            for (uint64_t j = 0; j < r2_size; ++j) {
                pair_transcripts[out_offset + j] = r2_list[j];
            }
            pair_sizes_out[i] = r2_size;
        } else {
            pair_sizes_out[i] = 0;
        }
        return;
    }

    if (r2_size == 0) {
        if (!read_had_mapped_kmers[r2_idx]) {
            const int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
            for (uint64_t j = 0; j < r1_size; ++j) {
                pair_transcripts[out_offset + j] = r1_list[j];
            }
            pair_sizes_out[i] = r1_size;
        } else {
            pair_sizes_out[i] = 0;
        }
        return;
    }

    // Both non-empty: sorted set intersection
    const int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
    const int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];

    uint64_t count = 0;
    uint64_t a = 0, b = 0;
    while (a < r1_size && b < r2_size) {
        int va = r1_list[a];
        int vb = r2_list[b];
        if (va == vb) {
            pair_transcripts[out_offset + count] = va;
            ++count;
            ++a;
            ++b;
        } else if (va < vb) {
            ++a;
        } else {
            ++b;
        }
    }

    pair_sizes_out[i] = count;
}