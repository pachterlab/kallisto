//
// kallisto.metal
// Apple Metal Shading Language (MSL) port of the kallisto GPU kernels.
// Targets Apple Silicon (M-series) with unified memory.
//

#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

// ============================================================================
// Shared constants
// ============================================================================

#define MAX_ECS_PER_READ 512
#define EMPTY_EC_VALUE   (-1)

// UINT64_MAX is not defined in MSL headers
#define METAL_UINT64_MAX (ulong(-1))

// ============================================================================
// Utility: base encoding
// ============================================================================

static inline uint8_t encode_base_or_invalid(char b) {
    if (b == 'A') return 0;
    if (b == 'C') return 1;
    if (b == 'G') return 2;
    if (b == 'T') return 3;
    return 4;  // invalid
}

// ============================================================================
// Utility: FNV-1a hash over sorted int array
// ============================================================================

static inline uint64_t hash_sorted_vector(const device int* vec, uint64_t size) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (uint64_t i = 0; i < size; ++i) {
        h ^= static_cast<uint64_t>(vec[i]);
        h *= 0x100000001b3ULL;
    }
    if (h == METAL_UINT64_MAX) h = METAL_UINT64_MAX - 2;
    return h;
}

static inline bool verify_transcript_lists_equal(
    const device int* list1, uint64_t size1,
    const device int* list2, uint64_t size2)
{
    if (size1 != size2) return false;
    for (uint64_t i = 0; i < size1; ++i) {
        if (list1[i] != list2[i]) return false;
    }
    return true;
}

// ============================================================================
// kmer_kernel — extract canonical k-mers from reads
// ============================================================================

kernel void kmer_kernel(
    const device char*     reads            [[ buffer(0) ]],
    const device uint64_t* offsets          [[ buffer(1) ]],
    const device uint32_t* lengths          [[ buffer(2) ]],
    const device uint32_t* kmers_counts     [[ buffer(3) ]],
    const device uint64_t* kmers_offsets    [[ buffer(4) ]],
          device uint64_t* out_kmers        [[ buffer(5) ]],
    constant uint64_t&     empty_kmer_value [[ buffer(6) ]],
    constant uint32_t&     num_reads        [[ buffer(7) ]],
    constant uint32_t&     k                [[ buffer(8) ]],
    uint                   gid              [[ thread_position_in_grid ]])
{
    uint32_t r = gid;
    if (r >= num_reads) return;

    const device char* read = reads + offsets[r];
    uint32_t L       = lengths[r];
    uint32_t n_kmers = kmers_counts[r];
    uint64_t out0    = kmers_offsets[r];

    if (L < k || n_kmers == 0) return;

    const uint32_t fwd_new_bit_pos = 2 * (32 - k);
    const uint32_t rev_high_bit_pos = 62;
    const uint64_t low_bits_mask = ~((uint64_t(1) << fwd_new_bit_pos) - 1);

    uint64_t fwd = 0, rev = 0;
    int bad = 0;

    for (uint32_t i = 0; i < k; ++i) {
        char c = read[i];
        uint8_t code = encode_base_or_invalid(c);
        if (code > 3) { ++bad; code = 0; }
        uint32_t fwd_bit_pos = 2 * (31 - i);
        fwd |= (uint64_t(code)) << fwd_bit_pos;
        uint8_t comp = code ^ 0x3;
        uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
        rev |= (uint64_t(comp)) << rev_bit_pos;
    }

    uint64_t canon0 = (fwd < rev) ? fwd : rev;
    out_kmers[out0] = (bad > 0) ? empty_kmer_value : canon0;

    for (uint32_t i = 1; i < n_kmers; ++i) {
        char c_out = read[i - 1];
        uint8_t code_out = encode_base_or_invalid(c_out);
        if (code_out > 3) --bad;

        char c_in = read[i + k - 1];
        uint8_t code_in = encode_base_or_invalid(c_in);
        if (code_in > 3) { ++bad; code_in = 0; }

        fwd = (fwd << 2) | (uint64_t(code_in) << fwd_new_bit_pos);
        uint8_t comp_in = code_in ^ 0x3;
        rev = ((rev >> 2) | (uint64_t(comp_in) << rev_high_bit_pos)) & low_bits_mask;

        uint64_t canon = (fwd < rev) ? fwd : rev;
        out_kmers[out0 + i] = (bad > 0) ? empty_kmer_value : canon;
    }
}

// ============================================================================
// contig_kmer_kernel — extract k-mers from contig EC blocks
// ============================================================================

kernel void contig_kmer_kernel(
    const device char*     seq                  [[ buffer(0) ]],
    const device uint64_t* block_seq_offsets    [[ buffer(1) ]],
    const device uint32_t* block_num_kmers      [[ buffer(2) ]],
    const device int*      block_ec_ids         [[ buffer(3) ]],
    const device uint64_t* block_output_offsets [[ buffer(4) ]],
          device uint64_t* out_kmers            [[ buffer(5) ]],
          device int*      out_ecs              [[ buffer(6) ]],
    constant uint64_t&     empty_kmer_value     [[ buffer(7) ]],
    constant int&          empty_ec_value       [[ buffer(8) ]],
    constant uint32_t&     k                    [[ buffer(9) ]],
    constant uint64_t&     num_blocks           [[ buffer(10) ]],
    uint                   gid                  [[ thread_position_in_grid ]])
{
    uint64_t block_id = gid;
    if (block_id >= num_blocks) return;

    const device char* block_seq = seq + block_seq_offsets[block_id];
    uint32_t n_kmers = block_num_kmers[block_id];
    int ec_id        = block_ec_ids[block_id];
    uint64_t out_base = block_output_offsets[block_id];

    if (n_kmers == 0) return;

    const uint32_t fwd_new_bit_pos  = 2 * (32 - k);
    const uint32_t rev_high_bit_pos = 62;
    const uint64_t low_bits_mask    = ~((uint64_t(1) << fwd_new_bit_pos) - 1);

    uint64_t fwd = 0, rev = 0;
    int bad = 0;

    for (uint32_t i = 0; i < k; ++i) {
        char c = block_seq[i];
        uint8_t code = encode_base_or_invalid(c);
        if (code > 3) { ++bad; code = 0; }
        uint32_t fwd_bit_pos = 2 * (31 - i);
        fwd |= uint64_t(code) << fwd_bit_pos;
        uint8_t comp = code ^ 0x3;
        uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
        rev |= uint64_t(comp) << rev_bit_pos;
    }

    uint64_t canon = (bad == 0) ? ((fwd < rev) ? fwd : rev) : empty_kmer_value;
    out_kmers[out_base] = canon;
    out_ecs[out_base]   = (bad == 0) ? ec_id : empty_ec_value;
    ++out_base;

    for (uint32_t i = 1; i < n_kmers; ++i) {
        char c_out = block_seq[i - 1];
        uint8_t code_out = encode_base_or_invalid(c_out);
        if (code_out > 3) --bad;

        char c_in = block_seq[i + k - 1];
        uint8_t code_in = encode_base_or_invalid(c_in);
        if (code_in > 3) { ++bad; code_in = 0; }

        fwd = (fwd << 2) | (uint64_t(code_in) << fwd_new_bit_pos);
        uint8_t comp_in = code_in ^ 0x3;
        rev = ((rev >> 2) | (uint64_t(comp_in) << rev_high_bit_pos)) & low_bits_mask;

        canon = (bad == 0) ? ((fwd < rev) ? fwd : rev) : empty_kmer_value;
        out_kmers[out_base] = canon;
        out_ecs[out_base]   = (bad == 0) ? ec_id : empty_ec_value;
        ++out_base;
    }
}

// ============================================================================
// kmer_lookup_kernel — bucket hash table lookup (4 slots per bucket/cache line)
//
// Layout: capacity/4 buckets, each with 4 consecutive KmerSlots (64 bytes =
// one GPU cache line).  Each probe loads ONE cache line and checks ALL 4 slots,
// giving 4× better cache-line utilisation than the previous per-slot scheme.
//
// Expected probes per lookup (86.5% load factor, Poisson bucket fill):
//   ~1.27 cache lines  vs  2.3 for the old double-hashing per-slot approach
// → ≈1.8× fewer DRAM round-trips.
//
// Hashing:
//   num_buckets = capacity / 4  (always a power of 2, stored in capacity>>2)
//   h1  = kmer  & bucket_mask           (bucket index, low bits)
//   h2  = ((kmer >> 17) & bucket_mask) | 1  (step, odd → coprime w/ 2^N)
//   base = ((h1 + probe * h2) & bucket_mask) * 4   (slot index in flat array)
// ============================================================================

// 16-byte aligned; must match MetalIndex.h KmerSlot exactly.
// fwd_run: # consecutive k-mers (including this one) sharing the same EC in
//          the forward direction within the de Bruijn block (for SENSE reads).
// bwd_run: same count going backward (for ANTISENSE reads).
// Min 1; capped at 255.  Enables process_reads_kernel to skip whole EC blocks.
struct KmerSlot {
    uint64_t key;      // canonical k-mer; empty_sentinel if empty
    int32_t  value;    // ec_id (or -1)
    uint8_t  fwd_run;  // forward run length
    uint8_t  bwd_run;  // backward run length
    int16_t  _pad;
};

kernel void kmer_lookup_kernel(
    const device KmerSlot* table          [[ buffer(0) ]],
    const device uint64_t* query_kmers    [[ buffer(1) ]],
          device int*      out_ecs        [[ buffer(2) ]],
    constant uint64_t&     capacity       [[ buffer(3) ]],
    constant uint64_t&     empty_sentinel [[ buffer(4) ]],
    constant uint64_t&     num_kmers      [[ buffer(5) ]],
    uint                   gid            [[ thread_position_in_grid ]])
{
    if (gid >= num_kmers) return;

    uint64_t kmer = query_kmers[gid];
    if (kmer == empty_sentinel) {
        out_ecs[gid] = -1;
        return;
    }

    // Bucket double hashing — probe in bucket space (4 slots = 64B per bucket)
    uint64_t num_buckets  = capacity >> 2;   // capacity / 4
    uint64_t bucket_mask  = num_buckets - 1;
    uint64_t h1 = kmer & bucket_mask;
    uint64_t h2 = ((kmer >> 17) & bucket_mask) | 1;

    for (uint64_t probe = 0; probe < num_buckets; ++probe) {
        uint64_t bucket = (h1 + probe * h2) & bucket_mask;
        uint64_t base   = bucket * 4;

        // Check all 4 slots in this cache-line-aligned bucket
        for (uint64_t s = 0; s < 4; ++s) {
            uint64_t slot_key = table[base + s].key;
            if (slot_key == kmer) {
                out_ecs[gid] = table[base + s].value;
                return;
            }
            if (slot_key == empty_sentinel) {
                // Empty slot → kmer not in table
                out_ecs[gid] = -1;
                return;
            }
        }
        // All 4 slots occupied by other kmers → probe next bucket
    }
    out_ecs[gid] = -1;
}

// ============================================================================
// process_reads_kernel — one thread per read
//
// Replaces three separate kernels:
//   kmer_kernel            (k-mer extraction)
//   kmer_lookup_kernel     (hash table lookup)
//   compute_ec_collapse_sizes_kernel  (EC dedup/sort per read)
//
// Design: incremental sliding window (O(1) per position, O(L) total) with
// run-length skipping of hash LOOKUPS only.  After a hit at position p with
// run R, positions p+1 … p+R-1 still update the sliding window state but
// skip the expensive hash probe.  This keeps k-mer extraction at the same
// cost as kmer_kernel while reducing hash table DRAM traffic from ~120
// probes/read to ~2–5 probes/read.
//
// Expected lookups per read: ~2–5 (vs ~120 without skipping).
// Working set per batch: 100K reads × 4 probes × 64B ≈ 25 MB (fits in SLC).
//
// Output: temp_ecs[gid × MAX_PROC_ECS_PER_READ..] + out_ec_counts[gid]
// scatter_proc_ecs_kernel then writes the ragged ReadECCollapser format.
// ============================================================================

#define MAX_PROC_ECS_PER_READ 32
#define TG_KMERS_PER_READ 512

// Bucket hash lookup: returns ec_id (or -1), sets out_fwd/out_bwd run lengths.
static inline int bucket_lookup(
    const device KmerSlot* table,
    uint64_t kmer, uint64_t num_buckets, uint64_t bucket_mask,
    uint64_t empty_sentinel,
    thread uint8_t& out_fwd, thread uint8_t& out_bwd)
{
    uint64_t h1 = kmer & bucket_mask;
    uint64_t h2 = ((kmer >> 17) & bucket_mask) | 1;

    for (uint64_t probe = 0; probe < num_buckets; ++probe) {
        uint64_t base = ((h1 + probe * h2) & bucket_mask) * 4;
        for (uint64_t s = 0; s < 4; ++s) {
            uint64_t k2 = table[base + s].key;
            if (k2 == kmer) {
                out_fwd = table[base + s].fwd_run;
                out_bwd = table[base + s].bwd_run;
                return table[base + s].value;
            }
            if (k2 == empty_sentinel) {
                out_fwd = 1; out_bwd = 1;
                return -1;
            }
        }
    }
    out_fwd = 1; out_bwd = 1;
    return -1;
}

// Insert ec into a sorted, deduplicated private array (binary search).
static inline void insert_ec_sorted_small(
    thread int* ecs, thread int& count, int ec)
{
    int lo = 0, hi = count;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (ecs[mid] < ec) lo = mid + 1;
        else hi = mid;
    }
    if (lo < count && ecs[lo] == ec) return;  // already present
    if (count < MAX_PROC_ECS_PER_READ) {
        for (int j = count; j > lo; --j) ecs[j] = ecs[j-1];
        ecs[lo] = ec;
        ++count;
    }
}

kernel void process_reads_kernel(
    const device char*     reads            [[ buffer(0) ]],
    const device uint64_t* read_offsets     [[ buffer(1) ]],
    const device uint32_t* read_lengths     [[ buffer(2) ]],
    const device KmerSlot* kmer_table       [[ buffer(3) ]],
          device int*      temp_ecs         [[ buffer(4) ]],   // [num_reads × MAX_PROC_ECS_PER_READ]
          device uint64_t* out_ec_counts    [[ buffer(5) ]],
    constant uint64_t&     num_reads        [[ buffer(6) ]],
    constant uint32_t&     k               [[ buffer(7) ]],
    constant uint64_t&     capacity         [[ buffer(8) ]],
    constant uint64_t&     empty_sentinel   [[ buffer(9) ]],
    constant uint64_t&     empty_kmer_value [[ buffer(10) ]],
    uint                   gid             [[ thread_position_in_grid ]])
{
    if (gid >= num_reads) return;

    uint32_t L = read_lengths[gid];
    if (L < k) { out_ec_counts[gid] = 0; return; }

    const device char* read = reads + read_offsets[gid];
    uint32_t n_kmers = L - k + 1;

    uint64_t num_buckets = capacity >> 2;
    uint64_t bucket_mask = num_buckets - 1;

    int unique_ecs[MAX_PROC_ECS_PER_READ];
    int count = 0;

    // Incremental sliding window — same as kmer_kernel but we skip lookups
    // for positions inside a run (next_lookup_pos tracks when to probe next).
    const uint32_t fwd_new_bit_pos  = 2 * (32 - k);
    const uint32_t rev_high_bit_pos = 62;
    const uint64_t low_bits_mask    = ~((uint64_t(1) << fwd_new_bit_pos) - 1);

    uint64_t fwd = 0, rev = 0;
    int bad = 0;

    // Build the initial k-mer at position 0
    for (uint32_t i = 0; i < k; ++i) {
        char c = read[i];
        uint8_t code = encode_base_or_invalid(c);
        if (code > 3) { ++bad; code = 0; }
        uint32_t fwd_bit_pos = 2 * (31 - i);
        fwd |= uint64_t(code) << fwd_bit_pos;
        uint8_t comp = code ^ 0x3;
        uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
        rev |= uint64_t(comp) << rev_bit_pos;
    }

    uint32_t next_lookup_pos = 0;  // first position to actually probe

    for (uint32_t pos = 0; pos < n_kmers; ++pos) {
        // Slide window for positions > 0
        if (pos > 0) {
            char c_out = read[pos - 1];
            uint8_t code_out = encode_base_or_invalid(c_out);
            if (code_out > 3) --bad;

            char c_in = read[pos + k - 1];
            uint8_t code_in = encode_base_or_invalid(c_in);
            if (code_in > 3) { ++bad; code_in = 0; }

            fwd = (fwd << 2) | (uint64_t(code_in) << fwd_new_bit_pos);
            uint8_t comp_in = code_in ^ 0x3;
            rev = ((rev >> 2) | (uint64_t(comp_in) << rev_high_bit_pos)) & low_bits_mask;
        }

        // Skip lookup if inside a run or if k-mer has bad bases
        if (pos < next_lookup_pos || bad > 0) continue;

        bool is_sense  = (fwd <= rev);
        uint64_t canon = is_sense ? fwd : rev;
        if (canon == empty_kmer_value) continue;

        uint8_t fwd_run, bwd_run;
        int ec_id = bucket_lookup(kmer_table, canon, num_buckets, bucket_mask,
                                  empty_sentinel, fwd_run, bwd_run);

        if (ec_id >= 0) {
            insert_ec_sorted_small(unique_ecs, count, ec_id);
            uint32_t run = is_sense ? uint32_t(fwd_run) : uint32_t(bwd_run);
            next_lookup_pos = pos + (run > 0 ? run : 1);
        }
        // else next_lookup_pos stays at pos (increments naturally with loop)
    }

    out_ec_counts[gid] = uint64_t(count);
    device int* my_out = temp_ecs + uint64_t(gid) * MAX_PROC_ECS_PER_READ;
    for (int i = 0; i < count; ++i) my_out[i] = unique_ecs[i];
}

// ============================================================================
// process_reads_threadgroup_kernel — one threadgroup per read
//
// Each lane handles a subset of k-mer positions within the read, so the random
// k-mer hash lookups happen in parallel across the group. Lane 0 then performs
// the small per-read EC dedup/writeback.
// ============================================================================

kernel void process_reads_threadgroup_kernel(
    const device char*     reads            [[ buffer(0) ]],
    const device uint64_t* read_offsets     [[ buffer(1) ]],
    const device uint32_t* read_lengths     [[ buffer(2) ]],
    const device KmerSlot* kmer_table       [[ buffer(3) ]],
          device int*      temp_ecs         [[ buffer(4) ]],
          device uint64_t* out_ec_counts    [[ buffer(5) ]],
    constant uint64_t&     num_reads        [[ buffer(6) ]],
    constant uint32_t&     k                [[ buffer(7) ]],
    constant uint64_t&     capacity         [[ buffer(8) ]],
    constant uint64_t&     empty_sentinel   [[ buffer(9) ]],
    constant uint64_t&     empty_kmer_value [[ buffer(10) ]],
    uint                   tid              [[ thread_index_in_threadgroup ]],
    uint                   group_id         [[ threadgroup_position_in_grid ]],
    uint                   group_size       [[ threads_per_threadgroup ]])
{
    uint64_t read_id = group_id;
    if (read_id >= num_reads) return;

    uint32_t L = read_lengths[read_id];
    if (L < k) {
        if (tid == 0) out_ec_counts[read_id] = 0;
        return;
    }

    const device char* read = reads + read_offsets[read_id];
    uint32_t n_kmers = L - k + 1;
    uint64_t num_buckets = capacity >> 2;
    uint64_t bucket_mask = num_buckets - 1;

    threadgroup int tg_ecs[TG_KMERS_PER_READ];

    for (uint32_t pos = tid; pos < min(n_kmers, (uint32_t)TG_KMERS_PER_READ); pos += group_size) {
        uint64_t fwd = 0, rev = 0;
        int bad = 0;

        for (uint32_t i = 0; i < k; ++i) {
            char c = read[pos + i];
            uint8_t code = encode_base_or_invalid(c);
            if (code > 3) { ++bad; code = 0; }
            uint32_t fwd_bit_pos = 2 * (31 - i);
            fwd |= uint64_t(code) << fwd_bit_pos;
            uint8_t comp = code ^ 0x3;
            uint32_t rev_bit_pos = 2 * (31 - (k - 1 - i));
            rev |= uint64_t(comp) << rev_bit_pos;
        }

        int ec_id = -1;
        if (bad == 0) {
            uint64_t canon = (fwd <= rev) ? fwd : rev;
            if (canon != empty_kmer_value) {
                uint8_t fwd_run, bwd_run;
                ec_id = bucket_lookup(kmer_table, canon, num_buckets, bucket_mask,
                                      empty_sentinel, fwd_run, bwd_run);
            }
        }
        tg_ecs[pos] = ec_id;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid == 0) {
        int unique_ecs[MAX_PROC_ECS_PER_READ];
        int count = 0;
        uint32_t limit = min(n_kmers, (uint32_t)TG_KMERS_PER_READ);
        for (uint32_t pos = 0; pos < limit; ++pos) {
            int ec = tg_ecs[pos];
            if (ec >= 0)
                insert_ec_sorted_small(unique_ecs, count, ec);
        }

        out_ec_counts[read_id] = uint64_t(count);
        device int* my_out = temp_ecs + uint64_t(read_id) * MAX_PROC_ECS_PER_READ;
        for (int i = 0; i < count; ++i) my_out[i] = unique_ecs[i];
    }
}

// ============================================================================
// scatter_proc_ecs_kernel — scatter fixed-stride temp_ecs into ragged read_ecs
// Stride is MAX_PROC_ECS_PER_READ = 32 (must match process_reads_kernel above)
// ============================================================================

kernel void scatter_proc_ecs_kernel(
    const device int*      temp_ecs       [[ buffer(0) ]],
    const device uint64_t* output_offsets [[ buffer(1) ]],
          device int*      output_ecs     [[ buffer(2) ]],
    constant uint64_t&     num_reads      [[ buffer(3) ]],
    uint                   gid            [[ thread_position_in_grid ]])
{
    if (gid >= num_reads) return;
    uint64_t start = output_offsets[gid];
    uint64_t count = output_offsets[gid + 1] - start;
    const device int* src = temp_ecs + uint64_t(gid) * MAX_PROC_ECS_PER_READ;
    for (uint64_t i = 0; i < count; ++i)
        output_ecs[start + i] = src[i];
}

// ============================================================================
// kmer_bsearch_kernel — binary search on sorted kmer table (replaces open-
// addressing hash lookup; SortedKmerTable has 1.39 GB vs 4.3 GB for hash).
// One thread per query kmer.
// ============================================================================

kernel void kmer_bsearch_kernel(
    const device uint64_t* keys           [[ buffer(0) ]],  // sorted kmer values
    const device int32_t*  values         [[ buffer(1) ]],  // parallel EC IDs
    const device uint64_t* query_kmers    [[ buffer(2) ]],  // queries
          device int*      out_ecs        [[ buffer(3) ]],  // output
    constant uint64_t&     table_size     [[ buffer(4) ]],  // # entries in sorted table
    constant uint64_t&     empty_sentinel [[ buffer(5) ]],  // deleted-kmer sentinel
    constant uint64_t&     num_queries    [[ buffer(6) ]],
    uint                   gid            [[ thread_position_in_grid ]])
{
    if (gid >= num_queries) return;

    uint64_t kmer = query_kmers[gid];
    if (kmer == empty_sentinel || table_size == 0) {
        out_ecs[gid] = -1;
        return;
    }

    uint64_t lo = 0, hi = table_size;
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (keys[mid] < kmer) lo = mid + 1;
        else                  hi = mid;
    }
    out_ecs[gid] = (lo < table_size && keys[lo] == kmer) ? values[lo] : -1;
}

// ============================================================================
// compute_ec_collapse_sizes_kernel
// ============================================================================

kernel void compute_ec_collapse_sizes_kernel(
    const device uint64_t* read_id_to_kmer_first [[ buffer(0) ]],
    const device uint32_t* read_id_to_kmer_count [[ buffer(1) ]],
    const device int*      ecs                   [[ buffer(2) ]],
          device uint64_t* output_sizes          [[ buffer(3) ]],
          device int*      temp_ecs              [[ buffer(4) ]],
    constant uint64_t&     num_reads             [[ buffer(5) ]],
    uint                   gid                   [[ thread_position_in_grid ]])
{
    uint64_t read_id = gid;
    if (read_id >= num_reads) return;

    uint64_t kmer_start = read_id_to_kmer_first[read_id];
    uint32_t kmer_count = read_id_to_kmer_count[read_id];

    if (kmer_count == 0) {
        output_sizes[read_id] = 0;
        return;
    }

    int unique_ecs[MAX_ECS_PER_READ];
    int count  = 0;
    int prev_ec = -1;

    for (uint32_t i = 0; i < kmer_count; ++i) {
        int ec = ecs[kmer_start + i];
        if (ec == -1) continue;
        if (ec == prev_ec) continue;
        prev_ec = ec;

        if (count < 8) {
            bool found = false;
            int insert_pos = count;
            for (int j = 0; j < count; ++j) {
                if (unique_ecs[j] == ec) { found = true; break; }
                if (unique_ecs[j] > ec && insert_pos == count) insert_pos = j;
            }
            if (found) continue;
            if (count < MAX_ECS_PER_READ) {
                for (int j = count; j > insert_pos; --j) unique_ecs[j] = unique_ecs[j-1];
                unique_ecs[insert_pos] = ec;
                count++;
            }
        } else {
            int lo = 0, hi = count;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (unique_ecs[mid] < ec) lo = mid + 1;
                else hi = mid;
            }
            if (lo < count && unique_ecs[lo] == ec) continue;
            if (count < MAX_ECS_PER_READ) {
                for (int j = count; j > lo; --j) unique_ecs[j] = unique_ecs[j-1];
                unique_ecs[lo] = ec;
                count++;
            }
        }
    }

    output_sizes[read_id] = count;

    device int* my_temp = temp_ecs + read_id * MAX_ECS_PER_READ;
    for (int i = 0; i < count; ++i) my_temp[i] = unique_ecs[i];
}

// ============================================================================
// collapse_ecs_per_read_kernel
// ============================================================================

kernel void collapse_ecs_per_read_kernel(
    const device int*      temp_ecs        [[ buffer(0) ]],
    const device uint64_t* output_offsets  [[ buffer(1) ]],
          device int*      output_ecs      [[ buffer(2) ]],
    constant uint64_t&     num_reads       [[ buffer(3) ]],
    uint                   gid             [[ thread_position_in_grid ]])
{
    uint64_t read_id = gid;
    if (read_id >= num_reads) return;

    uint64_t output_start = output_offsets[read_id];
    uint64_t output_size  = output_offsets[read_id + 1] - output_start;
    if (output_size == 0) return;

    const device int* my_temp = temp_ecs + read_id * MAX_ECS_PER_READ;
    for (uint64_t i = 0; i < output_size; ++i)
        output_ecs[output_start + i] = my_temp[i];
}

// ============================================================================
// compute_intersection_sizes_kernel
// ============================================================================

kernel void compute_intersection_sizes_kernel(
    const device int*      read_ecs              [[ buffer(0) ]],
    const device uint64_t* read_ec_offsets       [[ buffer(1) ]],
    const device uint64_t* ecmap_offsets         [[ buffer(2) ]],
          device uint64_t* output_sizes          [[ buffer(3) ]],
          device uint64_t* smallest_ec_indices   [[ buffer(4) ]],
    constant uint64_t&     num_reads             [[ buffer(5) ]],
    constant uint64_t&     num_ecs               [[ buffer(6) ]],
    constant uint64_t&     max_transcripts       [[ buffer(7) ]],
    constant uint64_t&     read_ecs_size         [[ buffer(8) ]],
    uint                   gid                   [[ thread_position_in_grid ]])
{
    uint64_t read_id = gid;
    if (read_id >= num_reads) return;

    output_sizes[read_id] = 0;
    smallest_ec_indices[read_id] = 0;

    uint64_t ec_start = read_ec_offsets[read_id];
    uint64_t ec_end   = read_ec_offsets[read_id + 1];

    if (ec_end <= ec_start || ec_start >= read_ecs_size) return;

    uint64_t ec_count = ec_end - ec_start;
    int first_ec = read_ecs[ec_start];
    uint64_t first_tx_start = ecmap_offsets[first_ec];
    uint64_t first_tx_end   = ecmap_offsets[first_ec + 1];

    if (ec_count == 1) {
        output_sizes[read_id]        = first_tx_end - first_tx_start;
        smallest_ec_indices[read_id] = ec_start;
        return;
    }

    uint64_t min_ec_size = first_tx_end - first_tx_start;
    uint64_t min_ec_idx  = ec_start;

    for (uint64_t i = ec_start + 1; i < ec_end; ++i) {
        int ec = read_ecs[i];
        uint64_t tx_start = ecmap_offsets[ec];
        uint64_t tx_end   = ecmap_offsets[ec + 1];
        uint64_t ec_size  = tx_end - tx_start;
        if (ec_size < min_ec_size) {
            min_ec_size = ec_size;
            min_ec_idx  = i;
        }
    }

    output_sizes[read_id]        = min_ec_size;
    smallest_ec_indices[read_id] = min_ec_idx;
}

// ============================================================================
// intersect_transcripts_kernel
// ============================================================================

kernel void intersect_transcripts_kernel(
    const device int*      read_ecs              [[ buffer(0) ]],
    const device uint64_t* read_ec_offsets       [[ buffer(1) ]],
    const device int*      ecmap_transcripts     [[ buffer(2) ]],
    const device uint64_t* ecmap_offsets         [[ buffer(3) ]],
    const device uint64_t* output_offsets        [[ buffer(4) ]],
          device int*      output_transcripts    [[ buffer(5) ]],
          device uint64_t* output_sizes          [[ buffer(6) ]],
    const device uint64_t* smallest_ec_indices   [[ buffer(7) ]],
    constant uint64_t&     num_reads             [[ buffer(8) ]],
    constant uint64_t&     num_ecs               [[ buffer(9) ]],
    constant uint64_t&     max_transcripts       [[ buffer(10) ]],
    constant uint64_t&     read_ecs_size         [[ buffer(11) ]],
    constant uint64_t&     max_output_transcripts[[ buffer(12) ]],
    uint                   gid                   [[ thread_position_in_grid ]])
{
    uint64_t read_id = gid;
    if (read_id >= num_reads) return;

    uint64_t ec_start = read_ec_offsets[read_id];
    uint64_t ec_end   = read_ec_offsets[read_id + 1];

    if (ec_end <= ec_start || ec_start >= read_ecs_size) {
        output_sizes[read_id] = 0;
        return;
    }

    uint64_t ec_count    = ec_end - ec_start;
    uint64_t output_start = output_offsets[read_id];
    uint64_t output_end   = output_offsets[read_id + 1];

    if (output_end <= output_start) { output_sizes[read_id] = 0; return; }

    uint64_t smallest_ec_idx = smallest_ec_indices[read_id];
    int      smallest_ec     = read_ecs[smallest_ec_idx];
    uint64_t smallest_tx_start = ecmap_offsets[smallest_ec];
    uint64_t smallest_tx_end   = ecmap_offsets[smallest_ec + 1];
    uint64_t smallest_ec_size  = smallest_tx_end - smallest_tx_start;

    uint64_t output_limit = min(output_end, output_start + max_output_transcripts);

    if (ec_count == 1) {
        uint64_t idx = 0;
        for (uint64_t i = smallest_tx_start; i < smallest_tx_end && output_start + idx < output_limit; ++i, ++idx)
            output_transcripts[output_start + idx] = ecmap_transcripts[i];
        output_sizes[read_id] = idx;
        return;
    }

    uint64_t current_size = smallest_ec_size;
    uint64_t copy_count = min(current_size, output_limit - output_start);
    for (uint64_t i = 0; i < copy_count; ++i)
        output_transcripts[output_start + i] = ecmap_transcripts[smallest_tx_start + i];
    current_size = copy_count;

    if (current_size == 0) { output_sizes[read_id] = 0; return; }

    for (uint64_t ec_idx = ec_start; ec_idx < ec_end && current_size > 0; ++ec_idx) {
        if (ec_idx == smallest_ec_idx) continue;

        int ec = read_ecs[ec_idx];
        uint64_t tx_start = ecmap_offsets[ec];
        uint64_t tx_end   = ecmap_offsets[ec + 1];
        uint64_t ec_size  = tx_end - tx_start;
        uint64_t write_idx = 0;

        int first_current = output_transcripts[output_start];
        int last_current  = output_transcripts[output_start + current_size - 1];
        int first_ec_tx   = ecmap_transcripts[tx_start];
        int last_ec_tx    = ecmap_transcripts[tx_end - 1];

        if (last_current < first_ec_tx || last_ec_tx < first_current) {
            current_size = 0;
            break;
        }

        if (ec_size > current_size * 5) {
            for (uint64_t i = 0; i < current_size && output_start + write_idx < output_limit; ++i) {
                int val = output_transcripts[output_start + i];
                int64_t lo = 0, hi = int64_t(ec_size) - 1;
                bool found = false;
                while (lo <= hi) {
                    int64_t mid = lo + (hi - lo) / 2;
                    int b = ecmap_transcripts[tx_start + mid];
                    if (b == val) { found = true; break; }
                    else if (b < val) lo = mid + 1;
                    else hi = mid - 1;
                }
                if (found) output_transcripts[output_start + write_idx++] = val;
            }
        } else {
            uint64_t read_i = 0, j = tx_start;
            while (read_i < current_size && j < tx_end && output_start + write_idx < output_limit) {
                int a = output_transcripts[output_start + read_i];
                int b = ecmap_transcripts[j];
                if (a < b) ++read_i;
                else if (b < a) ++j;
                else {
                    output_transcripts[output_start + write_idx] = a;
                    ++write_idx; ++read_i; ++j;
                }
            }
        }

        current_size = write_idx;
        if (current_size == 0) break;
    }

    output_sizes[read_id] = current_size;
}

// ============================================================================
// compute_pair_intersection_sizes_kernel
// ============================================================================

kernel void compute_pair_intersection_sizes_kernel(
    const device int*      read_transcripts        [[ buffer(0) ]],
    const device uint64_t* read_transcript_offsets [[ buffer(1) ]],
    const device uint64_t* read_transcript_sizes   [[ buffer(2) ]],
    const device uint8_t*  read_had_mapped_kmers   [[ buffer(3) ]],
          device uint64_t* pair_sizes              [[ buffer(4) ]],
    constant uint32_t&     r1_count                [[ buffer(5) ]],
    uint                   gid                     [[ thread_position_in_grid ]])
{
    uint32_t i = gid;
    if (i >= r1_count) return;

    uint32_t r1_idx = i;
    uint32_t r2_idx = r1_count + i;

    uint64_t r1_size = read_transcript_sizes[r1_idx];
    uint64_t r2_size = read_transcript_sizes[r2_idx];

    if (r1_size == 0 && r2_size == 0) { pair_sizes[i] = 0; return; }
    if (r1_size == 0) {
        pair_sizes[i] = read_had_mapped_kmers[r1_idx] ? 0 : r2_size;
        return;
    }
    if (r2_size == 0) {
        pair_sizes[i] = read_had_mapped_kmers[r2_idx] ? 0 : r1_size;
        return;
    }

    const device int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
    const device int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];

    uint64_t count = 0, a = 0, b = 0;
    while (a < r1_size && b < r2_size) {
        int va = r1_list[a], vb = r2_list[b];
        if (va == vb) { ++count; ++a; ++b; }
        else if (va < vb) ++a;
        else ++b;
    }
    pair_sizes[i] = count;
}

// ============================================================================
// intersect_pairs_kernel
// ============================================================================

kernel void intersect_pairs_kernel(
    const device int*      read_transcripts        [[ buffer(0) ]],
    const device uint64_t* read_transcript_offsets [[ buffer(1) ]],
    const device uint64_t* read_transcript_sizes   [[ buffer(2) ]],
    const device uint8_t*  read_had_mapped_kmers   [[ buffer(3) ]],
    const device uint64_t* pair_offsets            [[ buffer(4) ]],
          device int*      pair_transcripts        [[ buffer(5) ]],
          device uint64_t* pair_sizes_out          [[ buffer(6) ]],
    constant uint32_t&     r1_count                [[ buffer(7) ]],
    uint                   gid                     [[ thread_position_in_grid ]])
{
    uint32_t i = gid;
    if (i >= r1_count) return;

    uint32_t r1_idx = i, r2_idx = r1_count + i;
    uint64_t r1_size = read_transcript_sizes[r1_idx];
    uint64_t r2_size = read_transcript_sizes[r2_idx];
    uint64_t out_offset = pair_offsets[i];

    if (r1_size == 0 && r2_size == 0) { pair_sizes_out[i] = 0; return; }

    if (r1_size == 0) {
        if (!read_had_mapped_kmers[r1_idx]) {
            const device int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];
            for (uint64_t j = 0; j < r2_size; ++j) pair_transcripts[out_offset + j] = r2_list[j];
            pair_sizes_out[i] = r2_size;
        } else { pair_sizes_out[i] = 0; }
        return;
    }
    if (r2_size == 0) {
        if (!read_had_mapped_kmers[r2_idx]) {
            const device int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
            for (uint64_t j = 0; j < r1_size; ++j) pair_transcripts[out_offset + j] = r1_list[j];
            pair_sizes_out[i] = r1_size;
        } else { pair_sizes_out[i] = 0; }
        return;
    }

    const device int* r1_list = read_transcripts + read_transcript_offsets[r1_idx];
    const device int* r2_list = read_transcripts + read_transcript_offsets[r2_idx];

    uint64_t count = 0, a = 0, b = 0;
    while (a < r1_size && b < r2_size) {
        int va = r1_list[a], vb = r2_list[b];
        if (va == vb) { pair_transcripts[out_offset + count] = va; ++count; ++a; ++b; }
        else if (va < vb) ++a;
        else ++b;
    }
    pair_sizes_out[i] = count;
}

// ============================================================================
// hash_transcript_sets_kernel — FNV-1a per transcript set
// ============================================================================

kernel void hash_transcript_sets_kernel(
    const device int*      transcripts  [[ buffer(0) ]],
    const device uint64_t* offsets      [[ buffer(1) ]],
    const device uint64_t* sizes        [[ buffer(2) ]],
          device uint64_t* out_hashes   [[ buffer(3) ]],
    constant uint64_t&     max_tx       [[ buffer(4) ]],
    constant uint64_t&     count        [[ buffer(5) ]],
    uint                   gid          [[ thread_position_in_grid ]])
{
    if (gid >= count) return;

    uint64_t tx_start = offsets[gid];
    uint64_t tx_size  = sizes[gid];

    if (tx_size == 0 || tx_start + tx_size > max_tx) {
        out_hashes[gid] = METAL_UINT64_MAX;
        return;
    }

    out_hashes[gid] = hash_sorted_vector(transcripts + tx_start, tx_size);
}

// ============================================================================
// compute_denom_kernel — EM step 1: denominator per EC
// ============================================================================

kernel void compute_denom_kernel(
    const device float*    alpha      [[ buffer(0) ]],
    const device float*    eff_lens   [[ buffer(1) ]],
    const device int*      transcripts[[ buffer(2) ]],
    const device uint64_t* offsets    [[ buffer(3) ]],
    const device int*      ec_counts  [[ buffer(4) ]],
          device float*    denom      [[ buffer(5) ]],
    constant uint64_t&     num_ecs    [[ buffer(6) ]],
    uint                   gid        [[ thread_position_in_grid ]])
{
    size_t e = gid;
    if (e >= num_ecs) return;

    if (ec_counts[e] == 0) { denom[e] = 0.0f; return; }

    uint64_t start = offsets[e];
    uint64_t end   = offsets[e + 1];

    float d = 0.0f;
    for (uint64_t i = start; i < end; ++i) {
        int t = transcripts[i];
        d += alpha[t] / eff_lens[t];
    }
    denom[e] = d;
}

// ============================================================================
// em_gather_kernel — EM step 2: one thread per transcript (no atomicAdd)
// ============================================================================

kernel void em_gather_kernel(
    const device float*    alpha            [[ buffer(0) ]],
    const device float*    eff_lens         [[ buffer(1) ]],
    const device int*      trans_ecs        [[ buffer(2) ]],
    const device int*      trans_ec_offsets [[ buffer(3) ]],
    const device int*      ec_counts        [[ buffer(4) ]],
    const device float*    denom            [[ buffer(5) ]],
          device float*    alpha_next       [[ buffer(6) ]],
    constant uint64_t&     num_trans        [[ buffer(7) ]],
    constant float&        tolerance        [[ buffer(8) ]],
    uint                   gid              [[ thread_position_in_grid ]])
{
    size_t t = gid;
    if (t >= num_trans) return;

    int start = trans_ec_offsets[t];
    int end   = trans_ec_offsets[t + 1];

    if (start >= end) { alpha_next[t] = 0.0f; return; }

    float a_over_el = alpha[t] / eff_lens[t];
    float sum = 0.0f;

    for (int i = start; i < end; ++i) {
        int e     = trans_ecs[i];
        int count = ec_counts[e];
        if (count == 0) continue;
        float d = denom[e];
        if (d < tolerance) continue;
        sum += float(count) * a_over_el / d;
    }

    alpha_next[t] = sum;
}

// ============================================================================
// scatter_ec_counts_kernel — atomic scatter of EC id counts
// ============================================================================

kernel void scatter_ec_counts_kernel(
    const device int*    ec_ids     [[ buffer(0) ]],
          device atomic_int* ec_counts [[ buffer(1) ]],
    constant uint64_t&   num_ids    [[ buffer(2) ]],
    constant uint64_t&   num_ecs    [[ buffer(3) ]],
    uint                 gid        [[ thread_position_in_grid ]])
{
    if (gid >= num_ids) return;
    int ec_id = ec_ids[gid];
    if (ec_id >= 0 && uint64_t(ec_id) < num_ecs)
        atomic_fetch_add_explicit(ec_counts + ec_id, 1, memory_order_relaxed);
}

// ============================================================================
// Exclusive prefix scan for uint64_t arrays (Blelloch, threadgroup 1024)
// Two-phase: per-block scan, then add block offsets
// ============================================================================

kernel void prefix_scan_exclusive_u64(
    const device uint64_t* input        [[ buffer(0) ]],
          device uint64_t* output       [[ buffer(1) ]],
          device uint64_t* block_sums   [[ buffer(2) ]],
    constant uint64_t&     n            [[ buffer(3) ]],
    threadgroup uint64_t*  shared       [[ threadgroup(0) ]],
    uint                   tid          [[ thread_index_in_threadgroup ]],
    uint                   gid          [[ thread_position_in_grid ]],
    uint                   group_id     [[ threadgroup_position_in_grid ]],
    uint                   group_size   [[ threads_per_threadgroup ]])
{
    // Load
    shared[tid] = (gid < n) ? input[gid] : 0;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Up-sweep (reduce)
    for (uint stride = 1; stride < group_size; stride <<= 1) {
        uint idx = (tid + 1) * (stride * 2) - 1;
        if (idx < group_size)
            shared[idx] += shared[idx - stride];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    // Save block sum and clear last
    if (tid == group_size - 1) {
        block_sums[group_id] = shared[tid];
        shared[tid] = 0;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Down-sweep
    for (uint stride = group_size >> 1; stride > 0; stride >>= 1) {
        uint idx = (tid + 1) * (stride * 2) - 1;
        if (idx < group_size) {
            uint64_t t_val = shared[idx - stride];
            shared[idx - stride] = shared[idx];
            shared[idx] += t_val;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (gid < n) output[gid] = shared[tid];
}

kernel void prefix_scan_add_block_offsets_u64(
          device uint64_t* data         [[ buffer(0) ]],
    const device uint64_t* block_sums   [[ buffer(1) ]],
    constant uint64_t&     n            [[ buffer(2) ]],
    uint                   gid          [[ thread_position_in_grid ]],
    uint                   group_id     [[ threadgroup_position_in_grid ]])
{
    if (gid >= n) return;
    // group 0 gets offset 0; groups 1+ add block_sums[group_id]
    // block_sums itself must have been prefix-scanned beforehand
    data[gid] += block_sums[group_id];
}
