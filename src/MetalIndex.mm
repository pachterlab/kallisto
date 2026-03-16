#import "MetalIndex.h"
#include "BenchmarkStats.h"

#include <algorithm>
#include <iostream>
#include <chrono>
#include <cstdio>
#include <unordered_set>

// ============================================================================
// Hash utilities
// ============================================================================

uint64_t hash_sorted_vector_cpu(const int* vec, size_t n) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < n; ++i) {
        h ^= static_cast<uint64_t>(vec[i]);
        h *= 0x100000001b3ULL;
    }
    if (h == UINT64_MAX) h = UINT64_MAX - 2;
    return h;
}

uint64_t hash_sorted_roaring_cpu(const Roaring& r) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (uint32_t x : r) {
        h ^= static_cast<uint64_t>(x);
        h *= 0x100000001b3ULL;
    }
    if (h == UINT64_MAX) h = UINT64_MAX - 2;
    return h;
}

// ============================================================================
// MetalHashTable — open-addressing hash table with double hashing
// ============================================================================

// KmerHashEntry: (kmer, ec_id, fwd_run, bwd_run)
// fwd_run = how many consecutive k-mers from this position (forward in block) share EC
// bwd_run = how many consecutive k-mers from this position (backward in block) share EC
// Both capped at 255 (uint8_t max).
struct KmerHashEntry {
    uint64_t kmer;
    int32_t  ec_id;
    uint8_t  fwd_run;
    uint8_t  bwd_run;
};

// Helper: collect KmerHashEntry records from a KmerIndex contig traversal
static void collect_kmer_pairs(const KmerIndex& index,
                                std::vector<KmerHashEntry>& out) {
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
            auto ec_it = index.ecmapinv.find(trs);
            int block_kmers = (int)(mc.second - mc.first);  // k-mers in this block
            if (ec_it != index.ecmapinv.end()) {
                KmerIterator kit(blockseq.c_str()), kit_end;
                int pos = 0;
                for (auto it = kit; it != kit_end; ++it, ++pos) {
                    uint8_t fwd = (uint8_t)std::min(block_kmers - pos, 255);
                    uint8_t bwd = (uint8_t)std::min(pos + 1, 255);
                    out.push_back({to_ullong_metal(it->first.rep()),
                                   ec_it->second, fwd, bwd});
                }
            }
            contigpos = mc.second;
            ++j;
        }
    }
}

// Build a bucket hash table: capacity/4 buckets, 4 slots per bucket.
// Each bucket = 64 bytes = one GPU cache line.
//
// Capacity = 2^27 = 134,217,728 total slots → 2^25 = 33,554,432 buckets.
// Load factor ≈ 116M / 134M = 86.5% slots; avg 3.46 kmers per bucket.
// Expected P(bucket overflow) ≈ 26.8%, so avg ~1.27 cache lines per lookup
// vs ~2.3 for per-slot double hashing → ≈1.8× fewer DRAM round-trips.
//
// Hash:  bucket_mask = (capacity/4) - 1 = 2^25 - 1
//        h1  = kmer & bucket_mask
//        h2  = ((kmer >> 17) & bucket_mask) | 1  (odd, coprime with 2^25)
// Probe: bucket_idx = (h1 + probe * h2) & bucket_mask
//        base_slot  = bucket_idx * 4
static void finalize_hash_table(MetalHashTable& t,
                                 std::vector<KmerHashEntry>& entries)
{
    // Deduplicate by kmer key — sort by kmer, keep first occurrence
    // (same kmer always maps to the same ec_id; fwd/bwd_run from first hit).
    // Sort via parallel (kmer, idx) pairs to avoid libc++ issues with custom types.
    {
        std::vector<std::pair<uint64_t, size_t>> keys;
        keys.reserve(entries.size());
        for (size_t i = 0; i < entries.size(); ++i)
            keys.push_back({entries[i].kmer, i});
        std::sort(keys.begin(), keys.end());  // sorts by kmer (first element)

        std::vector<KmerHashEntry> deduped;
        deduped.reserve(entries.size());
        uint64_t prev = ~(uint64_t)0;
        for (auto& kv : keys) {
            if (kv.first == prev) continue;
            prev = kv.first;
            deduped.push_back(entries[kv.second]);
        }
        entries = std::move(deduped);
    }

    Kmer empty_km;
    empty_km.set_deleted();
    t.empty_sentinel = to_ullong_metal(empty_km);
    t.num_kmers      = entries.size();

    // Total slot capacity = 2^27.  Num buckets = capacity/4 = 2^25.
    t.capacity = (size_t)1 << 27;  // 134,217,728 total slots

    // Allocate and initialize table (fill with empty sentinel)
    t.slots.resize(t.capacity);
    KmerSlot* s = t.slots.data();
    for (size_t i = 0; i < t.capacity; ++i) {
        s[i].key     = t.empty_sentinel;
        s[i].value   = -1;
        s[i].fwd_run = 1;
        s[i].bwd_run = 1;
        s[i]._pad    = 0;
    }

    const size_t   num_buckets = t.capacity >> 2;   // capacity / 4
    const uint64_t bucket_mask = (uint64_t)num_buckets - 1;

    for (const auto& e : entries) {
        uint64_t k = e.kmer;
        if (k == t.empty_sentinel) continue;

        uint64_t h1 = k & bucket_mask;
        uint64_t h2 = ((k >> 17) & bucket_mask) | 1;

        bool inserted = false;
        for (size_t probe = 0; !inserted; ++probe) {
            uint64_t bucket = (h1 + probe * h2) & bucket_mask;
            size_t   base   = (size_t)bucket * 4;
            for (int slot = 0; slot < 4; ++slot) {
                if (s[base + slot].key == t.empty_sentinel) {
                    s[base + slot].key     = k;
                    s[base + slot].value   = e.ec_id;
                    s[base + slot].fwd_run = e.fwd_run;
                    s[base + slot].bwd_run = e.bwd_run;
                    inserted = true;
                    break;
                }
                if (s[base + slot].key == k) {
                    inserted = true;  // duplicate
                    break;
                }
            }
        }
    }

    std::cerr << "[metal] bucket hash table: " << num_buckets
              << " buckets, " << t.num_kmers << " kmers, LF="
              << (100.0 * t.num_kmers / t.capacity) << "%\n";
}

void MetalHashTable::build(const KmerIndex& index) {
    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<KmerHashEntry> entries;
    collect_kmer_pairs(index, entries);
    finalize_hash_table(*this, entries);

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.setup_kmer_to_ec_map_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;

    std::cerr << "[metal] kmer table (hash, double-hash): " << num_kmers
              << " kmers, capacity " << capacity
              << ", " << (capacity * 16 / 1024 / 1024) << " MB" << std::endl;
}

void MetalHashTable::build(const GPUIndex& gpu_index) {
    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<KmerHashEntry> entries;
    for (size_t c = 0; c < gpu_index.num_contigs; ++c) {
        const std::string& seq = gpu_index.contig_sequences[c];
        for (const auto& blk : gpu_index.contig_ec_blocks[c]) {
            uint32_t nk = blk.end - blk.start;
            if (nk == 0) continue;
            std::string block_seq = seq.substr(blk.start, nk + gpu_index.k - 1);
            KmerIterator kit(block_seq.c_str()), kit_end;
            int pos = 0;
            for (auto it = kit; it != kit_end; ++it, ++pos) {
                uint8_t fwd = (uint8_t)std::min((int)nk - pos, 255);
                uint8_t bwd = (uint8_t)std::min(pos + 1, 255);
                entries.push_back({to_ullong_metal(it->first.rep()), blk.ec_id, fwd, bwd});
            }
        }
    }
    finalize_hash_table(*this, entries);

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.setup_kmer_to_ec_map_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;

    std::cerr << "[metal] kmer table (hash, gpu_index): " << num_kmers << " kmers" << std::endl;
}

int MetalHashTable::lookup(uint64_t kmer) const {
    if (kmer == empty_sentinel || capacity == 0) return -1;
    const size_t   num_buckets = capacity >> 2;
    const uint64_t bucket_mask = (uint64_t)num_buckets - 1;
    uint64_t h1 = kmer & bucket_mask;
    uint64_t h2 = ((kmer >> 17) & bucket_mask) | 1;
    const KmerSlot* s = slots.data();
    for (size_t probe = 0; probe < num_buckets; ++probe) {
        size_t base = ((size_t)((h1 + probe * h2) & bucket_mask)) * 4;
        for (int slot = 0; slot < 4; ++slot) {
            if (s[base + slot].key == kmer)           return s[base + slot].value;
            if (s[base + slot].key == empty_sentinel) return -1;
        }
    }
    return -1;
}

// ============================================================================
// MetalHashTable — cache save/load
// Format: uint64_t magic, uint64_t capacity, uint64_t num_kmers,
//         uint64_t empty_sentinel, then capacity KmerSlot entries
// ============================================================================

// Magic v4: bucket hash + fwd_run/bwd_run in every slot — invalidates v3 caches.
static const uint64_t HASH_TABLE_CACHE_MAGIC = 0x4B4D42484B54424DULL; // v4

bool MetalHashTable::save_cache(const std::string& path) const {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    uint64_t hdr[4] = { HASH_TABLE_CACHE_MAGIC,
                         (uint64_t)capacity,
                         (uint64_t)num_kmers,
                         empty_sentinel };
    bool ok = (fwrite(hdr, sizeof(uint64_t), 4, f) == 4) &&
              (fwrite(slots.data(), sizeof(KmerSlot), capacity, f) == capacity);
    fclose(f);
    if (!ok) remove(path.c_str());
    return ok;
}

bool MetalHashTable::load_cache(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    uint64_t hdr[4];
    if (fread(hdr, sizeof(uint64_t), 4, f) != 4 || hdr[0] != HASH_TABLE_CACHE_MAGIC) {
        fclose(f); return false;
    }
    capacity       = (size_t)hdr[1];
    num_kmers      = (size_t)hdr[2];
    empty_sentinel = hdr[3];
    slots.resize(capacity);
    bool ok = (fread(slots.data(), sizeof(KmerSlot), capacity, f) == capacity);
    fclose(f);
    if (!ok) { capacity = 0; num_kmers = 0; return false; }
    return true;
}

// ============================================================================
// MetalECMap
// ============================================================================

static void build_ecmap_impl(
    const EcMapInv& ecmapinv,
    MetalBuffer<int>& transcripts,
    MetalBuffer<uint64_t>& offsets_buf,
    size_t& num_ecs)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    int32_t max_ec = -1;
    for (const auto& e : ecmapinv)
        if (e.second > max_ec) max_ec = e.second;

    size_t num_ec = (max_ec >= 0) ? static_cast<size_t>(max_ec + 1) : 0;

    std::vector<const Roaring*> ec_to_trs(num_ec, nullptr);
    for (const auto& e : ecmapinv)
        ec_to_trs[e.second] = &e.first;

    std::vector<int> h_tx;
    std::vector<uint64_t> h_off(num_ec + 1, 0);
    uint64_t off = 0;
    for (size_t ec = 0; ec < num_ec; ++ec) {
        h_off[ec] = off;
        if (ec_to_trs[ec]) {
            for (uint32_t tr : *ec_to_trs[ec]) h_tx.push_back((int)tr);
            off += ec_to_trs[ec]->cardinality();
        }
    }
    h_off[num_ec] = off;
    num_ecs = num_ec;

    transcripts.from_host(h_tx);
    offsets_buf.from_host(h_off);

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.setup_gpu_ecmap_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

void MetalECMap::build(const KmerIndex& index) {
    build_ecmap_impl(index.ecmapinv, transcripts, offsets, num_ecs);
    std::cerr << "[metal] ecmap: " << num_ecs << " ECs, "
              << transcripts.size() << " transcript entries" << std::endl;
}

void MetalECMap::build(const GPUIndex& index) {
    build_ecmap_impl(index.ecmapinv, transcripts, offsets, num_ecs);
    std::cerr << "[metal] ecmap (gpu_index): " << num_ecs << " ECs, "
              << transcripts.size() << " transcript entries" << std::endl;
}

void MetalECMap::append_ecs(const std::vector<int>& new_tx,
                            const std::vector<uint64_t>& new_offsets,
                            int num_new_ecs)
{
    if (num_new_ecs == 0) return;

    size_t old_tx_size  = transcripts.size();
    size_t old_num_ecs  = num_ecs;

    // Grow transcript buffer
    size_t new_tx_size = old_tx_size + new_tx.size();
    transcripts.resize(new_tx_size);
    std::memcpy(transcripts.data() + old_tx_size, new_tx.data(), new_tx.size() * sizeof(int));

    // Grow offsets buffer (shifted by old_tx_size)
    size_t new_num_ecs = old_num_ecs + num_new_ecs;
    offsets.resize(new_num_ecs + 1);
    uint64_t* off_ptr = offsets.data() + old_num_ecs;
    for (int i = 0; i <= num_new_ecs; ++i)
        off_ptr[i] = new_offsets[i] + static_cast<uint64_t>(old_tx_size);

    num_ecs = new_num_ecs;
}

// ============================================================================
// MetalECMapInv
// ============================================================================

void MetalECMapInv::build(const KmerIndex& index) {
    auto t0 = std::chrono::high_resolution_clock::now();

    map.reserve(index.ecmapinv.size() * 2);
    std::unordered_set<uint64_t> seen;
    size_t dups = 0;

    for (const auto& e : index.ecmapinv) {
        uint64_t h = hash_sorted_roaring_cpu(e.first);
        if (!seen.insert(h).second) ++dups;
        map[h] = e.second;
    }
    num_ecs = index.ecmapinv.size();

    if (dups > 0)
        std::cerr << "  WARNING: " << dups << " hash collisions in EC inverse map" << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.setup_gpu_ecmapinv_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

void MetalECMapInv::build(const GPUIndex& index) {
    auto t0 = std::chrono::high_resolution_clock::now();

    map.reserve(index.ecmapinv.size() * 2);
    for (const auto& e : index.ecmapinv) {
        uint64_t h = hash_sorted_roaring_cpu(e.first);
        map[h] = e.second;
    }
    num_ecs = index.ecmapinv.size();

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.setup_gpu_ecmapinv_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

int MetalECMapInv::find(uint64_t hash) const {
    auto it = map.find(hash);
    return (it != map.end()) ? it->second : -1;
}

void MetalECMapInv::insert(uint64_t hash, int ec_id) {
    map[hash] = ec_id;
    ++num_ecs;
}
