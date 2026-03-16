#import "MetalPipeline.h"

#include <cstring>
#include <fstream>
#include <stdexcept>

// ============================================================================
// Exclusive prefix scan over uint64_t array (CPU fallback for small arrays,
// GPU two-phase Blelloch scan for larger ones)
// ============================================================================

// CPU exclusive prefix scan (fast for small counts)
static void exclusive_scan_u64_cpu(const uint64_t* in, uint64_t* out, size_t n) {
    uint64_t running = 0;
    for (size_t i = 0; i < n; ++i) {
        out[i] = running;
        running += in[i];
    }
}

// Compute total from prefix-scanned output array + last element of input
static uint64_t prefix_scan_total(const uint64_t* in_sizes,
                                  const uint64_t* out_offsets,
                                  size_t n)
{
    if (n == 0) return 0;
    return out_offsets[n - 1] + in_sizes[n - 1];
}

// ============================================================================
// HostReadLoader
// ============================================================================

HostReadLoader::HostReadLoader(const ProgramOptions& opt,
                               uint32_t pair_limit_,
                               uint32_t read_length_hint)
    : pair_limit(pair_limit_),
      is_paired(!opt.single_end && opt.files.size() >= 2)
{
    uint32_t max_reads = is_paired ? (pair_limit * 2) : pair_limit;
    read_id_to_offset.resize(max_reads);
    read_id_to_length.resize(max_reads);
    read_id_to_kmer_first.resize(max_reads);
    read_id_to_kmer_count.resize(max_reads);
    reads.reserve(static_cast<size_t>(max_reads) * read_length_hint);
    // fetch_buf is the buffer passed to fetchSequences. It fills until
    // bufpos would exceed limit. We use 8 MB (same as MasterProcessor bufsize)
    // so each fetchSequences call returns ~40K pairs. The carryover logic in
    // load() handles the case where a batch needs multiple fetchSequences calls.
    fetch_buf.resize(1ULL << 23);  // 8 MB
}

bool HostReadLoader::load(FastqSequenceReader& reader) {
    auto t0 = std::chrono::high_resolution_clock::now();

    reads.clear();
    read_count = 0; r1_count = 0; kmer_count = 0;

    // If all pending sequences are consumed and no more data is coming, we're done.
    if (pending_offset >= pending_seqs.size() && !pending_has_more) return false;

    struct ReadEntry { uint32_t length; uint64_t temp_offset; };
    std::vector<ReadEntry> r1_entries, r2_entries;
    r1_entries.reserve(pair_limit);
    if (is_paired) r2_entries.reserve(pair_limit);

    // names/quals/umis/flags are only needed for the fetchSequences call itself
    std::vector<std::pair<const char*, int>> names, quals;
    std::vector<uint32_t> flags;
    std::vector<std::string> umis;

    uint32_t pairs_loaded = 0;

    while (pairs_loaded < pair_limit) {
        // Refill pending_seqs when exhausted
        if (pending_offset >= pending_seqs.size()) {
            if (!pending_has_more) break;
            pending_has_more = reader.fetchSequences(
                fetch_buf.data(), (int)fetch_buf.size(),
                pending_seqs, names, quals, flags, umis, readbatch_id,
                false, false);
            pending_offset = 0;
            if (pending_seqs.empty()) {
                if (!pending_has_more) break;
                continue;
            }
        }

        // Consume sequences from pending_seqs[pending_offset ...]
        if (is_paired) {
            while (pending_offset + 1 < pending_seqs.size() && pairs_loaded < pair_limit) {
                size_t j = pending_offset;
                uint64_t off1 = reads.size();
                reads.insert(reads.end(), pending_seqs[j].first,
                             pending_seqs[j].first + pending_seqs[j].second);
                r1_entries.push_back({(uint32_t)pending_seqs[j].second, off1});

                uint64_t off2 = reads.size();
                reads.insert(reads.end(), pending_seqs[j+1].first,
                             pending_seqs[j+1].first + pending_seqs[j+1].second);
                r2_entries.push_back({(uint32_t)pending_seqs[j+1].second, off2});

                pending_offset += 2;
                ++pairs_loaded;
            }
            // If exactly one orphan sequence remains at the end of the batch, skip it.
            // (Shouldn't happen with well-formed paired FASTQ; if pairs_loaded==pair_limit
            // and complete pairs remain, pending_offset is preserved for the next load().)
            if (pending_offset + 1 == pending_seqs.size())
                pending_offset = pending_seqs.size();
        } else {
            while (pending_offset < pending_seqs.size() && pairs_loaded < pair_limit) {
                size_t j = pending_offset;
                uint64_t off1 = reads.size();
                reads.insert(reads.end(), pending_seqs[j].first,
                             pending_seqs[j].first + pending_seqs[j].second);
                r1_entries.push_back({(uint32_t)pending_seqs[j].second, off1});
                ++pending_offset;
                ++pairs_loaded;
            }
        }
    }

    uint32_t n_r1 = (uint32_t)r1_entries.size();
    uint32_t n_r2 = (uint32_t)r2_entries.size();
    uint32_t total = n_r1 + n_r2;
    if (total == 0) return false;

    read_id_to_offset.resize(total);
    read_id_to_length.resize(total);
    read_id_to_kmer_first.resize(total);
    read_id_to_kmer_count.resize(total);

    uint32_t k = Kmer::k;
    uint64_t running_kmer = 0;

    for (uint32_t i = 0; i < n_r1; ++i) {
        read_id_to_offset[i] = r1_entries[i].temp_offset;
        read_id_to_length[i] = r1_entries[i].length;
        uint32_t nk = (r1_entries[i].length >= k) ? (r1_entries[i].length - k + 1) : 0;
        read_id_to_kmer_first[i] = running_kmer;
        read_id_to_kmer_count[i] = nk;
        running_kmer += nk;
    }
    for (uint32_t i = 0; i < n_r2; ++i) {
        uint32_t idx = n_r1 + i;
        read_id_to_offset[idx] = r2_entries[i].temp_offset;
        read_id_to_length[idx] = r2_entries[i].length;
        uint32_t nk = (r2_entries[i].length >= k) ? (r2_entries[i].length - k + 1) : 0;
        read_id_to_kmer_first[idx] = running_kmer;
        read_id_to_kmer_count[idx] = nk;
        running_kmer += nk;
    }

    read_count = total;
    r1_count   = n_r1;
    kmer_count = running_kmer;

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.io_decompress_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;

    return !reads.empty();
}

// ============================================================================
// DeviceKmerLoader
// ============================================================================

void DeviceKmerLoader::load(const HostReadLoader& h) {
    auto t0 = std::chrono::high_resolution_clock::now();

    reads.from_host(h.reads);
    read_id_to_offset.from_host(h.read_id_to_offset);
    read_id_to_length.from_host(h.read_id_to_length);
    read_id_to_kmer_first.from_host(h.read_id_to_kmer_first);
    read_id_to_kmer_count.from_host(h.read_id_to_kmer_count);
    kmers.resize(h.kmer_count);
    read_count = h.read_count;
    r1_count   = h.r1_count;

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.host_h2d_copy_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

void DeviceKmerLoader::run(uint64_t empty_kmer_value, uint32_t k) {
    if (read_count == 0) return;

    MetalContext::get().dispatch(
        "kmer_kernel",
        (NSUInteger)read_count,
        {
            reads.metalBuffer(),
            read_id_to_offset.metalBuffer(),
            read_id_to_length.metalBuffer(),
            read_id_to_kmer_count.metalBuffer(),
            read_id_to_kmer_first.metalBuffer(),
            kmers.metalBuffer()
        },
        {
            as_constant(empty_kmer_value),
            as_constant((uint32_t)read_count),
            as_constant(k)
        }
    );
}

// ============================================================================
// ReadECCollapser
// ============================================================================

void ReadECCollapser::collapse(const DeviceKmerLoader& loader,
                               const MetalBuffer<int>& d_ecs)
{
    read_count = loader.read_count;

    if (read_count == 0) {
        read_ec_offsets.resize(1, 0);
        read_ecs.resize(1);
        return;
    }

    constexpr int MAX_ECS_PER_READ = 512;
    size_t temp_size = read_count * MAX_ECS_PER_READ;
    if (temp_ecs_.size() < temp_size) temp_ecs_.resize(temp_size);

    // Sizes buffer (output of first kernel)
    MetalBuffer<uint64_t> output_sizes(read_count);

    auto t0 = std::chrono::high_resolution_clock::now();

    MetalContext::get().dispatch(
        "compute_ec_collapse_sizes_kernel",
        (NSUInteger)read_count,
        {
            loader.read_id_to_kmer_first.metalBuffer(),
            loader.read_id_to_kmer_count.metalBuffer(),
            d_ecs.metalBuffer(),
            output_sizes.metalBuffer(),
            temp_ecs_.metalBuffer()
        },
        { as_constant((uint64_t)read_count) }
    );

    // CPU prefix scan
    read_ec_offsets.resize(read_count + 1);
    exclusive_scan_u64_cpu(output_sizes.data(), read_ec_offsets.data(), read_count);
    uint64_t total = prefix_scan_total(output_sizes.data(), read_ec_offsets.data(), read_count);
    read_ec_offsets[read_count] = total;

    read_ecs.resize(total == 0 ? 1 : total);

    MetalContext::get().dispatch(
        "collapse_ecs_per_read_kernel",
        (NSUInteger)read_count,
        {
            temp_ecs_.metalBuffer(),
            read_ec_offsets.metalBuffer(),
            read_ecs.metalBuffer()
        },
        { as_constant((uint64_t)read_count) }
    );

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_ec_collapse_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

// ============================================================================
// ReadECCollapser::process_reads — run-length skipping path
//
// Dispatches process_reads_kernel (one thread per read; does k-mer extraction
// + bucket hash lookup + run-length skip + EC dedup/sort all in one pass).
// Then CPU-prefix-scans the per-read counts and scatters to ragged format via
// scatter_proc_ecs_kernel.
// ============================================================================

void ReadECCollapser::process_reads(const DeviceKmerLoader& loader,
                                     const MetalHashTable& kmer_table,
                                     uint32_t k)
{
    read_count = loader.read_count;

    if (read_count == 0) {
        read_ec_offsets.resize(1, 0);
        read_ecs.resize(1);
        return;
    }

    // MAX_PROC_ECS_PER_READ = 32 must match the Metal define.
    constexpr size_t PROC_STRIDE = 32;
    size_t temp_size = (size_t)read_count * PROC_STRIDE;
    if (temp_ecs_.size() < temp_size) temp_ecs_.resize(temp_size);

    MetalBuffer<uint64_t> output_sizes(read_count);

    auto t0 = std::chrono::high_resolution_clock::now();

    uint64_t cap          = (uint64_t)kmer_table.capacity;
    uint64_t empty_sent   = kmer_table.empty_sentinel;
    Kmer empty_km; empty_km.set_empty();
    uint64_t empty_kv     = to_ullong_metal(empty_km);
    uint64_t n_reads      = (uint64_t)read_count;

    MetalContext::get().dispatch(
        "process_reads_kernel",
        (NSUInteger)read_count,
        {
            loader.reads.metalBuffer(),            // buffer(0): flat read chars
            loader.read_id_to_offset.metalBuffer(),// buffer(1): byte offsets
            loader.read_id_to_length.metalBuffer(),// buffer(2): read lengths
            kmer_table.slots.metalBuffer(),        // buffer(3): hash table
            temp_ecs_.metalBuffer(),               // buffer(4): fixed-stride ECs out
            output_sizes.metalBuffer()             // buffer(5): per-read EC count
        },
        {
            as_constant(n_reads),                  // buffer(6)
            as_constant(k),                        // buffer(7)
            as_constant(cap),                      // buffer(8)
            as_constant(empty_sent),               // buffer(9)
            as_constant(empty_kv)                  // buffer(10)
        }
    );

    // CPU prefix scan: counts → offsets
    read_ec_offsets.resize(read_count + 1);
    exclusive_scan_u64_cpu(output_sizes.data(), read_ec_offsets.data(), read_count);
    uint64_t total = prefix_scan_total(output_sizes.data(), read_ec_offsets.data(), read_count);
    read_ec_offsets[read_count] = total;

    read_ecs.resize(total == 0 ? 1 : total);

    // Scatter fixed-stride temp into ragged read_ecs
    MetalContext::get().dispatch(
        "scatter_proc_ecs_kernel",
        (NSUInteger)read_count,
        {
            temp_ecs_.metalBuffer(),               // buffer(0)
            read_ec_offsets.metalBuffer(),         // buffer(1)
            read_ecs.metalBuffer()                 // buffer(2)
        },
        { as_constant(n_reads) }                   // buffer(3)
    );

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_kmer_lookup_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

// ============================================================================
// ReadTranscriptIntersector
// ============================================================================

void ReadTranscriptIntersector::intersect(const ReadECCollapser& collapser,
                                          const MetalECMap& ecmap)
{
    read_count = collapser.read_count;

    if (read_count == 0) {
        read_transcript_offsets.resize(1, 0);
        read_transcript_sizes.resize(1, 0);
        read_transcripts.resize(1);
        read_had_mapped_kmers.resize(1, 0);
        return;
    }

    MetalBuffer<uint64_t> output_sizes(read_count);
    MetalBuffer<uint64_t> smallest_ec_indices(read_count);

    auto t0 = std::chrono::high_resolution_clock::now();

    uint64_t num_reads    = read_count;
    uint64_t num_ecs      = ecmap.num_ecs;
    uint64_t max_tx       = ecmap.transcripts.size();
    uint64_t read_ecs_sz  = collapser.read_ecs.size();

    MetalContext::get().dispatch(
        "compute_intersection_sizes_kernel",
        (NSUInteger)read_count,
        {
            collapser.read_ecs.metalBuffer(),
            collapser.read_ec_offsets.metalBuffer(),
            ecmap.offsets.metalBuffer(),
            output_sizes.metalBuffer(),
            smallest_ec_indices.metalBuffer()
        },
        {
            as_constant(num_reads),
            as_constant(num_ecs),
            as_constant(max_tx),
            as_constant(read_ecs_sz)
        }
    );

    // Prefix scan
    read_transcript_offsets.resize(read_count + 1);
    exclusive_scan_u64_cpu(output_sizes.data(), read_transcript_offsets.data(), read_count);
    uint64_t total = prefix_scan_total(output_sizes.data(), read_transcript_offsets.data(), read_count);
    read_transcript_offsets[read_count] = total;

    read_transcripts.resize(total == 0 ? 1 : total);

    MetalContext::get().dispatch(
        "intersect_transcripts_kernel",
        (NSUInteger)read_count,
        {
            collapser.read_ecs.metalBuffer(),
            collapser.read_ec_offsets.metalBuffer(),
            ecmap.transcripts.metalBuffer(),
            ecmap.offsets.metalBuffer(),
            read_transcript_offsets.metalBuffer(),
            read_transcripts.metalBuffer(),
            output_sizes.metalBuffer(),
            smallest_ec_indices.metalBuffer()
        },
        {
            as_constant(num_reads),
            as_constant(num_ecs),
            as_constant(max_tx),
            as_constant(read_ecs_sz),
            as_constant(total == 0 ? (uint64_t)1 : total)
        }
    );

    read_transcript_sizes = std::move(output_sizes);

    // Build read_had_mapped_kmers: 1 if read had at least one mapped kmer (ec count > 0)
    read_had_mapped_kmers.resize(read_count);
    const uint64_t* ec_off = collapser.read_ec_offsets.data();
    uint8_t* mapped = read_had_mapped_kmers.data();
    for (size_t i = 0; i < read_count; ++i)
        mapped[i] = (ec_off[i + 1] > ec_off[i]) ? 1 : 0;

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_transcript_intersection_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

// ============================================================================
// PairIntersector
// ============================================================================

void PairIntersector::intersect_pairs(const ReadTranscriptIntersector& intersector,
                                      uint32_t r1_count_)
{
    pair_count = r1_count_;

    if (pair_count == 0) {
        pair_transcript_offsets.resize(1, 0);
        pair_transcript_sizes.resize(1, 0);
        pair_transcripts.resize(1);
        return;
    }

    MetalBuffer<uint64_t> sizes(pair_count);

    MetalContext::get().dispatch(
        "compute_pair_intersection_sizes_kernel",
        (NSUInteger)pair_count,
        {
            intersector.read_transcripts.metalBuffer(),
            intersector.read_transcript_offsets.metalBuffer(),
            intersector.read_transcript_sizes.metalBuffer(),
            intersector.read_had_mapped_kmers.metalBuffer(),
            sizes.metalBuffer()
        },
        { as_constant((uint32_t)pair_count) }
    );

    pair_transcript_offsets.resize(pair_count + 1);
    exclusive_scan_u64_cpu(sizes.data(), pair_transcript_offsets.data(), pair_count);
    uint64_t total = prefix_scan_total(sizes.data(), pair_transcript_offsets.data(), pair_count);
    pair_transcript_offsets[pair_count] = total;

    pair_transcripts.resize(total == 0 ? 1 : total);
    pair_transcript_sizes.resize(pair_count);

    MetalContext::get().dispatch(
        "intersect_pairs_kernel",
        (NSUInteger)pair_count,
        {
            intersector.read_transcripts.metalBuffer(),
            intersector.read_transcript_offsets.metalBuffer(),
            intersector.read_transcript_sizes.metalBuffer(),
            intersector.read_had_mapped_kmers.metalBuffer(),
            pair_transcript_offsets.metalBuffer(),
            pair_transcripts.metalBuffer(),
            pair_transcript_sizes.metalBuffer()
        },
        { as_constant((uint32_t)pair_count) }
    );
}

// ============================================================================
// ReadECLookup
// ============================================================================

void ReadECLookup::lookup_impl(const MetalBuffer<int>&      transcripts,
                                const MetalBuffer<uint64_t>& offsets,
                                const MetalBuffer<uint64_t>& sizes,
                                uint64_t count,
                                MetalECMapInv& ecmapinv,
                                const MetalECMap& ecmap)
{
    read_count = count;

    if (count == 0) {
        read_ecs_final.resize(1);
        read_hashes.resize(1);
        return;
    }

    read_hashes.resize(count);
    read_ecs_final.resize(count);

    // Step 1: GPU computes FNV-1a hashes for each transcript set
    uint64_t max_tx = transcripts.size();
    MetalContext::get().dispatch(
        "hash_transcript_sets_kernel",
        (NSUInteger)count,
        {
            transcripts.metalBuffer(),  // buffer(0)
            offsets.metalBuffer(),      // buffer(1)
            sizes.metalBuffer(),        // buffer(2)
            read_hashes.metalBuffer()   // buffer(3): out_hashes
        },
        {
            as_constant(max_tx),        // buffer(4): max_tx
            as_constant(count)          // buffer(5): count
        }
    );

    // Step 2: CPU looks up each hash in ecmapinv, then verifies
    const uint64_t* hashes = read_hashes.data();
    const uint64_t* tx_off = offsets.data();
    const uint64_t* tx_sz  = sizes.data();
    const int*      tx_buf = transcripts.data();
    const int*      em_tx  = ecmap.transcripts.data();
    const uint64_t* em_off = ecmap.offsets.data();
    size_t em_num_ecs      = ecmap.num_ecs;
    size_t em_max_tx       = ecmap.transcripts.size();
    int*   ecs_out         = read_ecs_final.data();

    auto t0 = std::chrono::high_resolution_clock::now();

    for (uint64_t id = 0; id < count; ++id) {
        uint64_t h    = hashes[id];
        uint64_t start = tx_off[id];
        uint64_t sz   = tx_sz[id];

        if (sz == 0 || start + sz > max_tx || h == UINT64_MAX) {
            ecs_out[id] = -1;
            continue;
        }

        int ec_id = ecmapinv.find(h);
        if (ec_id < 0 || (size_t)ec_id >= em_num_ecs) {
            ecs_out[id] = -1;
            continue;
        }

        // Verify against forward ecmap
        uint64_t stored_start = em_off[ec_id];
        uint64_t stored_end   = em_off[ec_id + 1];
        uint64_t stored_size  = stored_end - stored_start;

        if (stored_start + stored_size > em_max_tx || stored_size != sz) {
            ecs_out[id] = -1;
            continue;
        }

        bool match = true;
        for (uint64_t i = 0; i < sz; ++i) {
            if (tx_buf[start + i] != em_tx[stored_start + i]) { match = false; break; }
        }
        ecs_out[id] = match ? ec_id : -1;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_ec_lookup_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

void ReadECLookup::lookup(const ReadTranscriptIntersector& intersector,
                          MetalECMapInv& ecmapinv,
                          const MetalECMap& ecmap)
{
    lookup_impl(intersector.read_transcripts,
                intersector.read_transcript_offsets,
                intersector.read_transcript_sizes,
                intersector.read_count,
                ecmapinv, ecmap);
}

void ReadECLookup::lookup_pairs(const PairIntersector& pairs,
                                MetalECMapInv& ecmapinv,
                                const MetalECMap& ecmap)
{
    lookup_impl(pairs.pair_transcripts,
                pairs.pair_transcript_offsets,
                pairs.pair_transcript_sizes,
                pairs.pair_count,
                ecmapinv, ecmap);
}

// ============================================================================
// ECCounter
// ============================================================================

ECCounter::ECCounter(size_t n) : num_ecs(n) {
    ec_counts.resize(n, 0);
}

void ECCounter::ensure_capacity(size_t new_num_ecs) {
    if (new_num_ecs > num_ecs) {
        ec_counts.resize(new_num_ecs, 0);
        num_ecs = new_num_ecs;
    }
}

void ECCounter::count_batch(const ReadECLookup& lookup) {
    if (lookup.read_count == 0) return;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Simple CPU scatter — on Apple Silicon the buffer is already in RAM
    const int* ecs = lookup.read_ecs_final.data();
    int* counts    = ec_counts.data();
    size_t n_ecs   = num_ecs;

    for (size_t i = 0; i < lookup.read_count; ++i) {
        int ec = ecs[i];
        if (ec >= 0 && (size_t)ec < n_ecs)
            ++counts[ec];
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_ec_counting_ms +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
}

void ECCounter::write_counts(const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open())
        throw std::runtime_error("ECCounter::write_counts: cannot open " + filename);

    const int* counts = ec_counts.data();
    for (size_t ec_id = 0; ec_id < num_ecs; ++ec_id)
        outfile << ec_id << "\t" << counts[ec_id] << "\n";
}

// ============================================================================
// NewECHandler
// ============================================================================

void NewECHandler::handle_batch(const PairIntersector& pairs,
                                ReadECLookup&  ec_lookup,
                                MetalECMapInv& ecmapinv,
                                MetalECMap&    ecmap,
                                ECCounter&     ec_counter)
{
    if (ec_lookup.read_count == 0 || pairs.pair_count == 0) return;

    uint64_t n = ec_lookup.read_count;
    const int*      ecs      = ec_lookup.read_ecs_final.data();
    const uint64_t* tx_sizes = pairs.pair_transcript_sizes.data();
    const uint64_t* hashes   = ec_lookup.read_hashes.data();
    const uint64_t* tx_off   = pairs.pair_transcript_offsets.data();
    const int*      tx_buf   = pairs.pair_transcripts.data();

    // Step 1: collect novel pair indices (ec==-1 AND tx_size>0)
    std::vector<uint64_t> novel_indices;
    for (uint64_t i = 0; i < n; ++i) {
        if (ecs[i] == -1 && tx_sizes[i] > 0)
            novel_indices.push_back(i);
    }
    if (novel_indices.empty()) return;
    total_new_read_count += (int)novel_indices.size();

    // Step 2: deduplicate hashes
    std::unordered_map<uint64_t, uint64_t> hash_to_first;
    for (uint64_t idx : novel_indices) {
        uint64_t h = hashes[idx];
        if (hash_to_first.find(h) == hash_to_first.end())
            hash_to_first[h] = idx;
    }

    // Step 3: identify truly new hashes (not in ecmapinv)
    std::vector<std::pair<uint64_t, uint64_t>> truly_new; // (hash, first_pair_idx)
    for (const auto& kv : hash_to_first) {
        if (ecmapinv.find(kv.first) < 0)
            truly_new.push_back(kv);
    }

    if (!truly_new.empty()) {
        // Step 4: insert into ecmapinv + ecmap
        std::vector<int>      all_new_tx;
        std::vector<uint64_t> new_offsets;
        new_offsets.reserve(truly_new.size() + 1);
        uint64_t tx_offset = 0;

        for (const auto& kv : truly_new) {
            uint64_t h    = kv.first;
            uint64_t pidx = kv.second;
            int new_id    = next_ec_id++;

            ecmapinv.insert(h, new_id);

            uint64_t start = tx_off[pidx];
            uint64_t sz    = tx_sizes[pidx];
            new_offsets.push_back(tx_offset);
            all_new_tx.insert(all_new_tx.end(), tx_buf + start, tx_buf + start + sz);
            tx_offset += sz;
        }
        new_offsets.push_back(tx_offset);

        ecmap.append_ecs(all_new_tx, new_offsets, (int)truly_new.size());
        ec_counter.ensure_capacity(next_ec_id);
        total_new_ec_count += (int)truly_new.size();
    }

    // Step 5: re-look-up all novel pairs and write back into read_ecs_final
    int* ecs_out = ec_lookup.read_ecs_final.data();
    const int*      em_tx  = ecmap.transcripts.data();
    const uint64_t* em_off = ecmap.offsets.data();
    size_t em_num_ecs      = ecmap.num_ecs;
    size_t em_max_tx       = ecmap.transcripts.size();

    for (uint64_t idx : novel_indices) {
        uint64_t h  = hashes[idx];
        int ec_id   = ecmapinv.find(h);
        if (ec_id < 0 || (size_t)ec_id >= em_num_ecs) { ecs_out[idx] = -1; continue; }

        uint64_t start = tx_off[idx];
        uint64_t sz    = tx_sizes[idx];
        uint64_t stored_start = em_off[ec_id];
        uint64_t stored_end   = em_off[ec_id + 1];
        uint64_t stored_size  = stored_end - stored_start;

        if (stored_size != sz || stored_start + stored_size > em_max_tx) {
            ecs_out[idx] = -1; continue;
        }

        bool match = true;
        for (uint64_t i = 0; i < sz; ++i) {
            if (tx_buf[start + i] != em_tx[stored_start + i]) { match = false; break; }
        }
        ecs_out[idx] = match ? ec_id : -1;
    }
}
