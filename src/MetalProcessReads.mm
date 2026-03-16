#import "MetalProcessReads.h"
#import "MetalUtils.h"
#import "MetalIndex.h"
#import "MetalPipeline.h"
#import "MetalEM.h"
#include "BenchmarkStats.h"
#include "PlaintextWriter.h"
#include "EMAlgorithm.h"
#include "weights.h"
#include "MinCollector.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <algorithm>
#include <cmath>

// Global benchmark stats (declared in BenchmarkStats.h)
BenchmarkStats g_benchmark_stats;

// ============================================================================
// Internal: write EC counts + transcript lists (mirrors CUDA version's output)
// ============================================================================

static void write_ec_counts_with_transcripts(const std::string& filename,
                                              const ECCounter& ec_counter,
                                              const MetalECMap& ecmap)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Warning: could not open " << filename << std::endl;
        return;
    }

    const int*      counts = ec_counter.ec_counts.data();
    const int*      tx     = ecmap.transcripts.data();
    const uint64_t* off    = ecmap.offsets.data();
    size_t n_ecs = ec_counter.num_ecs;

    for (size_t ec_id = 0; ec_id < n_ecs; ++ec_id) {
        outfile << ec_id << "\t" << counts[ec_id] << "\t";
        uint64_t start = off[ec_id];
        uint64_t end   = off[ec_id + 1];
        for (uint64_t i = start; i < end; ++i) {
            if (i > start) outfile << ',';
            outfile << tx[i];
        }
        outfile << '\n';
    }
}

// ============================================================================
// Internal: per-batch processing
// ============================================================================

static void process_batch(DeviceKmerLoader&       loader,
                          const MetalKmerTable&    kmer_table,
                          MetalBuffer<int>&         d_ecs,
                          ReadECCollapser&          collapser,
                          ReadTranscriptIntersector& intersector,
                          PairIntersector&           pair_intersector,
                          ReadECLookup&              ec_lookup,
                          ECCounter&                 ec_counter,
                          NewECHandler&              new_ec_handler,
                          MetalECMap&                ecmap,
                          MetalECMapInv&             ecmapinv)
{
    if (loader.read_count == 0) return;

    // Step 1: Extract k-mers (kmer_kernel)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        loader.run(kmer_table.empty_sentinel, Kmer::k);
        auto t1 = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.gpu_kmer_extraction_ms +=
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    }

    // Step 2: Hash table lookup (kmer_lookup_kernel)
    {
        uint64_t total_kmers = loader.kmers.size();
        d_ecs.resize(total_kmers == 0 ? 1 : total_kmers);
        g_benchmark_stats.total_kmers += total_kmers;

        auto t0 = std::chrono::high_resolution_clock::now();
        MetalContext::get().dispatch(
            "kmer_lookup_kernel",
            (NSUInteger)total_kmers,
            {
                kmer_table.slots.metalBuffer(),
                loader.kmers.metalBuffer(),
                d_ecs.metalBuffer()
            },
            {
                as_constant(kmer_table.capacity),
                as_constant(kmer_table.empty_sentinel),
                as_constant(total_kmers)
            }
        );
        auto t1 = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.gpu_kmer_lookup_ms +=
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    }

    // Step 3: EC dedup per read (collapse_ecs_per_read_kernel)
    collapser.collapse(loader, d_ecs);

    // 4. Intersect with EC->transcript map
    intersector.intersect(collapser, ecmap);

    bool is_paired = (loader.r1_count > 0 &&
                      loader.read_count > loader.r1_count);

    if (is_paired) {
        // 5a. Paired-end intersection
        pair_intersector.intersect_pairs(intersector, loader.r1_count);

        // 5b. EC lookup for pairs
        ec_lookup.lookup_pairs(pair_intersector, ecmapinv, ecmap);

        // 5c. Handle novel ECs
        new_ec_handler.handle_batch(pair_intersector, ec_lookup, ecmapinv, ecmap, ec_counter);
    } else {
        // 5. Single-end EC lookup
        ec_lookup.lookup(intersector, ecmapinv, ecmap);
    }

    // 6. Count batch
    ec_counter.count_batch(ec_lookup);
}

// ============================================================================
// Fill run stats from accumulated counters
// ============================================================================

static void fill_run_stats(int64_t total_processed,
                           const ECCounter& ec_counter,
                           const MetalECMap& ecmap,
                           MetalRunStats* out_stats)
{
    if (!out_stats) return;

    out_stats->num_processed = total_processed;

    const int* counts = ec_counter.ec_counts.data();
    size_t n_ecs = ec_counter.num_ecs;

    int64_t n_aligned = 0;
    int64_t n_unique  = 0;
    for (size_t e = 0; e < n_ecs; ++e) {
        int c = counts[e];
        if (c > 0) {
            n_aligned += c;
            uint64_t start = ecmap.offsets[e];
            uint64_t end   = ecmap.offsets[e + 1];
            if (end - start == 1) n_unique += c;
        }
    }
    out_stats->num_pseudoaligned = n_aligned;
    out_stats->num_unique        = n_unique;
}

// ============================================================================
// Benchmark summary printout
// ============================================================================

static void print_benchmark_summary() {
    auto& s = g_benchmark_stats;
    std::cerr << "\n[benchmark] Metal pipeline summary\n";
    std::cerr << "  setup_index_load:         " << s.setup_index_load_ms << " ms\n";
    std::cerr << "  setup_kmer_table:         " << s.setup_kmer_to_ec_map_ms << " ms\n";
    std::cerr << "  setup_ecmap:              " << s.setup_gpu_ecmap_ms << " ms\n";
    std::cerr << "  setup_ecmapinv:           " << s.setup_gpu_ecmapinv_ms << " ms\n";
    std::cerr << "  kmer_extraction:          " << s.gpu_kmer_extraction_ms << " ms\n";
    std::cerr << "  kmer_lookup:              " << s.gpu_kmer_lookup_ms << " ms\n";
    std::cerr << "  ec_collapse:              " << s.gpu_ec_collapse_ms << " ms\n";
    std::cerr << "  transcript_intersection:  " << s.gpu_transcript_intersection_ms << " ms\n";
    std::cerr << "  ec_lookup:                " << s.gpu_ec_lookup_ms << " ms\n";
    std::cerr << "  ec_counting:              " << s.gpu_ec_counting_ms << " ms\n";
    std::cerr << "  em:                       " << s.gpu_em_ms << " ms\n";
    std::cerr << "  h2d_copy:                 " << s.host_h2d_copy_ms << " ms\n";
    std::cerr << "  io_decompress:            " << s.io_decompress_ms << " ms\n";
    std::cerr << "  wall_clock_pipeline:      " << s.wall_clock_pipeline_ms << " ms\n";
    std::cerr << "  wall_clock_total:         " << s.wall_clock_total_ms << " ms\n";
    std::cerr << "  batches:                  " << s.batch_count << "\n";
    if (s.batch_count > 0)
        std::cerr << "  reads_per_batch:          "
                  << s.total_kmers / s.batch_count << "\n";
}

// ============================================================================
// Effective length computation (mirrors CUDA version)
// ============================================================================

static std::vector<double> compute_eff_lens(const ProgramOptions& opt,
                                             const std::vector<uint32_t>& target_lens)
{
    size_t num_trans = target_lens.size();
    std::vector<double> eff_lens(num_trans);

    if (opt.fld > 0.0 && opt.sd > 0.0) {
        auto mean_fl_trunc = trunc_gaussian_fld(0, MAX_FRAG_LEN, opt.fld, opt.sd);
        auto fl_means      = get_frag_len_means(target_lens, mean_fl_trunc);
        eff_lens           = calc_eff_lens(target_lens, fl_means);
    } else {
        double mean_fld = opt.fld;
        for (size_t t = 0; t < num_trans; ++t) {
            double len = static_cast<double>(target_lens[t]);
            double eff = len - mean_fld + 1.0;
            eff_lens[t] = (eff < 1.0) ? len : eff;
        }
    }
    return eff_lens;
}

// ============================================================================
// Shared pipeline core
// ============================================================================

struct MetalPipelineData {
    MetalKmerTable    kmer_table;
    MetalECMap        ecmap;
    MetalECMapInv     ecmapinv;
    ECCounter         ec_counter;
    NewECHandler      new_ec_handler;
    MetalBuffer<int>  d_ecs;
    ReadECCollapser       collapser;
    ReadTranscriptIntersector intersector;
    PairIntersector   pair_intersector;
    ReadECLookup      ec_lookup;
    int64_t           total_processed = 0;

    MetalPipelineData(size_t n_ecs)
        : ec_counter(n_ecs), new_ec_handler((int)n_ecs) {}
};

// Double-buffered batch processing loop (identical logic to CUDA version)
static void run_batch_loop(ProgramOptions& opt,
                           MetalPipelineData& pd,
                           int64_t& total_processed,
                           std::chrono::high_resolution_clock::time_point& pipeline_start,
                           std::chrono::high_resolution_clock::time_point& pipeline_end)
{
    // Larger batches → more concurrent DRAM requests → better GPU latency hiding.
    // Sweet spot at 300K pairs (600K reads):
    //   kmers: 600K×120×8B = 576MB, d_ecs: 600K×120×4B = 288MB
    //   temp_ecs_: 600K×512×4B = 1.23GB, reads: ~90MB → ~2.2GB total (+ 2GB table = 4.2GB)
    // Beyond 300K, memory pressure increases without throughput gain (kmer_lookup saturates DRAM BW).
    const uint32_t BATCH_SIZE = 300000;

    // Single shared reader — the loader thread is the only one that calls it,
    // so no locking on the reader itself is needed.
    FastqSequenceReader shared_reader(opt);

    HostReadLoader h_loaders[2] = {
        HostReadLoader(opt, BATCH_SIZE),
        HostReadLoader(opt, BATCH_SIZE)
    };
    DeviceKmerLoader d_loaders[2];
    int curr = 0;

    std::mutex mtx;
    std::condition_variable cv;
    enum LoadState { IDLE, REQUESTED, DONE, SHUTDOWN };
    LoadState load_state = IDLE;
    int  load_buf      = 0;
    bool load_has_data = false;

    std::thread loader_thread([&]() {
        while (true) {
            {
                std::unique_lock<std::mutex> lock(mtx);
                cv.wait(lock, [&]{ return load_state == REQUESTED || load_state == SHUTDOWN; });
                if (load_state == SHUTDOWN) break;
            }
            bool ok = h_loaders[load_buf].load(shared_reader);
            if (ok) d_loaders[load_buf].load(h_loaders[load_buf]);
            {
                std::lock_guard<std::mutex> lock(mtx);
                load_has_data = ok;
                load_state    = DONE;
            }
            cv.notify_one();
        }
    });

    auto request_load = [&](int buf) {
        std::lock_guard<std::mutex> lock(mtx);
        load_buf   = buf;
        load_state = REQUESTED;
        cv.notify_one();
    };

    auto wait_load = [&]() -> bool {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return load_state == DONE; });
        load_state = IDLE;
        return load_has_data;
    };

    // Kick off first batch load (overlaps with index build above)
    request_load(curr);
    bool has_data = wait_load();

    pipeline_start = std::chrono::high_resolution_clock::now();

    while (has_data) {
        int next = 1 - curr;
        request_load(next);

        process_batch(d_loaders[curr], pd.kmer_table, pd.d_ecs,
                      pd.collapser, pd.intersector, pd.pair_intersector,
                      pd.ec_lookup, pd.ec_counter, pd.new_ec_handler,
                      pd.ecmap, pd.ecmapinv);

        auto& dl = d_loaders[curr];
        if (dl.r1_count > 0 && dl.read_count > dl.r1_count)
            total_processed += (int64_t)dl.r1_count;
        else
            total_processed += (int64_t)dl.read_count;

        g_benchmark_stats.batch_count++;

        has_data = wait_load();
        curr = next;
    }

    pipeline_end = std::chrono::high_resolution_clock::now();

    {
        std::lock_guard<std::mutex> lock(mtx);
        load_state = SHUTDOWN;
    }
    cv.notify_one();
    loader_thread.join();
}

// ============================================================================
// metal_run() — KmerIndex variant
// ============================================================================

void metal_run(ProgramOptions& opt, const KmerIndex& index, MetalRunStats* out_stats) {
    auto wall_start = std::chrono::high_resolution_clock::now();

    size_t n_ecs = index.ecmapinv.size();
    MetalPipelineData pd(n_ecs);

    std::string cache_path = opt.index + ".metalcache";
    std::cerr << "[metal] Building k-mer hash table..." << std::endl;
    if (!pd.kmer_table.load_cache(cache_path)) {
        pd.kmer_table.build(index);
        std::cerr << "[metal] Saving k-mer table cache to " << cache_path << " ..." << std::flush;
        if (pd.kmer_table.save_cache(cache_path))
            std::cerr << " done" << std::endl;
        else
            std::cerr << " FAILED (will rebuild next run)" << std::endl;
    } else {
        std::cerr << "[metal] k-mer table loaded from cache: " << cache_path << std::endl;
    }

    std::cerr << "[metal] Building EC map..." << std::endl;
    pd.ecmap.build(index);

    std::cerr << "[metal] Building EC inverse map..." << std::endl;
    pd.ecmapinv.build(index);

    std::chrono::high_resolution_clock::time_point pipeline_start, pipeline_end;
    run_batch_loop(opt, pd, pd.total_processed, pipeline_start, pipeline_end);

    // EM
    {
        size_t num_trans = index.num_trans;
        auto eff_lens = compute_eff_lens(opt, index.target_lens_);

        auto em_start = std::chrono::high_resolution_clock::now();
        MetalEM em(num_trans, pd.ecmap.num_ecs, eff_lens);
        em.build_transpose(pd.ecmap);
        em.run_transpose(pd.ecmap, pd.ec_counter.ec_counts, opt.iterations, 50);
        auto em_end = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.gpu_em_ms =
            std::chrono::duration_cast<std::chrono::microseconds>(em_end - em_start).count() / 1000.0;

        // Output
        if (!opt.output.empty()) {
            struct stat st;
            if (stat(opt.output.c_str(), &st) != 0)
                mkdir(opt.output.c_str(), 0777);

            write_ec_counts_with_transcripts(opt.output + "/counts_metal.txt",
                                             pd.ec_counter, pd.ecmap);

            if (pd.new_ec_handler.total_new_ecs() > 0)
                std::cerr << "  Novel ECs from paired-end: "
                          << pd.new_ec_handler.total_new_ecs()
                          << " (" << pd.new_ec_handler.total_new_reads() << " read pairs)\n";

            auto alpha_f = em.d_alpha.to_host();
            std::vector<double> alpha(alpha_f.begin(), alpha_f.end());
            plaintext_writer(opt.output + "/abundance.tsv",
                             index.target_names_, alpha, eff_lens, index.target_lens_);
        }

        fill_run_stats(pd.total_processed, pd.ec_counter, pd.ecmap, out_stats);
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.wall_clock_total_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start).count() / 1000.0;
    g_benchmark_stats.wall_clock_pipeline_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(pipeline_end - pipeline_start).count() / 1000.0;

    print_benchmark_summary();
}

// ============================================================================
// metal_run() — GPUIndex variant
// ============================================================================

void metal_run(ProgramOptions& opt, const GPUIndex& index, MetalRunStats* out_stats) {
    auto wall_start = std::chrono::high_resolution_clock::now();

    size_t n_ecs = index.ecmapinv.size();
    MetalPipelineData pd(n_ecs);

    std::cerr << "[metal] Building k-mer sorted table from GPU index..." << std::endl;
    pd.kmer_table.build(index);

    std::cerr << "[metal] Building EC map..." << std::endl;
    pd.ecmap.build(index);

    std::cerr << "[metal] Building EC inverse map..." << std::endl;
    pd.ecmapinv.build(index);

    std::chrono::high_resolution_clock::time_point pipeline_start, pipeline_end;
    run_batch_loop(opt, pd, pd.total_processed, pipeline_start, pipeline_end);

    // EM
    {
        size_t num_trans = index.num_transcripts;
        auto eff_lens = compute_eff_lens(opt, index.target_lens_);

        auto em_start = std::chrono::high_resolution_clock::now();
        MetalEM em(num_trans, pd.ecmap.num_ecs, eff_lens);
        em.build_transpose(pd.ecmap);
        em.run_transpose(pd.ecmap, pd.ec_counter.ec_counts, opt.iterations, 50);
        auto em_end = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.gpu_em_ms =
            std::chrono::duration_cast<std::chrono::microseconds>(em_end - em_start).count() / 1000.0;

        // Output
        if (!opt.output.empty()) {
            struct stat st;
            if (stat(opt.output.c_str(), &st) != 0)
                mkdir(opt.output.c_str(), 0777);

            write_ec_counts_with_transcripts(opt.output + "/counts_metal.txt",
                                             pd.ec_counter, pd.ecmap);

            if (pd.new_ec_handler.total_new_ecs() > 0)
                std::cerr << "  Novel ECs from paired-end: "
                          << pd.new_ec_handler.total_new_ecs()
                          << " (" << pd.new_ec_handler.total_new_reads() << " read pairs)\n";

            auto alpha_f = em.d_alpha.to_host();
            std::vector<double> alpha(alpha_f.begin(), alpha_f.end());
            plaintext_writer(opt.output + "/abundance.tsv",
                             index.target_names_, alpha, eff_lens, index.target_lens_);
        }

        fill_run_stats(pd.total_processed, pd.ec_counter, pd.ecmap, out_stats);
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.wall_clock_total_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start).count() / 1000.0;
    g_benchmark_stats.wall_clock_pipeline_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(pipeline_end - pipeline_start).count() / 1000.0;

    print_benchmark_summary();
}
