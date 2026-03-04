#include "GPUProcessReads.cuh"
#include "GPUIndex.cuh"
#include "GPUKernels.cuh"
#include "GPUPipeline.cuh"
#include "GPUReadLoader.cuh"
#include "GPUEM.cuh"
#include "PlaintextWriter.h"
#include "MinCollector.h"
#include "EMAlgorithm.h"
#include "weights.h"
// #include "Kmer.hpp"
#include <thrust/device_vector.h>
#include <thrust/count.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <map>
#include <set>
#include <cuda_runtime.h>

// Global benchmark statistics instance (defined here, declared in BenchmarkStats.h)
BenchmarkStats g_benchmark_stats;

// Write EC counts with transcript lists (matching CPU counts_cpu.txt format)
static void write_ec_counts_with_transcripts(
    const std::string& filename,
    const ECCounter& ec_counter,
    const GPUECMap& gpu_ecmap)
{
  size_t num_ecs = ec_counter.ec_counts.size();
  std::vector<int> h_counts(num_ecs);
  thrust::copy(ec_counter.ec_counts.begin(), ec_counter.ec_counts.end(), h_counts.begin());

  std::vector<int> h_transcripts(gpu_ecmap.transcripts.size());
  thrust::copy(gpu_ecmap.transcripts.begin(), gpu_ecmap.transcripts.end(), h_transcripts.begin());

  std::vector<uint64_t> h_offsets(gpu_ecmap.offsets.size());
  thrust::copy(gpu_ecmap.offsets.begin(), gpu_ecmap.offsets.end(), h_offsets.begin());

  std::ofstream outfile(filename);
  for (size_t ec_id = 0; ec_id < num_ecs; ++ec_id) {
    outfile << ec_id << "\t" << h_counts[ec_id] << "\t";
    uint64_t start = h_offsets[ec_id];
    uint64_t end   = h_offsets[ec_id + 1];
    for (uint64_t i = start; i < end; ++i) {
      if (i > start) outfile << ",";
      outfile << h_transcripts[i];
    }
    outfile << "\n";
  }
  outfile.close();
  std::cout << "  Wrote EC counts to: " << filename << std::endl;
}

// DeviceKmerLoader::load_from_gpu implementation
void DeviceKmerLoader::load_from_gpu(GPUReadLoader& gpu) {
    auto start_time = std::chrono::high_resolution_clock::now();

    read_count = gpu.read_count;
    r1_count = gpu.r1_count;
    if (read_count == 0) return;

    read_id_to_offset.resize(gpu.read_count);
    thrust::copy(gpu.read_id_to_offset.begin(), gpu.read_id_to_offset.end(), read_id_to_offset.begin());

    read_id_to_length.resize(gpu.read_count);
    thrust::copy(gpu.read_id_to_length.begin(), gpu.read_id_to_length.end(), read_id_to_length.begin());

    read_id_to_kmer_first.resize(gpu.read_count);
    thrust::copy(gpu.read_id_to_kmer_first.begin(), gpu.read_id_to_kmer_first.end(), read_id_to_kmer_first.begin());

    read_id_to_kmer_count.resize(gpu.read_count);
    thrust::copy(gpu.read_id_to_kmer_count.begin(), gpu.read_id_to_kmer_count.end(), read_id_to_kmer_count.begin());

    // Copy the reads data (FASTQ text, already on GPU - device-to-device copy)
    reads.resize(gpu.d_reads_size);
    cudaMemcpyAsync(reads.data().get(), gpu.d_reads, gpu.d_reads_size, cudaMemcpyDeviceToDevice, 0);

    kmers.resize(gpu.kmer_count);

    // Ensure all D2D copies on this thread's stream are complete before
    // the main thread starts reading from these buffers.
    // With --default-stream per-thread, stream 0 = this thread's default stream.
    cudaStreamSynchronize(0);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    g_benchmark_stats.host_h2d_copy_ms += duration.count() / 1000.0;
}

// Forward declaration
void print_benchmark_summary();

// Fill run stats from ec_counter and gpu_ecmap (num_processed set by caller)
static void fill_run_stats(int64_t num_processed,
    const ECCounter& ec_counter, const GPUECMap& gpu_ecmap, GPURunStats* out) {
  if (!out) return;
  out->num_processed = num_processed;
  size_t num_ecs = ec_counter.ec_counts.size();
  std::vector<int> h_counts(num_ecs);
  thrust::copy(ec_counter.ec_counts.begin(), ec_counter.ec_counts.end(), h_counts.begin());
  std::vector<uint64_t> h_offsets(gpu_ecmap.offsets.size());
  thrust::copy(gpu_ecmap.offsets.begin(), gpu_ecmap.offsets.end(), h_offsets.begin());
  int64_t pseudo = 0, uniq = 0;
  for (size_t ec_id = 0; ec_id < num_ecs; ++ec_id) {
    int c = h_counts[ec_id];
    pseudo += c;
    uint64_t num_tx = (ec_id + 1 < h_offsets.size()) ? (h_offsets[ec_id + 1] - h_offsets[ec_id]) : 0;
    if (num_tx == 1) uniq += c;
  }
  out->num_pseudoaligned = pseudo;
  out->num_unique = uniq;
}

// Process a single batch through the GPU pipeline
static void process_batch(
    DeviceKmerLoader& d_loader,
    thrust::device_vector<int>& d_ecs,
    cuco::static_map<uint64_t, int>& d_map,
    ReadECCollapser& collapser,
    ReadTranscriptIntersector& intersector,
    PairIntersector& pair_intersector,
    ReadECLookup& ec_lookup,
    ECCounter& ec_counter,
    NewECHandler& new_ec_handler,
    GPUECMap& gpu_ecmap,
    GPUECMapInv& gpu_ecmapinv
) {
    bool is_paired = (d_loader.r1_count > 0);

    // Kmer extraction
    cudaEvent_t kmer_start, kmer_stop;
    cudaEventCreate(&kmer_start);
    cudaEventCreate(&kmer_stop);
    cudaEventRecord(kmer_start);
    d_loader.run();
    cudaEventRecord(kmer_stop);
    cudaEventSynchronize(kmer_stop);
    float kmer_ms = 0;
    cudaEventElapsedTime(&kmer_ms, kmer_start, kmer_stop);
    g_benchmark_stats.gpu_kmer_extraction_ms += kmer_ms;
    cudaEventDestroy(kmer_start);
    cudaEventDestroy(kmer_stop);

    // Kmer lookup
    d_ecs.resize(d_loader.kmers.size(), -1);
    cudaEvent_t lookup_start, lookup_stop;
    cudaEventCreate(&lookup_start);
    cudaEventCreate(&lookup_stop);
    cudaEventRecord(lookup_start);
    d_map.find(d_loader.kmers.begin(),
               d_loader.kmers.end(),
               d_ecs.begin());
    cudaEventRecord(lookup_stop);
    cudaEventSynchronize(lookup_stop);
    float lookup_ms = 0;
    cudaEventElapsedTime(&lookup_ms, lookup_start, lookup_stop);
    g_benchmark_stats.gpu_kmer_lookup_ms += lookup_ms;
    cudaEventDestroy(lookup_start);
    cudaEventDestroy(lookup_stop);

    // EC collapse (internally syncs via cudaEventSynchronize)
    collapser.collapse(d_loader, d_ecs);

    // Transcript intersection per read (internally syncs via cudaEventSynchronize)
    intersector.intersect_gpu(collapser, gpu_ecmap);

    if (is_paired) {
        // Pair intersection: intersect R1[i] and R2[i] transcript sets
        pair_intersector.intersect_pairs(intersector, d_loader.r1_count);

        // EC lookup on pair results (with verification against GPUECMap)
        ec_lookup.lookup_pairs(pair_intersector, gpu_ecmapinv, gpu_ecmap);

        // Handle novel ECs: filter, dedup, insert into dynamic_map, re-find
        new_ec_handler.handle_batch(pair_intersector, ec_lookup,
                                     gpu_ecmapinv, gpu_ecmap, ec_counter);

        // EC counting (one count per pair, now includes novel ECs)
        ec_counter.count_batch(ec_lookup);
    } else {
        // Single-end: pass per-read transcript sets directly to EC lookup
        ec_lookup.lookup(intersector, gpu_ecmapinv, gpu_ecmap);
        ec_counter.count_batch(ec_lookup);
    }
}

void gpu_run(ProgramOptions& opt, const KmerIndex& index, GPURunStats* out_stats) {
  auto wall_start = std::chrono::high_resolution_clock::now();

  GPUReadLoader gpu_loader(opt);

  // Double-buffered pipeline with persistent loader thread.
  // With --default-stream per-thread, each thread gets its own CUDA stream,
  // so the loader's D2D copies overlap with the main thread's compute.
  DeviceKmerLoader d_loaders[2];
  int curr = 0;
  int64_t total_processed = 0;

  // Loader thread communication
  std::mutex mtx;
  std::condition_variable cv;
  enum LoadState { IDLE, REQUESTED, DONE, SHUTDOWN };
  LoadState load_state = IDLE;
  int load_buf = 0;
  bool load_has_data = false;

  std::thread loader_thread([&]() {
    while (true) {
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return load_state == REQUESTED || load_state == SHUTDOWN; });
        if (load_state == SHUTDOWN) break;
      }
      bool ok = gpu_loader.load();
      if (ok) {
        d_loaders[load_buf].load_from_gpu(gpu_loader);
      }
      {
        std::lock_guard<std::mutex> lock(mtx);
        load_has_data = ok;
        load_state = DONE;
      }
      cv.notify_one();
    }
  });

  auto request_load = [&](int buf) {
    std::lock_guard<std::mutex> lock(mtx);
    load_buf = buf;
    load_state = REQUESTED;
    cv.notify_one();
  };

  auto wait_load = [&]() -> bool {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [&]{ return load_state == DONE; });
    load_state = IDLE;
    return load_has_data;
  };

  // Start loading first batch in background - file reading overlaps with GPU setup below
  request_load(curr);

  // Build GPU data structures while files are being read from disk
  auto d_map = build_kmer_to_ec_map(index);

  std::cout << "[Building GPU EC Map]" << std::endl;
  GPUECMap gpu_ecmap(index);

  std::cout << "[Building GPU EC Map Inverse]" << std::endl;
  GPUECMapInv gpu_ecmapinv(index);

  ECCounter ec_counter(index.ecmapinv.size());
  NewECHandler new_ec_handler(static_cast<int>(index.ecmapinv.size()));

  thrust::device_vector<int> d_ecs;
  ReadECCollapser collapser;
  ReadTranscriptIntersector intersector;
  PairIntersector pair_intersector;
  ReadECLookup ec_lookup;

  // Wait for first batch (file reading should have overlapped with GPU setup)
  bool has_data = wait_load();

  auto pipeline_start = std::chrono::high_resolution_clock::now();

  while (has_data) {
    int next = 1 - curr;

    // Request next batch loading (overlaps with compute below)
    request_load(next);

    // Process current batch on this thread's default stream
    process_batch(d_loaders[curr], d_ecs, d_map, collapser, intersector,
                  pair_intersector, ec_lookup, ec_counter, new_ec_handler,
                  gpu_ecmap, gpu_ecmapinv);
    // Count pairs for paired-end, reads for single-end (n_processed = fragments)
    if (d_loaders[curr].r1_count > 0 && d_loaders[curr].read_count > d_loaders[curr].r1_count) {
      total_processed += static_cast<int64_t>(d_loaders[curr].r1_count);
    } else {
      total_processed += static_cast<int64_t>(d_loaders[curr].read_count);
    }
    g_benchmark_stats.batch_count++;

    // Wait for next batch loading to complete
    has_data = wait_load();
    curr = next;
  }

  auto pipeline_end = std::chrono::high_resolution_clock::now();

  // Shutdown loader thread
  {
    std::lock_guard<std::mutex> lock(mtx);
    load_state = SHUTDOWN;
  }
  cv.notify_one();
  loader_thread.join();

  // EM algorithm
  {
    size_t num_trans = index.num_trans;

    // Compute effective lengths using same truncated gaussian FLD as CPU
    // (MinCollector::init_mean_fl_trunc + get_frag_len_means + calc_eff_lens)
    std::vector<double> eff_lens(num_trans);
    if (opt.fld > 0.0 && opt.sd > 0.0) {
      auto mean_fl_trunc = trunc_gaussian_fld(0, MAX_FRAG_LEN, opt.fld, opt.sd);
      auto fl_means = get_frag_len_means(index.target_lens_, mean_fl_trunc);
      eff_lens = calc_eff_lens(index.target_lens_, fl_means);
    } else {
      double mean_fld = opt.fld;
      for (size_t t = 0; t < num_trans; ++t) {
        double len = static_cast<double>(index.target_lens_[t]);
        double eff = len - mean_fld + 1.0;
        eff_lens[t] = (eff < 1.0) ? len : eff;
      }
    }

    auto em_start = std::chrono::high_resolution_clock::now();

    GPUEM em(num_trans, gpu_ecmap.num_ecs, eff_lens);
    em.build_transpose(gpu_ecmap);
    int em_rounds = em.run_transpose(gpu_ecmap, ec_counter.ec_counts, opt.iterations, 50);

    auto em_end = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_em_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(em_end - em_start).count() / 1000.0;

    // Copy alpha (estimated counts) to host
    std::vector<double> alpha(num_trans);
    thrust::copy(em.d_alpha.begin(), em.d_alpha.end(), alpha.begin());

    // Write output
    if (!opt.output.empty()) {
      struct stat st;
      if (stat(opt.output.c_str(), &st) != 0) {
        #ifdef _WIN32
          _mkdir(opt.output.c_str());
        #else
          mkdir(opt.output.c_str(), 0777);
        #endif
      }

      write_ec_counts_with_transcripts(opt.output + "/counts_gpu.txt",
                                       ec_counter, gpu_ecmap);

      if (new_ec_handler.total_new_ecs() > 0) {
        std::cout << "  New ECs from paired-end: " << new_ec_handler.total_new_ecs()
                  << " (covering " << new_ec_handler.total_new_reads() << " read pairs)" << std::endl;
      }

      plaintext_writer(opt.output + "/abundance.tsv",
                       index.target_names_, alpha, eff_lens, index.target_lens_);
    }
  }

  fill_run_stats(total_processed, ec_counter, gpu_ecmap, out_stats);

  auto wall_end = std::chrono::high_resolution_clock::now();
  g_benchmark_stats.wall_clock_total_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start).count() / 1000.0;
  g_benchmark_stats.wall_clock_pipeline_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(pipeline_end - pipeline_start).count() / 1000.0;

  print_benchmark_summary();
}

void gpu_run(ProgramOptions& opt, const GPUIndex& index, GPURunStats* out_stats) {
  auto wall_start = std::chrono::high_resolution_clock::now();

  GPUReadLoader gpu_loader(opt);

  // Double-buffered pipeline with persistent loader thread
  DeviceKmerLoader d_loaders[2];
  int curr = 0;

  // Loader thread communication
  std::mutex mtx;
  std::condition_variable cv;
  enum LoadState { IDLE, REQUESTED, DONE, SHUTDOWN };
  LoadState load_state = IDLE;
  int load_buf = 0;
  bool load_has_data = false;

  std::thread loader_thread([&]() {
    while (true) {
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return load_state == REQUESTED || load_state == SHUTDOWN; });
        if (load_state == SHUTDOWN) break;
      }
      bool ok = gpu_loader.load();
      if (ok) {
        d_loaders[load_buf].load_from_gpu(gpu_loader);
      }
      {
        std::lock_guard<std::mutex> lock(mtx);
        load_has_data = ok;
        load_state = DONE;
      }
      cv.notify_one();
    }
  });

  auto request_load = [&](int buf) {
    std::lock_guard<std::mutex> lock(mtx);
    load_buf = buf;
    load_state = REQUESTED;
    cv.notify_one();
  };

  auto wait_load = [&]() -> bool {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [&]{ return load_state == DONE; });
    load_state = IDLE;
    return load_has_data;
  };

  // Start loading first batch in background - file reading overlaps with GPU setup below
  request_load(curr);

  int64_t total_processed = 0;

  // Build GPU data structures while files are being read from disk
  std::cout << "[Building GPU k-mer to EC map from contigs]" << std::endl;
  auto d_map = build_kmer_to_ec_map_from_contigs(index);

  std::cout << "[Building GPU EC Map]" << std::endl;
  GPUECMap gpu_ecmap(index);

  std::cout << "[Building GPU EC Map Inverse]" << std::endl;
  GPUECMapInv gpu_ecmapinv(index);

  ECCounter ec_counter(index.ecmapinv.size());
  NewECHandler new_ec_handler(static_cast<int>(index.ecmapinv.size()));

  thrust::device_vector<int> d_ecs;
  ReadECCollapser collapser;
  ReadTranscriptIntersector intersector;
  PairIntersector pair_intersector;
  ReadECLookup ec_lookup;

  // Wait for first batch (file reading should have overlapped with GPU setup)
  bool has_data = wait_load();

  auto pipeline_start = std::chrono::high_resolution_clock::now();

  while (has_data) {
    int next = 1 - curr;

    // Request next batch loading (overlaps with compute below)
    request_load(next);

    // Process current batch on this thread's default stream
    process_batch(d_loaders[curr], d_ecs, d_map, collapser, intersector,
                  pair_intersector, ec_lookup, ec_counter, new_ec_handler,
                  gpu_ecmap, gpu_ecmapinv);
    // Count pairs for paired-end, reads for single-end (n_processed = fragments)
    if (d_loaders[curr].r1_count > 0 && d_loaders[curr].read_count > d_loaders[curr].r1_count) {
      total_processed += static_cast<int64_t>(d_loaders[curr].r1_count);
    } else {
      total_processed += static_cast<int64_t>(d_loaders[curr].read_count);
    }
    g_benchmark_stats.batch_count++;

    // Wait for next batch loading to complete
    has_data = wait_load();
    curr = next;
  }

  auto pipeline_end = std::chrono::high_resolution_clock::now();

  // Shutdown loader thread
  {
    std::lock_guard<std::mutex> lock(mtx);
    load_state = SHUTDOWN;
  }
  cv.notify_one();
  loader_thread.join();

  // EM algorithm
  {
    size_t num_trans = index.num_transcripts;

    // Compute effective lengths using same truncated gaussian FLD as CPU
    // (MinCollector::init_mean_fl_trunc + get_frag_len_means + calc_eff_lens)
    std::vector<double> eff_lens(num_trans);
    if (opt.fld > 0.0 && opt.sd > 0.0) {
      auto mean_fl_trunc = trunc_gaussian_fld(0, MAX_FRAG_LEN, opt.fld, opt.sd);
      auto fl_means = get_frag_len_means(index.target_lens_, mean_fl_trunc);
      eff_lens = calc_eff_lens(index.target_lens_, fl_means);
    } else {
      double mean_fld = opt.fld;
      for (size_t t = 0; t < num_trans; ++t) {
        double len = static_cast<double>(index.target_lens_[t]);
        double eff = len - mean_fld + 1.0;
        eff_lens[t] = (eff < 1.0) ? len : eff;
      }
    }

    auto em_start = std::chrono::high_resolution_clock::now();

    GPUEM em(num_trans, gpu_ecmap.num_ecs, eff_lens);
    em.build_transpose(gpu_ecmap);
    int em_rounds = em.run_transpose(gpu_ecmap, ec_counter.ec_counts, opt.iterations, 50);

    auto em_end = std::chrono::high_resolution_clock::now();
    g_benchmark_stats.gpu_em_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(em_end - em_start).count() / 1000.0;

    // Copy alpha (estimated counts) to host
    std::vector<double> alpha(num_trans);
    thrust::copy(em.d_alpha.begin(), em.d_alpha.end(), alpha.begin());

    // Write output
    if (!opt.output.empty()) {
      struct stat st;
      if (stat(opt.output.c_str(), &st) != 0) {
        #ifdef _WIN32
          _mkdir(opt.output.c_str());
        #else
          mkdir(opt.output.c_str(), 0777);
        #endif
      }

      write_ec_counts_with_transcripts(opt.output + "/counts_gpu.txt",
                                       ec_counter, gpu_ecmap);

      if (new_ec_handler.total_new_ecs() > 0) {
        std::cout << "  New ECs from paired-end: " << new_ec_handler.total_new_ecs()
                  << " (covering " << new_ec_handler.total_new_reads() << " read pairs)" << std::endl;
      }

      plaintext_writer(opt.output + "/abundance.tsv",
                       index.target_names_, alpha, eff_lens, index.target_lens_);
    }
  }

  fill_run_stats(total_processed, ec_counter, gpu_ecmap, out_stats);

  auto wall_end = std::chrono::high_resolution_clock::now();
  g_benchmark_stats.wall_clock_total_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start).count() / 1000.0;
  g_benchmark_stats.wall_clock_pipeline_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(pipeline_end - pipeline_start).count() / 1000.0;

  print_benchmark_summary();
}

// Helper: convert Roaring bitmap to sorted vector of ints
static std::vector<int> roaring_to_sorted_vec(const Roaring& r) {
  std::vector<int> v(r.cardinality());
  if (!v.empty()) {
    std::vector<uint32_t> tmp(r.cardinality());
    r.toUint32Array(tmp.data());
    for (size_t i = 0; i < tmp.size(); ++i) {
      v[i] = static_cast<int>(tmp[i]);
    }
  }
  return v;
}

// Helper: format transcript vector as string (for file output)
static std::string format_tx_vec(const std::vector<int>& v, size_t max_show = 100) {
  std::ostringstream os;
  os << "{";
  for (size_t j = 0; j < v.size(); ++j) {
    if (j > 0) os << ",";
    if (j >= max_show) { os << "...(" << v.size() << " total)"; break; }
    os << v[j];
  }
  os << "}";
  return os.str();
}

// Look up a single k-mer in the CPU index (bifrost graph) and return its EC id.
// Uses the same approach as build_kmer_to_ec_map: find kmer in dbg, get transcript
// set from Node::ec[um.dist], look up integer EC id via ecmapinv.
// Returns -1 if not found in graph, -2 if found but no ecmapinv entry.
static int cpu_kmer_ec_lookup(KmerIndex& index, const Kmer& km) {
  const_UnitigMap<Node> um = index.dbg.find(km);
  if (um.isEmpty) return -1;
  const Roaring& trs = um.getData()->ec[um.dist].getIndices();
  auto ec_it = index.ecmapinv.find(trs);
  return (ec_it != index.ecmapinv.end()) ? ec_it->second : -2;
}

// Get transcript set for an EC id from the GPU ecmap host copy.
static std::vector<int> ec_transcripts(int ec_id,
    const std::vector<int>& h_gpu_transcripts,
    const std::vector<uint64_t>& h_gpu_offsets,
    size_t gpu_num_ecs)
{
  if (ec_id < 0 || ec_id >= (int)gpu_num_ecs) return {};
  uint64_t s = h_gpu_offsets[ec_id];
  uint64_t e = h_gpu_offsets[ec_id + 1];
  std::vector<int> v(h_gpu_transcripts.begin() + s, h_gpu_transcripts.begin() + e);
  std::sort(v.begin(), v.end());
  return v;
}

// Naive CPU intersection: look up every k-mer individually in the bifrost graph
// (no skipping/jumping like index.match). Returns the transcript intersection as
// a Roaring bitmap. Sets had_mapped_kmer to true if at least one k-mer was found.
static Roaring naive_cpu_read_intersect(KmerIndex& index, const std::string& seq, bool& had_mapped_kmer) {
  int k = index.k;
  int L = static_cast<int>(seq.size());
  int nkmers = (L >= k) ? (L - k + 1) : 0;
  Roaring running;
  bool first = true;
  had_mapped_kmer = false;

  for (int p = 0; p < nkmers; p++) {
    bool has_n = false;
    for (int j = 0; j < k; j++) {
      char c = std::toupper(seq[p + j]);
      if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { has_n = true; break; }
    }
    if (has_n) continue;

    Kmer km(seq.c_str() + p);
    const_UnitigMap<Node> um = index.dbg.find(km);
    if (um.isEmpty) continue;
    const Roaring& trs = um.getData()->ec[um.dist].getIndices();
    if (trs.isEmpty()) continue;

    had_mapped_kmer = true;
    if (first) {
      running = trs;
      first = false;
    } else {
      running &= trs;
    }
    if (running.isEmpty()) break;
  }
  return running;
}

// Naive CPU pair intersection: mirrors the GPU's updated pair intersection logic.
// Falls back to the other mate only if the empty mate had no mapped k-mers at all.
static Roaring naive_cpu_pair_intersect(
    const Roaring& r1, bool r1_had_mapped,
    const Roaring& r2, bool r2_had_mapped)
{
  if (r1.isEmpty() && r2.isEmpty()) return {};

  if (r1.isEmpty()) {
    if (r1_had_mapped) return {};
    return r2;
  }
  if (r2.isEmpty()) {
    if (r2_had_mapped) return {};
    return r1;
  }
  return r1 & r2;
}

// Dump detailed per-k-mer analysis for a single read to the output stream.
static void dump_read_kmer_detail(
    std::ostream& out,
    KmerIndex& index,
    const std::string& seq,
    const std::string& label,
    uint32_t read_idx,
    const std::vector<uint64_t>* h_kmers_ptr,
    const std::vector<uint64_t>* h_kmer_first_ptr,
    const std::vector<uint32_t>* h_kmer_count_ptr,
    uint64_t empty_kmer_value,
    const std::vector<int>* h_gpu_ecs_ptr,
    const std::vector<int>& h_gpu_transcripts,
    const std::vector<uint64_t>& h_gpu_offsets,
    size_t gpu_num_ecs,
    const std::vector<int>& gpu_tx_result,
    const std::vector<int>& cpu_tx_result,
    const std::vector<std::pair<const_UnitigMap<Node>, int>>& cpu_matches)
{
  int k = index.k;
  int L = static_cast<int>(seq.size());
  int nkmers = (L >= k) ? (L - k + 1) : 0;
  bool have_gpu = (h_kmers_ptr && h_kmer_first_ptr && h_kmer_count_ptr && h_gpu_ecs_ptr &&
                   read_idx < h_kmer_count_ptr->size());

  out << "\n========== " << label << " DETAILED K-MER ANALYSIS ==========\n";
  out << "Sequence (" << L << " bp): " << seq << "\n\n";

  // --- Section 1: Diagonal k-mer view with EC for each position ---
  // Format: indent by position, show kmer, show GPU ec and CPU ec
  out << "--- K-mer diagonal (k=" << k << ", " << nkmers << " positions) ---\n";
  out << seq << "\n";

  // Collect per-position data for later use
  struct KmerInfo {
    int gpu_ec = -1;
    int cpu_ec = -1;     // from raw dbg lookup (same as GPU map construction)
    bool gpu_empty = false;
    bool has_n = false;
    bool cpu_in_match = false;  // was this position visited by index.match?
    int match_ec = -1;          // ec from index.match at this position
  };
  std::vector<KmerInfo> kinfo(nkmers);

  // Build set of positions visited by index.match
  for (auto& [um, pos] : cpu_matches) {
    if (pos >= 0 && pos < nkmers) {
      const Roaring& trs = um.getData()->ec[um.dist].getIndices();
      auto ec_it = index.ecmapinv.find(trs);
      kinfo[pos].cpu_in_match = true;
      kinfo[pos].match_ec = (ec_it != index.ecmapinv.end()) ? ec_it->second : -2;
    }
  }

  // Fill GPU and CPU ec for each position
  for (int p = 0; p < nkmers; p++) {
    auto& ki = kinfo[p];
    // Check for N
    for (int j = 0; j < k; j++) {
      char c = std::toupper(seq[p + j]);
      if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { ki.has_n = true; break; }
    }

    // GPU
    if (have_gpu) {
      uint64_t start = (*h_kmer_first_ptr)[read_idx];
      uint32_t nk = (*h_kmer_count_ptr)[read_idx];
      if (p < (int)nk && start + p < h_kmers_ptr->size()) {
        uint64_t kval = (*h_kmers_ptr)[start + p];
        ki.gpu_empty = (kval == empty_kmer_value);
        if (!ki.gpu_empty && start + p < h_gpu_ecs_ptr->size()) {
          ki.gpu_ec = (*h_gpu_ecs_ptr)[start + p];
        }
      }
    }

    // CPU: raw kmer lookup in the graph (same method as GPU map construction)
    if (!ki.has_n) {
      Kmer km(seq.c_str() + p);
      ki.cpu_ec = cpu_kmer_ec_lookup(index, km);
    }
  }

  // Print diagonal
  for (int p = 0; p < nkmers; p++) {
    auto& ki = kinfo[p];
    for (int s = 0; s < p; s++) out << ' ';
    out << seq.substr(p, k);

    // GPU ec
    if (ki.gpu_empty)       out << "  gpu:*";
    else if (ki.gpu_ec < 0) out << "  gpu:*";
    else                    out << "  gpu:" << ki.gpu_ec;

    // CPU raw ec
    if (ki.has_n)           out << "  cpu:N";
    else if (ki.cpu_ec == -1) out << "  cpu:*";
    else if (ki.cpu_ec == -2) out << "  cpu:?";
    else                    out << "  cpu:" << ki.cpu_ec;

    // Was this position visited by index.match?
    if (ki.cpu_in_match)    out << "  match:" << ki.match_ec;

    // Flag divergence
    if (!ki.has_n && !ki.gpu_empty && ki.gpu_ec >= 0 && ki.cpu_ec >= 0 && ki.gpu_ec != ki.cpu_ec)
      out << "  !!DIFF!!";

    out << "\n";
  }

  // --- Section 2: Unique ECs and their transcript sets ---
  out << "\n--- Equivalence classes for " << label << " ---\n";
  std::set<int> gpu_ecs, cpu_ecs, match_ecs;
  for (int p = 0; p < nkmers; p++) {
    if (kinfo[p].gpu_ec >= 0) gpu_ecs.insert(kinfo[p].gpu_ec);
    if (kinfo[p].cpu_ec >= 0) cpu_ecs.insert(kinfo[p].cpu_ec);
    if (kinfo[p].match_ec >= 0) match_ecs.insert(kinfo[p].match_ec);
  }
  std::set<int> all_ecs;
  all_ecs.insert(gpu_ecs.begin(), gpu_ecs.end());
  all_ecs.insert(cpu_ecs.begin(), cpu_ecs.end());
  all_ecs.insert(match_ecs.begin(), match_ecs.end());
  for (int ec : all_ecs) {
    auto tx = ec_transcripts(ec, h_gpu_transcripts, h_gpu_offsets, gpu_num_ecs);
    out << "  EC " << ec;
    if (gpu_ecs.count(ec)) out << " [gpu]";
    if (cpu_ecs.count(ec)) out << " [cpu]";
    if (match_ecs.count(ec)) out << " [match]";
    out << ": {";
    for (size_t j = 0; j < tx.size(); j++) {
      if (j) out << ",";
      out << tx[j];
    }
    out << "}\n";
  }

  // --- Section 3: Running intersections ---
  auto print_intersection_walk = [&](const std::string& method,
      std::function<int(int)> get_ec, bool skip_n, bool skip_empty) {
    out << "\n--- " << method << " intersection for " << label << " ---\n";
    Roaring running;
    bool first = true;
    int step = 0;
    for (int p = 0; p < nkmers; p++) {
      if (skip_n && kinfo[p].has_n) continue;
      int ec = get_ec(p);
      if (ec < 0) continue;
      if (skip_empty && kinfo[p].gpu_empty) continue;
      auto tx = ec_transcripts(ec, h_gpu_transcripts, h_gpu_offsets, gpu_num_ecs);
      Roaring tx_set;
      for (int t : tx) tx_set.add(static_cast<uint32_t>(t));
      if (first) {
        running = tx_set;
        first = false;
      } else {
        running &= tx_set;
      }
      out << "  pos=" << p << " ec=" << ec
          << " |ec|=" << tx.size()
          << " |intersection|=" << running.cardinality() << "\n";
      if (running.isEmpty()) {
        out << "  >>> intersection became EMPTY at pos " << p << " <<<\n";
        break;
      }
      step++;
    }
    auto final_tx = roaring_to_sorted_vec(running);
    auto it = index.ecmapinv.find(running);
    int final_ec = (it != index.ecmapinv.end()) ? it->second : -1;
    out << "  RESULT: ec=" << final_ec << " transcripts={";
    for (size_t j = 0; j < final_tx.size(); j++) {
      if (j) out << ",";
      out << final_tx[j];
    }
    out << "}\n";
  };

  // GPU: use GPU ec values per position
  print_intersection_walk("GPU per-kmer",
      [&](int p) { return kinfo[p].gpu_ec; }, false, true);

  // CPU raw: look up every kmer individually (same as GPU map build)
  print_intersection_walk("CPU raw per-kmer (all kmers)",
      [&](int p) { return kinfo[p].cpu_ec; }, true, false);

  // CPU index.match: only positions visited by the match algorithm
  print_intersection_walk("CPU index.match (jumping/skipping)",
      [&](int p) { return kinfo[p].cpu_in_match ? kinfo[p].match_ec : -1; }, false, false);

  // Final summary
  out << "\n--- " << label << " final comparison ---\n";
  out << "  GPU pipeline result:     " << format_tx_vec(gpu_tx_result) << "\n";
  out << "  CPU index.match result:  " << format_tx_vec(cpu_tx_result) << "\n";
  out << "  Match: " << (gpu_tx_result == cpu_tx_result ? "YES" : "NO") << "\n";
}

// Core comparison logic: given host-side GPU pipeline results and host reads, compare with CPU.
// If mismatch_file is non-empty, all mismatches are written there (one block per read) and we don't stop early.
static int compare_gpu_vs_cpu(
    KmerIndex& index,
    const ProgramOptions& opt,
    bool is_paired,
    uint32_t r1_count,
    uint32_t total_read_count,
    const std::vector<char>& h_reads,
    const std::vector<uint64_t>& h_read_offsets,
    const std::vector<uint32_t>& h_read_lengths,
    const std::vector<int>& h_read_transcripts,
    const std::vector<uint64_t>& h_read_tx_offsets,
    const std::vector<uint64_t>& h_read_tx_sizes,
    const std::vector<int>& h_pair_transcripts,
    const std::vector<uint64_t>& h_pair_tx_offsets,
    const std::vector<uint64_t>& h_pair_tx_sizes,
    const std::vector<int>& h_ecs_final,
    const std::string& mismatch_file,
    const std::vector<uint64_t>* h_kmers_ptr,
    const std::vector<uint64_t>* h_kmer_first_ptr,
    const std::vector<uint32_t>* h_kmer_count_ptr,
    uint64_t empty_kmer_value,
    const std::vector<int>* h_gpu_ecs_ptr,
    const std::vector<int>& h_gpu_transcripts,
    const std::vector<uint64_t>& h_gpu_offsets,
    size_t gpu_num_ecs)
{
  MinCollector mc(index, opt);

  uint32_t num_pairs = is_paired ? r1_count : total_read_count;
  int gpu_vs_naive_mismatch = 0;
  int gpu_vs_naive_match = 0;
  int gpu_vs_match_mismatch = 0;
  int match_only_diff = 0;  // differs from match-CPU but agrees with naive-CPU

  std::ofstream mismatch_out;
  if (!mismatch_file.empty()) {
    mismatch_out.open(mismatch_file);
    if (!mismatch_out) {
      std::cerr << "[debug-gpu] Cannot open mismatch file: " << mismatch_file << std::endl;
    }
  }

  std::cout << "[debug-gpu] Three-way comparison (GPU / naive-CPU / match-CPU) on "
            << num_pairs << (is_paired ? " read pairs" : " reads") << " ..." << std::endl;

  for (uint32_t i = 0; i < num_pairs; ++i) {
    uint32_t r1_idx = i;
    uint64_t r1_off = h_read_offsets[r1_idx];
    uint32_t r1_len = h_read_lengths[r1_idx];
    std::string r1_seq(h_reads.data() + r1_off, r1_len);

    std::string r2_seq;
    uint32_t r2_idx = 0;
    uint32_t r2_len = 0;
    if (is_paired) {
      r2_idx = r1_count + i;
      uint64_t r2_off = h_read_offsets[r2_idx];
      r2_len = h_read_lengths[r2_idx];
      r2_seq.assign(h_reads.data() + r2_off, r2_len);
    }

    // --- GPU results ---
    std::vector<int> gpu_r1_tx;
    if (r1_idx < h_read_tx_sizes.size()) {
      uint64_t off = h_read_tx_offsets[r1_idx];
      uint64_t sz = h_read_tx_sizes[r1_idx];
      gpu_r1_tx.assign(h_read_transcripts.begin() + off,
                       h_read_transcripts.begin() + off + sz);
      std::sort(gpu_r1_tx.begin(), gpu_r1_tx.end());
    }

    std::vector<int> gpu_r2_tx;
    if (is_paired && r2_idx < h_read_tx_sizes.size()) {
      uint64_t off = h_read_tx_offsets[r2_idx];
      uint64_t sz = h_read_tx_sizes[r2_idx];
      gpu_r2_tx.assign(h_read_transcripts.begin() + off,
                       h_read_transcripts.begin() + off + sz);
      std::sort(gpu_r2_tx.begin(), gpu_r2_tx.end());
    }

    std::vector<int> gpu_pair_tx;
    if (is_paired && i < h_pair_tx_sizes.size()) {
      uint64_t off = h_pair_tx_offsets[i];
      uint64_t sz = h_pair_tx_sizes[i];
      gpu_pair_tx.assign(h_pair_transcripts.begin() + off,
                         h_pair_transcripts.begin() + off + sz);
      std::sort(gpu_pair_tx.begin(), gpu_pair_tx.end());
    }

    // --- Naive CPU: check every k-mer individually (no skipping) ---
    bool naive_r1_had_mapped = false;
    Roaring naive_r1_ec = naive_cpu_read_intersect(index, r1_seq, naive_r1_had_mapped);
    std::vector<int> naive_r1_tx = roaring_to_sorted_vec(naive_r1_ec);

    bool naive_r2_had_mapped = false;
    Roaring naive_r2_ec;
    std::vector<int> naive_r2_tx;
    if (is_paired) {
      naive_r2_ec = naive_cpu_read_intersect(index, r2_seq, naive_r2_had_mapped);
      naive_r2_tx = roaring_to_sorted_vec(naive_r2_ec);
    }

    std::vector<int> naive_pair_tx;
    Roaring naive_pair_ec;
    if (is_paired) {
      naive_pair_ec = naive_cpu_pair_intersect(naive_r1_ec, naive_r1_had_mapped,
                                               naive_r2_ec, naive_r2_had_mapped);
      naive_pair_tx = roaring_to_sorted_vec(naive_pair_ec);
    }

    // --- Match CPU: index.match with skipping ---
    std::vector<std::pair<const_UnitigMap<Node>, int>> v1, v2;
    index.match(r1_seq.c_str(), r1_len, v1);
    if (is_paired) {
      index.match(r2_seq.c_str(), r2_len, v2);
    }

    Roaring match_r1_ec = mc.intersectECs(v1);
    std::vector<int> match_r1_tx = roaring_to_sorted_vec(match_r1_ec);

    std::vector<int> match_r2_tx;
    Roaring match_pair_ec;
    std::vector<int> match_pair_tx;
    if (is_paired) {
      Roaring match_r2_ec = mc.intersectECs(v2);
      match_r2_tx = roaring_to_sorted_vec(match_r2_ec);
      mc.intersectKmers(v1, v2, false, match_pair_ec);
      match_pair_tx = roaring_to_sorted_vec(match_pair_ec);
    }

    // --- Compare GPU vs naive CPU (the primary check) ---
    bool naive_r1_ok = (gpu_r1_tx == naive_r1_tx);
    bool naive_r2_ok = !is_paired || (gpu_r2_tx == naive_r2_tx);
    bool naive_pair_ok = !is_paired || (gpu_pair_tx == naive_pair_tx);
    bool gpu_naive_agree = naive_r1_ok && naive_r2_ok && naive_pair_ok;

    // --- Compare GPU vs match CPU (for statistics) ---
    bool match_pair_ok = !is_paired || (gpu_pair_tx == match_pair_tx);

    if (gpu_naive_agree) {
      gpu_vs_naive_match++;
      if (!match_pair_ok) match_only_diff++;
    } else {
      gpu_vs_naive_mismatch++;
    }
    if (!match_pair_ok) gpu_vs_match_mismatch++;

    // Only dump details for GPU vs naive-CPU mismatches
    if (!gpu_naive_agree) {
      std::ostream& out = mismatch_out.is_open() ? (std::ostream&)mismatch_out : std::cout;

      out << "\n################################################################\n";
      out << "### READ PAIR " << i << " [GPU != naive-CPU] ###\n";
      out << "################################################################\n";
      out << "R1_seq (" << r1_len << " bp): " << r1_seq << "\n";
      if (is_paired) out << "R2_seq (" << r2_len << " bp): " << r2_seq << "\n";

      if (opt.verbose) {
        dump_read_kmer_detail(out, index, r1_seq, "R1", r1_idx,
                              h_kmers_ptr, h_kmer_first_ptr, h_kmer_count_ptr,
                              empty_kmer_value, h_gpu_ecs_ptr,
                              h_gpu_transcripts, h_gpu_offsets, gpu_num_ecs,
                              gpu_r1_tx, match_r1_tx, v1);
        if (is_paired) {
          dump_read_kmer_detail(out, index, r2_seq, "R2", r2_idx,
                                h_kmers_ptr, h_kmer_first_ptr, h_kmer_count_ptr,
                                empty_kmer_value, h_gpu_ecs_ptr,
                                h_gpu_transcripts, h_gpu_offsets, gpu_num_ecs,
                                gpu_r2_tx, match_r2_tx, v2);
        }
      }

      out << "\n--- COMPARISON for pair " << i << " ---\n";
      out << "  GPU  R1: " << format_tx_vec(gpu_r1_tx) << "\n";
      out << "  Naive R1: " << format_tx_vec(naive_r1_tx)
          << (naive_r1_had_mapped ? "" : " (no k-mers mapped)") << "\n";
      out << "  Match R1: " << format_tx_vec(match_r1_tx) << "\n";
      if (!naive_r1_ok) out << "  >>> R1 DIFFERS: GPU vs naive-CPU <<<\n";

      if (is_paired) {
        out << "  GPU  R2: " << format_tx_vec(gpu_r2_tx) << "\n";
        out << "  Naive R2: " << format_tx_vec(naive_r2_tx)
            << (naive_r2_had_mapped ? "" : " (no k-mers mapped)") << "\n";
        out << "  Match R2: " << format_tx_vec(match_r2_tx) << "\n";
        if (!naive_r2_ok) out << "  >>> R2 DIFFERS: GPU vs naive-CPU <<<\n";

        int gpu_final_ec = (i < h_ecs_final.size()) ? h_ecs_final[i] : -1;
        auto naive_pair_ec_it = index.ecmapinv.find(naive_pair_ec);
        int naive_pair_ec_id = (naive_pair_ec_it != index.ecmapinv.end()) ? naive_pair_ec_it->second : -1;

        out << "  GPU  pair: ec=" << gpu_final_ec << " " << format_tx_vec(gpu_pair_tx) << "\n";
        out << "  Naive pair: ec=" << naive_pair_ec_id << " " << format_tx_vec(naive_pair_tx) << "\n";
        out << "  Match pair: " << format_tx_vec(match_pair_tx) << "\n";
        if (!naive_pair_ok) out << "  >>> PAIR DIFFERS: GPU vs naive-CPU <<<\n";
      }
    }
  }

  std::cout << "\n[debug-gpu] Three-way comparison results (" << num_pairs
            << (is_paired ? " pairs" : " reads") << "):\n"
            << "  GPU vs naive-CPU: " << gpu_vs_naive_match << " agree, "
            << gpu_vs_naive_mismatch << " disagree\n"
            << "  GPU vs match-CPU: " << (num_pairs - gpu_vs_match_mismatch) << " agree, "
            << gpu_vs_match_mismatch << " disagree\n"
            << "  match-only diffs (GPU==naive but GPU!=match): " << match_only_diff << "\n";
  if (gpu_vs_naive_mismatch == 0) {
    std::cout << "[debug-gpu] GPU and naive-CPU agree on all reads! All differences are due to index.match skipping.\n";
  }
  if (gpu_vs_naive_mismatch > 0 && mismatch_out.is_open()) {
    std::cout << "[debug-gpu] Wrote " << gpu_vs_naive_mismatch
              << " GPU-vs-naive mismatches to: " << mismatch_file << "\n";
  }
  return gpu_vs_naive_mismatch;
}

// Copy GPU pipeline results to host vectors for comparison
static void copy_pipeline_results_to_host(
    const ReadTranscriptIntersector& intersector,
    const PairIntersector& pair_intersector,
    const ReadECLookup& ec_lookup,
    bool is_paired,
    std::vector<int>& h_read_transcripts,
    std::vector<uint64_t>& h_read_tx_offsets,
    std::vector<uint64_t>& h_read_tx_sizes,
    std::vector<int>& h_pair_transcripts,
    std::vector<uint64_t>& h_pair_tx_offsets,
    std::vector<uint64_t>& h_pair_tx_sizes,
    std::vector<int>& h_ecs_final)
{
  h_read_transcripts.resize(intersector.read_transcripts.size());
  thrust::copy(intersector.read_transcripts.begin(), intersector.read_transcripts.end(),
               h_read_transcripts.begin());
  h_read_tx_offsets.resize(intersector.read_transcript_offsets.size());
  thrust::copy(intersector.read_transcript_offsets.begin(), intersector.read_transcript_offsets.end(),
               h_read_tx_offsets.begin());
  h_read_tx_sizes.resize(intersector.read_transcript_sizes.size());
  thrust::copy(intersector.read_transcript_sizes.begin(), intersector.read_transcript_sizes.end(),
               h_read_tx_sizes.begin());

  if (is_paired) {
    h_pair_transcripts.resize(pair_intersector.pair_transcripts.size());
    thrust::copy(pair_intersector.pair_transcripts.begin(), pair_intersector.pair_transcripts.end(),
                 h_pair_transcripts.begin());
    h_pair_tx_offsets.resize(pair_intersector.pair_transcript_offsets.size());
    thrust::copy(pair_intersector.pair_transcript_offsets.begin(), pair_intersector.pair_transcript_offsets.end(),
                 h_pair_tx_offsets.begin());
    h_pair_tx_sizes.resize(pair_intersector.pair_transcript_sizes.size());
    thrust::copy(pair_intersector.pair_transcript_sizes.begin(), pair_intersector.pair_transcript_sizes.end(),
                 h_pair_tx_sizes.begin());
  }

  h_ecs_final.resize(ec_lookup.read_ecs_final.size());
  thrust::copy(ec_lookup.read_ecs_final.begin(), ec_lookup.read_ecs_final.end(),
               h_ecs_final.begin());
}

// Load EC counts from counts_gpu.txt into host vectors.
// Format: ec_id\tcount\ttranscript1,transcript2,...
static bool load_counts_from_file(
    const std::string& filename,
    std::vector<int>& h_counts,
    std::vector<int>& h_transcripts,
    std::vector<uint64_t>& h_offsets)
{
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "[debug-gpu] Cannot open " << filename << std::endl;
    return false;
  }

  struct ECEntry {
    int ec_id;
    int count;
    std::vector<int> transcripts;
  };
  std::vector<ECEntry> entries;

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    ECEntry e;
    std::istringstream iss(line);
    std::string count_str, tx_str;
    iss >> e.ec_id >> e.count;
    std::getline(iss >> std::ws, tx_str);
    if (!tx_str.empty()) {
      std::istringstream txss(tx_str);
      std::string tok;
      while (std::getline(txss, tok, ',')) {
        if (!tok.empty()) e.transcripts.push_back(std::stoi(tok));
      }
    }
    entries.push_back(std::move(e));
  }
  infile.close();

  size_t num_ecs = entries.size();
  h_counts.resize(num_ecs);
  h_offsets.resize(num_ecs + 1);
  h_transcripts.clear();

  uint64_t offset = 0;
  for (size_t i = 0; i < num_ecs; i++) {
    h_counts[i] = entries[i].count;
    h_offsets[i] = offset;
    for (int t : entries[i].transcripts) {
      h_transcripts.push_back(t);
    }
    offset += entries[i].transcripts.size();
  }
  h_offsets[num_ecs] = offset;

  std::cout << "[debug-gpu] Loaded " << num_ecs << " ECs from " << filename << std::endl;
  return true;
}

// Build a map from sorted transcript list -> count for GPU EC counts.
static std::map<std::vector<int>, int> build_transcript_count_map_gpu(
    const std::vector<int>& h_counts,
    const std::vector<int>& h_transcripts,
    const std::vector<uint64_t>& h_offsets,
    size_t num_ecs)
{
  std::map<std::vector<int>, int> result;
  for (size_t ec = 0; ec < num_ecs; ec++) {
    if (h_counts[ec] == 0) continue;
    uint64_t start = h_offsets[ec];
    uint64_t end = h_offsets[ec + 1];
    std::vector<int> tx_list(h_transcripts.begin() + start, h_transcripts.begin() + end);
    std::sort(tx_list.begin(), tx_list.end());
    result[tx_list] += h_counts[ec];
  }
  return result;
}

// Build a map from sorted transcript list -> count for CPU EC counts.
static std::map<std::vector<int>, int> build_transcript_count_map_cpu(
    const MinCollector& tc,
    const KmerIndex& index)
{
  std::map<std::vector<int>, int> result;
  for (const auto& it : index.ecmapinv) {
    int ec_id = it.second;
    if (ec_id < 0 || ec_id >= (int)tc.counts.size()) continue;
    if (tc.counts[ec_id] == 0) continue;
    std::vector<int> tx_list;
    for (uint32_t t : it.first) tx_list.push_back(static_cast<int>(t));
    std::sort(tx_list.begin(), tx_list.end());
    result[tx_list] += tc.counts[ec_id];
  }
  return result;
}



void print_benchmark_summary() {
  const BenchmarkStats& stats = g_benchmark_stats;
  
  double total_setup = stats.setup_index_load_ms + stats.setup_kmer_to_ec_map_ms +
                       stats.setup_gpu_ecmap_ms + stats.setup_gpu_ecmapinv_ms;
  double total_gpu = stats.gpu_kmer_extraction_ms + stats.gpu_kmer_lookup_ms +
                     stats.gpu_ec_collapse_ms + stats.gpu_transcript_intersection_ms +
                     stats.gpu_ec_lookup_ms + stats.gpu_ec_counting_ms;
  double total_host = stats.host_h2d_copy_ms + stats.host_thrust_ops_ms + 
                      stats.host_file_write_ms;
  double component_sum = total_setup + total_gpu + total_host + stats.io_decompress_ms + stats.gpu_em_ms;
  double wall_total = stats.wall_clock_total_ms;
  double wall_pipeline = stats.wall_clock_pipeline_ms;

  // Use wall-clock time for percentages when available
  double ref_time = (wall_total > 0) ? wall_total : component_sum;
  
  std::cout << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << "GPU Pipeline Benchmark Summary" << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::fixed << std::setprecision(2);
  
  std::cout << std::endl << "Setup Operations:" << std::endl;
  std::cout << "  Index Load:             " << std::setw(10) << stats.setup_index_load_ms << " ms" << std::endl;
  std::cout << "  Kmer-to-EC Map Build:   " << std::setw(10) << stats.setup_kmer_to_ec_map_ms << " ms" << std::endl;
  std::cout << "  GPU EC Map Build:        " << std::setw(10) << stats.setup_gpu_ecmap_ms << " ms" << std::endl;
  std::cout << "  GPU EC Map Inv Build:    " << std::setw(10) << stats.setup_gpu_ecmapinv_ms << " ms" << std::endl;
  std::cout << "  Total Setup:             " << std::setw(10) << total_setup << " ms" << std::endl;
  
  std::cout << std::endl << "Pipeline (" << stats.batch_count << " batches):" << std::endl;
  std::cout << "  GPU Kmer Extraction:    " << std::setw(10) << stats.gpu_kmer_extraction_ms << " ms" << std::endl;
  std::cout << "  GPU Kmer Lookup:        " << std::setw(10) << stats.gpu_kmer_lookup_ms << " ms" << std::endl;
  std::cout << "  GPU EC Collapse:        " << std::setw(10) << stats.gpu_ec_collapse_ms << " ms" << std::endl;
  std::cout << "  GPU Transcript Intersect:" << std::setw(9) << stats.gpu_transcript_intersection_ms << " ms" << std::endl;
  std::cout << "  GPU EC Lookup:          " << std::setw(10) << stats.gpu_ec_lookup_ms << " ms" << std::endl;
  std::cout << "  GPU EC Counting:        " << std::setw(10) << stats.gpu_ec_counting_ms << " ms" << std::endl;
  std::cout << "  Total GPU compute:      " << std::setw(10) << total_gpu << " ms" << std::endl;
  std::cout << "  I/O & Decompression:    " << std::setw(10) << stats.io_decompress_ms << " ms" << std::endl;
  std::cout << "  D2D Copy:               " << std::setw(10) << stats.host_h2d_copy_ms << " ms" << std::endl;
  std::cout << "  Host Thrust Operations: " << std::setw(10) << stats.host_thrust_ops_ms << " ms" << std::endl;
  double pipeline_overlap = (total_gpu + stats.io_decompress_ms + stats.host_h2d_copy_ms + stats.host_thrust_ops_ms) - wall_pipeline;
  std::cout << "  Pipeline wall-clock:    " << std::setw(10) << wall_pipeline << " ms" << std::endl;
  if (pipeline_overlap > 0) {
    std::cout << "  Overlap saved:          " << std::setw(10) << pipeline_overlap << " ms" << std::endl;
  }

  std::cout << std::endl << "EM Algorithm:             " << std::setw(10) << stats.gpu_em_ms << " ms" << std::endl;

  std::cout << std::endl << "File Writing:             " << std::setw(10) << stats.host_file_write_ms << " ms" << std::endl;

  std::cout << std::endl << "Wall-clock total:         " << std::setw(10) << wall_total << " ms" << std::endl;
  std::cout << "  (Component sum:         " << std::setw(10) << component_sum << " ms)" << std::endl;
  std::cout << "=========================================" << std::endl;
}
