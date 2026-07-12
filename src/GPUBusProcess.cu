#include "GPUBusProcess.cuh"

#include <cuda_runtime.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>

#include <chrono>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "BUSData.h"
#include "BUSTools.h"
#include "GPUIndex.cuh"
#include "GPUKernels.cuh"
#include "GPUPipeline.cuh"
#include "GPUReadLoader.cuh"
#include "PlaintextWriter.h"

// Forward decl from GPUProcessReads.cu
void print_benchmark_summary();

// 2-bit encoding of a DNA segment matching the CPU stringToBinary in BUSData.cpp.
// Encodes seq[start .. start+span) and writes flag bits in flag_out.
__device__ inline uint64_t encode_segment(const char* seq, int start, int span,
                                          uint32_t& flag_out) {
  flag_out = 0;
  uint64_t r = 0;
  int numN = 0;
  int posN = 0;
  int k = span;
  if (k > 32) k = 32;
  for (int i = 0; i < k; ++i) {
    unsigned char b = static_cast<unsigned char>(seq[start + i]);
    uint64_t x = ((uint64_t)(b & 4)) >> 1;
    if ((b & 3) == 2) {
      if (numN == 0) posN = i;
      ++numN;
    }
    r = (r << 2) | (x + ((x ^ (uint64_t)(b & 2)) >> 1));
  }
  if (numN > 0) {
    if (numN > 3) numN = 3;
    flag_out = (uint32_t)((numN & 3) | ((posN & 31) << 2));
  }
  return r;
}

__global__ void extract_bc_umi_kernel(
    const char* __restrict__ d_reads,
    const uint64_t* __restrict__ d_offsets,
    const uint32_t* __restrict__ d_lengths,
    uint64_t* __restrict__ d_bc,
    uint64_t* __restrict__ d_umi,
    uint32_t* __restrict__ d_flags,
    uint8_t* __restrict__ d_bad,
    int bc_start, int bc_span,
    int umi_start, int umi_span,
    uint32_t num_pairs) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_pairs) return;
  uint64_t off = d_offsets[i];
  uint32_t len = d_lengths[i];
  bool bad = false;
  int need_len = (bc_start + bc_span > umi_start + umi_span)
                     ? bc_start + bc_span
                     : umi_start + umi_span;
  if ((int)len < need_len) bad = true;
  uint32_t bc_flag = 0;
  uint32_t umi_flag = 0;
  uint64_t bc_val = 0;
  uint64_t umi_val = 0;
  if (!bad) {
    const char* seq = d_reads + off;
    bc_val = encode_segment(seq, bc_start, bc_span, bc_flag);
    umi_val = encode_segment(seq, umi_start, umi_span, umi_flag);
  }
  d_bc[i] = bc_val;
  d_umi[i] = umi_val;
  d_flags[i] = (bc_flag & 0xff) | ((umi_flag & 0xff) << 8);
  d_bad[i] = bad ? 1 : 0;
}

__global__ void assemble_busdata_kernel(
    const uint64_t* __restrict__ d_bc,
    const uint64_t* __restrict__ d_umi,
    const uint32_t* __restrict__ d_flags,
    const uint8_t* __restrict__ d_bad,
    const int* __restrict__ d_ec_per_read,
    uint32_t r1_count,
    uint32_t num_pairs,
    BUSData* __restrict__ d_bus_out) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_pairs) return;
  uint32_t cdna_idx = r1_count + i;
  int ec = d_ec_per_read[cdna_idx];
  BUSData* b = d_bus_out + i;
  b->barcode = d_bc[i];
  b->UMI = d_umi[i];
  b->flags = d_flags[i];
  b->count = 1;
  b->pad = 0;
  b->ec = (ec < 0 || d_bad[i]) ? -1 : ec;
}

struct EcNonNegative {
  __device__ bool operator()(const BUSData& b) const { return b.ec >= 0; }
};

// Single batch of bus pipeline: kmer + lookup + collapse + intersect +
// EC lookup + handle_batch + count + bc/umi extract + assemble + compact + write.
// Returns number of compacted BUSData records emitted (already written to file).
static size_t process_bus_batch(
    DeviceKmerLoader& d_loader,
    thrust::device_vector<int>& d_ecs,
    cuco::static_map<uint64_t, int>& d_map,
    ReadECCollapser& collapser,
    ReadTranscriptIntersector& intersector,
    ReadECLookup& ec_lookup,
    ECCounter& ec_counter,
    NewECHandler& new_ec_handler,
    GPUECMap& gpu_ecmap,
    GPUECMapInv& gpu_ecmapinv,
    int bc_start, int bc_span, int umi_start, int umi_span,
    thrust::device_vector<uint64_t>& d_bc,
    thrust::device_vector<uint64_t>& d_umi,
    thrust::device_vector<uint32_t>& d_flags,
    thrust::device_vector<uint8_t>& d_bad,
    thrust::device_vector<BUSData>& d_bus,
    thrust::device_vector<BUSData>& d_bus_compact,
    BUSData* h_bus_pinned, size_t pinned_capacity,
    std::ofstream& bus_out) {
  d_loader.run();
  d_ecs.resize(d_loader.kmers.size(), -1);
  d_map.find(d_loader.kmers.begin(), d_loader.kmers.end(), d_ecs.begin());
  collapser.collapse(d_loader, d_ecs);
  intersector.intersect_gpu(collapser, gpu_ecmap);
  ec_lookup.lookup(intersector, gpu_ecmapinv, gpu_ecmap);
  new_ec_handler.handle_batch(intersector, ec_lookup, gpu_ecmapinv, gpu_ecmap,
                              ec_counter);
  ec_counter.count_batch(ec_lookup);

  uint32_t num_pairs = d_loader.r1_count;
  if (num_pairs == 0) return 0;

  d_bc.resize(num_pairs);
  d_umi.resize(num_pairs);
  d_flags.resize(num_pairs);
  d_bad.resize(num_pairs);
  d_bus.resize(num_pairs);
  d_bus_compact.resize(num_pairs);

  const int threads = 256;
  int blocks = (int)((num_pairs + threads - 1) / threads);

  extract_bc_umi_kernel<<<blocks, threads>>>(
      d_loader.reads.data().get(),
      d_loader.read_id_to_offset.data().get(),
      d_loader.read_id_to_length.data().get(),
      d_bc.data().get(), d_umi.data().get(), d_flags.data().get(),
      d_bad.data().get(),
      bc_start, bc_span, umi_start, umi_span, num_pairs);

  assemble_busdata_kernel<<<blocks, threads>>>(
      d_bc.data().get(), d_umi.data().get(), d_flags.data().get(),
      d_bad.data().get(),
      ec_lookup.read_ecs_final.data().get(),
      d_loader.r1_count, num_pairs,
      d_bus.data().get());

  auto compact_end = thrust::copy_if(
      thrust::device,
      d_bus.begin(), d_bus.begin() + num_pairs,
      d_bus_compact.begin(),
      EcNonNegative());
  size_t n_emit = compact_end - d_bus_compact.begin();
  if (n_emit == 0) return 0;

  if (n_emit > pinned_capacity) {
    // Fallback: chunked writes (shouldn't happen if pinned buffer >= num_pairs)
    std::vector<BUSData> h_bus(n_emit);
    thrust::copy(d_bus_compact.begin(), d_bus_compact.begin() + n_emit,
                 h_bus.begin());
    bus_out.write(reinterpret_cast<char*>(h_bus.data()),
                  n_emit * sizeof(BUSData));
  } else {
    cudaMemcpyAsync(h_bus_pinned, d_bus_compact.data().get(),
                    n_emit * sizeof(BUSData), cudaMemcpyDeviceToHost, 0);
    cudaStreamSynchronize(0);
    bus_out.write(reinterpret_cast<char*>(h_bus_pinned),
                  n_emit * sizeof(BUSData));
  }
  return n_emit;
}

// Write matrix.ec from GPU EC map (works for ECs minted on the GPU too).
static void write_ec_list_from_gpu_ecmap(const std::string& filename,
                                          const GPUECMap& gpu_ecmap) {
  std::vector<int> h_tx(gpu_ecmap.transcripts.size());
  thrust::copy(gpu_ecmap.transcripts.begin(), gpu_ecmap.transcripts.end(),
               h_tx.begin());
  std::vector<uint64_t> h_off(gpu_ecmap.offsets.size());
  thrust::copy(gpu_ecmap.offsets.begin(), gpu_ecmap.offsets.end(),
               h_off.begin());
  std::ofstream f(filename);
  for (size_t ec = 0; ec < gpu_ecmap.num_ecs; ++ec) {
    f << ec << "\t";
    uint64_t s = h_off[ec];
    uint64_t e = h_off[ec + 1];
    bool first = true;
    for (uint64_t i = s; i < e; ++i) {
      if (!first) f << ",";
      first = false;
      f << h_tx[i];
    }
    f << "\n";
  }
}

static void write_transcripts_file(const std::string& filename,
                                   const std::vector<std::string>& names,
                                   size_t count) {
  std::ofstream f(filename);
  for (size_t i = 0; i < count; ++i) {
    f << names[i] << "\n";
  }
}

static std::string make_call_string(int argc, char** argv) {
  std::string s;
  for (int i = 0; i < argc; ++i) {
    if (i > 0) s += " ";
    s += argv[i];
  }
  return s;
}

template <typename IndexT>
static cuco::static_map<uint64_t, int> build_kmer_map_for(const IndexT& index);

template <>
cuco::static_map<uint64_t, int> build_kmer_map_for<KmerIndex>(
    const KmerIndex& index) {
  return build_kmer_to_ec_map(index);
}

template <>
cuco::static_map<uint64_t, int> build_kmer_map_for<GPUIndex>(
    const GPUIndex& index) {
  return build_kmer_to_ec_map_from_contigs(index);
}

template <typename IndexT>
static size_t index_num_transcripts(const IndexT& index);

template <>
size_t index_num_transcripts<KmerIndex>(const KmerIndex& index) {
  return index.num_trans;
}

template <>
size_t index_num_transcripts<GPUIndex>(const GPUIndex& index) {
  return index.num_transcripts;
}

template <typename IndexT>
static int index_version_int(const IndexT& index);

template <>
int index_version_int<KmerIndex>(const KmerIndex& index) {
  return index.INDEX_VERSION;
}

template <>
int index_version_int<GPUIndex>(const GPUIndex& /*index*/) {
  return GPU_INDEX_FORMAT_VERSION;
}

template <typename IndexT>
static int index_k(const IndexT& index);

template <>
int index_k<KmerIndex>(const KmerIndex& index) { return index.k; }

template <>
int index_k<GPUIndex>(const GPUIndex& index) { return index.k; }

// Shared implementation for both KmerIndex and GPUIndex.
template <typename IndexT>
static void gpu_bus_run_impl(ProgramOptions& opt, const IndexT& index,
                             const std::string& start_time, int argc,
                             char** argv, GPURunStats* out_stats) {
  auto wall_start = std::chrono::high_resolution_clock::now();

  // Output directory
  if (!opt.output.empty()) {
    struct stat st;
    if (stat(opt.output.c_str(), &st) != 0) {
#ifdef _WIN32
      _mkdir(opt.output.c_str());
#else
      mkdir(opt.output.c_str(), 0777);
#endif
    }
  }

  // Open output.bus and write header eagerly: bclen and umilen are known from -x.
  std::ofstream busf_out(opt.output + "/output.bus", std::ios::binary);
  if (!busf_out.is_open()) {
    std::cerr << "Error: could not open " << opt.output << "/output.bus" << std::endl;
    exit(1);
  }
  int bclen = opt.busOptions.getBCLength();
  int umilen = opt.busOptions.getUMILength();
  writeBUSHeader(busf_out, bclen, umilen);

  // Compute BC/UMI offsets within file-0 read (v1: single piece each).
  const auto& bc = opt.busOptions.bc[0];
  const auto& umi = opt.busOptions.umi[0];
  int bc_start = bc.start;
  int bc_span = bc.stop - bc.start;
  int umi_start = umi.start;
  int umi_span = umi.stop - umi.start;

  GPUReadLoader gpu_loader(opt);

  DeviceKmerLoader d_loaders[2];
  int curr = 0;

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
        cv.wait(lock, [&] {
          return load_state == REQUESTED || load_state == SHUTDOWN;
        });
        if (load_state == SHUTDOWN) break;
      }
      bool ok = gpu_loader.load();
      if (ok) d_loaders[load_buf].load_from_gpu(gpu_loader);
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
    cv.wait(lock, [&] { return load_state == DONE; });
    load_state = IDLE;
    return load_has_data;
  };

  std::cout << "[Building GPU k-mer to EC map]" << std::endl;
  auto d_map = build_kmer_map_for<IndexT>(index);

  std::cout << "[Building GPU EC Map]" << std::endl;
  GPUECMap gpu_ecmap(index);

  std::cout << "[Building GPU EC Map Inverse]" << std::endl;
  GPUECMapInv gpu_ecmapinv(index);

  ECCounter ec_counter(index.ecmapinv.size());
  NewECHandler new_ec_handler(static_cast<int>(index.ecmapinv.size()));

  thrust::device_vector<int> d_ecs;
  ReadECCollapser collapser;
  ReadTranscriptIntersector intersector;
  ReadECLookup ec_lookup;

  thrust::device_vector<uint64_t> d_bc, d_umi;
  thrust::device_vector<uint32_t> d_flags;
  thrust::device_vector<uint8_t> d_bad;
  thrust::device_vector<BUSData> d_bus, d_bus_compact;

  // Pinned host buffer for streamed writes.
  size_t pinned_capacity = 1u << 20;  // 1M records ~32MB
  BUSData* h_bus_pinned = nullptr;
  cudaMallocHost(reinterpret_cast<void**>(&h_bus_pinned),
                 pinned_capacity * sizeof(BUSData));

  request_load(curr);

  bool has_data = wait_load();
  auto pipeline_start = std::chrono::high_resolution_clock::now();
  int64_t total_processed = 0;
  size_t total_emitted = 0;

  while (has_data) {
    int next = 1 - curr;
    request_load(next);

    size_t n_emit = process_bus_batch(
        d_loaders[curr], d_ecs, d_map, collapser, intersector, ec_lookup,
        ec_counter, new_ec_handler, gpu_ecmap, gpu_ecmapinv, bc_start, bc_span,
        umi_start, umi_span, d_bc, d_umi, d_flags, d_bad, d_bus, d_bus_compact,
        h_bus_pinned, pinned_capacity, busf_out);
    total_emitted += n_emit;

    if (d_loaders[curr].r1_count > 0 &&
        d_loaders[curr].read_count > d_loaders[curr].r1_count) {
      total_processed += static_cast<int64_t>(d_loaders[curr].r1_count);
    } else {
      total_processed += static_cast<int64_t>(d_loaders[curr].read_count);
    }
    g_benchmark_stats.batch_count++;

    has_data = wait_load();
    curr = next;
  }
  auto pipeline_end = std::chrono::high_resolution_clock::now();

  {
    std::lock_guard<std::mutex> lock(mtx);
    load_state = SHUTDOWN;
  }
  cv.notify_one();
  loader_thread.join();

  busf_out.close();
  if (h_bus_pinned) cudaFreeHost(h_bus_pinned);

  // Write matrix.ec, transcripts.txt, run_info.json
  write_ec_list_from_gpu_ecmap(opt.output + "/matrix.ec", gpu_ecmap);
  size_t num_transcripts = index_num_transcripts(index);
  write_transcripts_file(opt.output + "/transcripts.txt", index.target_names_,
                         num_transcripts);

  // run_info.json
  int64_t num_pseudoaligned = 0;
  int64_t num_unique = 0;
  {
    size_t num_ecs = ec_counter.ec_counts.size();
    std::vector<int> h_counts(num_ecs);
    thrust::copy(ec_counter.ec_counts.begin(), ec_counter.ec_counts.end(),
                 h_counts.begin());
    std::vector<uint64_t> h_off(gpu_ecmap.offsets.size());
    thrust::copy(gpu_ecmap.offsets.begin(), gpu_ecmap.offsets.end(),
                 h_off.begin());
    for (size_t ec = 0; ec < num_ecs; ++ec) {
      int c = h_counts[ec];
      num_pseudoaligned += c;
      uint64_t card = (ec + 1 < h_off.size()) ? h_off[ec + 1] - h_off[ec] : 0;
      if (card == 1) num_unique += c;
    }
  }

  std::string call = make_call_string(argc, argv);
  plaintext_aux(opt.output + "/run_info.json",
                std::to_string(num_transcripts),
                std::to_string(0),
                std::to_string(total_processed),
                std::to_string(num_pseudoaligned),
                std::to_string(num_unique),
                KALLISTO_VERSION,
                std::to_string(index_version_int(index)),
                std::to_string(index_k(index)),
                start_time, call);

  if (out_stats) {
    out_stats->num_processed = total_processed;
    out_stats->num_pseudoaligned = num_pseudoaligned;
    out_stats->num_unique = num_unique;
  }

  if (new_ec_handler.total_new_ecs() > 0) {
    std::cout << "  Novel ECs minted on GPU: " << new_ec_handler.total_new_ecs()
              << " (covering " << new_ec_handler.total_new_reads()
              << " reads)" << std::endl;
  }
  std::cout << "  BUS records written: " << total_emitted << std::endl;

  auto wall_end = std::chrono::high_resolution_clock::now();
  g_benchmark_stats.wall_clock_total_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start)
          .count() /
      1000.0;
  g_benchmark_stats.wall_clock_pipeline_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(pipeline_end -
                                                            pipeline_start)
          .count() /
      1000.0;

  print_benchmark_summary();
}

void gpu_bus_run(ProgramOptions& opt, const KmerIndex& index,
                 const std::string& start_time, int argc, char** argv,
                 GPURunStats* out_stats) {
  gpu_bus_run_impl<KmerIndex>(opt, index, start_time, argc, argv, out_stats);
}

void gpu_bus_run(ProgramOptions& opt, const GPUIndex& index,
                 const std::string& start_time, int argc, char** argv,
                 GPURunStats* out_stats) {
  gpu_bus_run_impl<GPUIndex>(opt, index, start_time, argc, argv, out_stats);
}
