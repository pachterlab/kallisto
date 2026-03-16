#ifndef BENCHMARK_STATS_H
#define BENCHMARK_STATS_H

// Benchmark statistics structure - shared between GPU pipeline components
struct BenchmarkStats {
  // Setup timings (milliseconds)
  double setup_index_load_ms = 0.0;
  double setup_kmer_to_ec_map_ms = 0.0;
  double setup_gpu_ecmap_ms = 0.0;
  double setup_gpu_ecmapinv_ms = 0.0;
  
  // GPU timings (milliseconds)
  double gpu_kmer_extraction_ms = 0.0;
  double gpu_kmer_lookup_ms = 0.0;
  double gpu_ec_collapse_ms = 0.0;
  double gpu_transcript_intersection_ms = 0.0;
  double gpu_ec_lookup_ms = 0.0;
  double gpu_ec_counting_ms = 0.0;
  double gpu_em_ms = 0.0;
  
  // Host timings (milliseconds)
  double host_h2d_copy_ms = 0.0;
  double host_thrust_ops_ms = 0.0;
  double host_file_write_ms = 0.0;
  
  // I/O timings (milliseconds)
  double io_decompress_ms = 0.0;

  // Wall-clock timings (actual elapsed time, accounts for overlap)
  double wall_clock_total_ms = 0.0;
  double wall_clock_pipeline_ms = 0.0;
  int batch_count = 0;

  // Kmer count tracking
  uint64_t total_kmers = 0;
};

// Global benchmark statistics (defined in GPUProcessReads.cu)
extern BenchmarkStats g_benchmark_stats;

#endif // BENCHMARK_STATS_H
