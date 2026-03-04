#ifndef GPU_PROCESS_READS_CUH
#define GPU_PROCESS_READS_CUH

#include "common.h"
#include "KmerIndex.h"
#include "GPUIndexFormat.h"

#include <vector>
#include <cuco/static_map.cuh>
#include <cstdint>

struct GPURunStats {
  int64_t num_processed = 0;
  int64_t num_pseudoaligned = 0;
  int64_t num_unique = 0;
};

void gpu_run(ProgramOptions& opt, const KmerIndex& index, GPURunStats* out_stats = nullptr);
void gpu_run(ProgramOptions& opt, const GPUIndex& index, GPURunStats* out_stats = nullptr);


#endif // GPU_PROCESS_READS_CUH
