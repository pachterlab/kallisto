#ifndef GPU_BUS_PROCESS_CUH
#define GPU_BUS_PROCESS_CUH

#include "GPUIndexFormat.h"
#include "GPUProcessReads.cuh"  // GPURunStats
#include "KmerIndex.h"
#include "common.h"

void gpu_bus_run(ProgramOptions& opt, const KmerIndex& index,
                 const std::string& start_time, int argc, char** argv,
                 GPURunStats* out_stats = nullptr);
void gpu_bus_run(ProgramOptions& opt, const GPUIndex& index,
                 const std::string& start_time, int argc, char** argv,
                 GPURunStats* out_stats = nullptr);

#endif  // GPU_BUS_PROCESS_CUH
