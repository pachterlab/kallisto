#pragma once

#include "common.h"
#include "KmerIndex.h"
#include "GPUIndexFormat.h"

#include <cstdint>
#include <vector>

struct MetalRunStats {
    int64_t num_processed    = 0;
    int64_t num_pseudoaligned = 0;
    int64_t num_unique        = 0;
};

void metal_run(ProgramOptions& opt, const KmerIndex& index,
               MetalRunStats* out_stats = nullptr);

void metal_run(ProgramOptions& opt, const GPUIndex& index,
               MetalRunStats* out_stats = nullptr);
