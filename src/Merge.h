#ifndef KALLISTO_MERGE_H
#define KALLISTO_MERGE_H

#include "common.h"

bool MergeBatchDirectories(const ProgramOptions &opt, int& num_targets, int64_t& num_processed, 
                            int64_t& num_pseudoaligned, int64_t&  num_unique, int& index_version);


#endif // KALLISTO_MERGE_H