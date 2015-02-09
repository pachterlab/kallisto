#ifndef PY_KALLISTO_H
#define PY_KALLISTO_H

#include "KmerIndex.h"

// counts gets overwritten with the counts per ecid
KmerIndex read_kal_py_index(const std::string& fname,
        const std::string& fasta,
        std::string& py_counts,
        std::vector<int>& counts);

#endif
