#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include <iostream>
#include <vector>

template <typename Index>
struct EMAlgorithm {

    // counts is vector from collector, with indices corresponding to ec ids
	EMAlgorithm(const ProgramOptions& opt, const Index& idx,
	        const std::vector<int>& counts,
	        const std::vector<double>& eff_lens) :
	    idx_(idx),
	    counts_(counts),
	    num_trans_(idx.num_trans),
	    eff_lens_(eff_lens),
	    alpha_(idx.num_trans, 1/idx.num_trans) // uniform distribution over transcripts
	{}

	void run(size_t n_iter = 500) {

	}

	int num_trans_;
	const Index &idx_;
	const std::vector<int> &counts_;
	const std::vector<double>& eff_lens_;
    std::vector<double> alpha_;
};

#endif // KALLISTO_EMALGORITHM_H
