#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "weights.h"

#include <iostream>
#include <vector>

template <typename Index>
struct EMAlgorithm {

    // counts is vector from collector, with indices corresponding to ec ids
    // TODO: initialize alpha a bit more intelligently
    // TODO: refactor to remove dependence on Index
	EMAlgorithm(const ProgramOptions& opt, const Index& idx,
	        const std::vector<int>& counts,
	        const std::vector<double>& eff_lens) :
	    idx_(idx),
	    counts_(counts),
	    num_trans_(idx.num_trans),
	    eff_lens_(eff_lens),
	    alpha_(idx.num_trans, 1/idx.num_trans) // uniform distribution over transcripts
	{}

	void run(const WeightMap& weight_map, size_t n_iter = 500) {
        std::vector<double> next_alpha(alpha_.size(), 0.0);

        assert(weight_map.size() == counts_.size());

        // there is a denominator normalizer for every single equivalence
        // class
        std::vector<double> denom(counts_.size(), 0.0);

	    for (auto i = 0; i < n_iter; ++i) {

            for (auto& ec_kv : idx_.ecmap ) {

                // first, compute the denominator: a normalizer
                // iterate over transcripts in EC map
                auto w_search = weight_map.find(ec_kv.first);

                // everything in ecmap should be in weight_map
                assert( w_search != weight_map.end() );
                assert( w_search->second.size() == ec_kv.second.size() );

                for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
                    denom[ec_kv.first] += ec_kv.second[t_it] *
                        w_search->second[t_it];
                }

                // next, compute the update step

                for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
                    // should we exp(log(.) this?
                    next_alpha[t_it] += counts_[ec_kv.first] *
                        w_search->second[t_it] * alpha_[t_it] / denom[ec_kv.first];
                }
            }

            // reassign alpha_ to next_alpha
            alpha_.swap(next_alpha);

            // clear all next_alpha values 0 for next iteration
            for (auto& a : next_alpha) {
                a = 0.0;
            }

            // clear all denominators
            for (auto& d : denom) {
                d = 0.0;
            }
	    }
	}

	int num_trans_;
	const Index &idx_;
	const std::vector<int> &counts_;
	const std::vector<double>& eff_lens_;
    std::vector<double> alpha_;
};

#endif // KALLISTO_EMALGORITHM_H
