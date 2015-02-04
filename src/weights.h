#ifndef KALLISTO_WEIGHTS_H
#define KALLISTO_WEIGHTS_H

#include "KmerIndex.h"

#include <unordered_map>
#include <utility>
#include <vector>

using WeightMap = std::unordered_map<int, std::vector<double>>;

std::vector<double> calc_eff_lens(const std::vector<int>& lengths, double mean)
{
    // for now do the total naive thing and subtract mean frag length
    std::vector<double> eff_lens;
    eff_lens.reserve(lengths.size());

    for (auto& cur_len: lengths)
    {
        eff_lens.push_back( static_cast<double>(cur_len) - mean + 1.0 );
    }

    return eff_lens;
}


WeightMap calc_weights(
        const std::vector<int>& counts,
        const EcMap& ecmap,
        const std::vector<double>& eff_lens)
{

    // TODO: throw some assertions in here to make sure the length of counts
    // and ec map are correct... as well as eff_lens size is reasonable

    // weights are stored _exactly_ in the same orientation as the ec map
    WeightMap weights;

    for (auto& kv : ecmap) {

			//std::cout << kv.first;
        std::vector<double> trans_weights;
        trans_weights.reserve(kv.second.size());

        for (auto& trans_id : kv.second) {
            trans_weights.push_back( static_cast<double>(counts[kv.first]) /
                    eff_lens[trans_id] );
        }

        weights.insert( {kv.first, trans_weights} );
    }


    return weights;
}

#endif // KALLISTO_WEIGHTS_H
