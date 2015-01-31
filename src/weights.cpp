#include "weights.h"

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
