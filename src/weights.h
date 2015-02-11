#ifndef KALLISTO_WEIGHTS_H
#define KALLISTO_WEIGHTS_H

#include "KmerIndex.h"

#include <unordered_map>
#include <utility>
#include <vector>

using WeightMap = std::unordered_map<int, std::vector<double>>;

std::vector<double> calc_eff_lens(const std::vector<int>& lengths, double mean);


WeightMap calc_weights(
  const std::vector<int>& counts,
  const EcMap& ecmap,
  const std::vector<double>& eff_lens);

#endif // KALLISTO_WEIGHTS_H
