#ifndef KALLISTO_WEIGHTS_H
#define KALLISTO_WEIGHTS_H

#include "KmerIndex.h"
#include "MinCollector.h"
#include <unordered_map>
#include <utility>
#include <vector>

using WeightMap = std::vector<std::vector<double>>;

// this function takes the 'mean_fl_trunc' from MinCollector and simply gives
// you back a 'mean fragment length' for every single transcript. this avoids
// you having to check the length every single time
std::vector<double> get_frag_len_means(const std::vector<int>& lengths,
    const std::vector<double>& mean_frag_len_trunc);

std::vector<double> calc_eff_lens(const std::vector<int>& lengths, double mean);
std::vector<double> calc_eff_lens(const std::vector<int>& lengths,
    const std::vector<double>& means);
std::vector<double> update_eff_lens(double mean, const MinCollector& tc, const KmerIndex &index, const std::vector<double> alpha, const std::vector<double> eff_lens);


WeightMap calc_weights(
  const std::vector<int>& counts,
  const EcMap& ecmap,
  const std::vector<double>& eff_lens);

#endif // KALLISTO_WEIGHTS_H
