#ifndef KALLISTO_WEIGHTS_H
#define KALLISTO_WEIGHTS_H

#include "KmerIndex.h"
#include "MinCollector.h"
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

struct MinCollector;

using WeightMap = std::vector<std::vector<double>>;

// this function takes the 'mean_fl_trunc' from MinCollector and simply gives
// you back a 'mean fragment length' for every single transcript. this avoids
// you having to check the length every single time
std::vector<double> get_frag_len_means(const std::vector<int>& lengths,
    const std::vector<double>& mean_frag_len_trunc);

// XXX: DEPRECATED. See overloaded function below.
std::vector<double> calc_eff_lens(const std::vector<int>& lengths, double mean);

// @param lengths the lengths of all the targets
// @param means the mean frag len of every transcript
std::vector<double> calc_eff_lens(const std::vector<int>& lengths,
    const std::vector<double>& means);

std::vector<double> update_eff_lens(const std::vector<double>& means,
    const MinCollector& tc,
    const KmerIndex &index, const std::vector<double>& alpha,
    const std::vector<double>& eff_lens, std::vector<double>& post_bias,
    const ProgramOptions& opt );


WeightMap calc_weights(
  const std::vector<int>& counts,
  const EcMap& ecmap,
  const std::vector<double>& eff_lens);


// truncated gaussian fragment length distribution
//
// use this for single-end reads since we don't have an empirical distribution
//
// start: inclusive
// stop: exclusive
// let's pretend all the input is sane
std::vector<double> trunc_gaussian_fld(int start, int stop, double mean,
    double sd);

std::vector<int> trunc_gaussian_counts(int start, int stop, double mean,
        double sd, int total_count);


#endif // KALLISTO_WEIGHTS_H
