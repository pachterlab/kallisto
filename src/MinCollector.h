#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "KmerIndex.h"

const int MAX_FRAG_LEN = 1000;

struct MinCollector {

  MinCollector(KmerIndex& ind, const ProgramOptions& opt)
    :
      index(ind),
      counts(index.ecmap.size(), 0),
      flens(MAX_FRAG_LEN),
      bias3(4096),
      bias5(4096),
      min_range(opt.min_range),
      k(opt.k),
      mean_fl(0.0),
      has_mean_fl(false),
      mean_fl_trunc(MAX_FRAG_LEN, 0.0),
      has_mean_fl_trunc(false)
      // eff_len_cache(MAX_FRAG_LEN, 0.0)
       {}

  int collect(std::vector<std::pair<KmerEntry,int>>& v1,
              std::vector<std::pair<KmerEntry,int>>& v2,
              bool nonpaired=false);

  int collect(std::vector<std::pair<KmerEntry,int>>& v1) {
    std::vector<std::pair<KmerEntry,int>> dummy;
    return collect(v1,dummy,true);

  }
  int increaseCount(const std::vector<int>& u);
  int decreaseCount(const int ec);

  std::vector<int> intersectECs(std::vector<std::pair<KmerEntry,int>>& v) const;


  void write(std::ostream& o) {
    for (int id = 0; id < counts.size(); id++) {
      o << id << "\t" << counts[id] << "\n";
    }
  }
  void loadCounts(ProgramOptions& opt);

  bool countBias(const char *s1, const char *s2, const std::vector<std::pair<KmerEntry,int>> v1, const std::vector<std::pair<KmerEntry,int>> v2, bool paired);

  double get_mean_frag_len() const;

  void get_mean_frag_lens_trunc() const;

  // this function is a bit of a hack. should only be used with SE data
  void init_mean(double mean);

  /* double eff_len(int len) const; */

  KmerIndex& index;
  std::vector<int> counts;
  std::vector<int> flens;
  std::vector<int> bias3, bias5;
  int min_range;
  int k;

  double mean_fl;
  bool has_mean_fl;

  std::vector<double> mean_fl_trunc;
  bool has_mean_fl_trunc;

  // std::vector<double> eff_len_cache;
};

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y);

// compute the mean fragment length from a min_collector


int hexamerToInt(const char *s, bool revcomp);

#endif // KALLISTO_MINCOLLECTOR_H
