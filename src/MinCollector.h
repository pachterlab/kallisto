#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "KmerIndex.h"


struct MinCollector {

  MinCollector(KmerIndex& ind, const ProgramOptions& opt) : index(ind), counts(index.ecmap.size(), 0), flens(1000), bias3(4096), bias5(4096), min_range(opt.min_range), k(opt.k), mean_fl(0.0), has_mean_fl(false) {}

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
  
  
  KmerIndex& index;
  std::vector<int> counts;
  std::vector<int> flens;
  std::vector<int> bias3, bias5;
  int min_range;
  int k;
  double mean_fl;
  bool has_mean_fl;
};

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y);

// compute the mean fragment length from a min_collector


int hexamerToInt(const char *s, bool revcomp);

#endif // KALLISTO_MINCOLLECTOR_H
