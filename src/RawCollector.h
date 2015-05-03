#ifndef KALLISTO_RAWCOLLECTOR_H
#define KALLISTO_RAWCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "KmerIndex.h"


struct RawCollector {

  RawCollector(KmerIndex& ind, const ProgramOptions& opt) : index(ind), counts(index.ecmap.size(), 0), flens(1000), min_range(opt.min_range), k(opt.k) {}

  int collect(std::vector<std::pair<int,int>>& v1,
              std::vector<std::pair<int,int>>& v2,
              bool nonpaired=false);

  int collect(std::vector<std::pair<int,int>>& v1) {
    std::vector<std::pair<int,int>> dummy;
    return collect(v1,dummy,true);

  }
  int increaseCount(const std::vector<int>& u);
  int decreaseCount(const int ec);

  void write(std::ostream& o) {
    for (int id = 0; id < counts.size(); id++) {
      o << id << "\t" << counts[id] << "\n";
    }
  }
  void loadCounts(ProgramOptions& opt);

  KmerIndex& index;
  std::vector<int> counts;
  std::vector<int> flens;
  int min_range;
  int k;

};


// compute the mean fragment length from a raw_collector
double get_mean_frag_len(const RawCollector& mc);

#endif // KALLISTO_RAWCOLLECTOR_H
