#include "RawCollector.h"
#include <algorithm>

// utility functions



int RawCollector::collect(std::vector<std::pair<int,int>>& v1,
                          std::vector<std::pair<int,int>>& v2, bool nonpaired) {

  for (auto x : v1) {
    ++counts[x.first];
  }

  if (!nonpaired) {
    for (auto x : v2) {
      ++counts[x.first];
    }
  }

  return -1;
}

int RawCollector::increaseCount(const std::vector<int>& u) {
  // not implemented
  return -1;
}

int RawCollector::decreaseCount(const int ec) {
  assert(ec >= 0 && ec <= index.ecmap.size());
  --counts[ec];
  return ec;
}

struct ComparePairsBySecond {
  bool operator()(std::pair<int,int> a, std::pair<int,int> b) {
    return a.second < b.second;
  }
};

void RawCollector::loadCounts(ProgramOptions& opt) {
  int num_ecs = counts.size();
  counts.clear();
  std::ifstream in((opt.output + "/counts.txt"));
  int i = 0;
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int j,c;
      ss >> j;
      ss >> c;
      if (j != i) {
        std::cerr << "Error: equivalence class does not match index. Found "
                  << j << ", expected " << i << std::endl;
        exit(1);
      }
      counts.push_back(c);
      i++;
    }

    if (i != num_ecs) {
      std::cerr << "Error: number of equivalence classes does not match index. Found "
                << i << ", expected " << num_ecs << std::endl;
      exit(1);
    }
  } else {
    std::cerr << "Error: could not open file " << opt.output << "/counts.txt" << std::endl;
    exit(1);

  }

}

double get_mean_frag_len(const RawCollector& mc) {
  auto total_counts = 0;
  double total_mass = 0.0;

  for ( size_t i = 0 ; i < mc.flens.size(); ++i ) {
    total_counts += mc.flens[i];
    total_mass += static_cast<double>(mc.flens[i] * i);
  }

  return total_mass / static_cast<double>(total_counts);
}
