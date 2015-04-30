#include "MinCollector.h"
#include <algorithm>

// utility functions

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y) {
  std::vector<int> v;
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      v.push_back(*a);
      ++a;
      ++b;
    }
  }
  return v;
}


int MinCollector::collect(std::vector<std::pair<int,int>>& v1,
                          std::vector<std::pair<int,int>>& v2, bool nonpaired) {

  if (v1.empty()) {
    return -1;
  } else if (!nonpaired && v2.empty()) {
    return -1;
  }

  std::vector<int> u1 = intersectECs(v1);
  std::vector<int> u2 = intersectECs(v2);

  std::vector<int> u;

  if (u1.empty() && u2.empty()) {
    return -1;
  }

  // non-strict intersection.
  if (u1.empty()) {
    u = u2;
  } else if (u2.empty()) {
    u = u1;
  } else {
    u = intersect(u1,u2);
  }

  if (u.empty()) {
    return -1;
  }
  return increaseCount(u);
}

int MinCollector::increaseCount(const std::vector<int>& u) {
  if (u.empty()) {
    return -1;
  }
  if (u.size() == 1) {
    int ec = u[0];
    ++counts[ec];
    return ec;
  }
  auto search = index.ecmapinv.find(u);
  if (search != index.ecmapinv.end()) {
    // ec class already exists, update count
    ++counts[search->second];
    return search->second;
  } else {
    // new ec class, update the index and count
    auto necs = counts.size();
    //index.ecmap.insert({necs,u});
    index.ecmap.push_back(u);
    index.ecmapinv.insert({u,necs});
    counts.push_back(1);
    return necs;
  }
}

int MinCollector::decreaseCount(const int ec) {
  assert(ec >= 0 && ec <= index.ecmap.size());
  --counts[ec];
  return ec;
}

struct ComparePairsBySecond {
  bool operator()(std::pair<int,int> a, std::pair<int,int> b) {
    return a.second < b.second;
  }
};

std::vector<int> MinCollector::intersectECs(std::vector<std::pair<int,int>>& v) const {
  if (v.empty()) {
    return {};
  }
  sort(v.begin(), v.end()); // sort by ec, and then first position

  /*
  std::vector<std::pair<int,int>> vp;
  vp.push_back(v[0]);
  for (int i = 1; i < v.size(); i++) {
    if (v[i].first != v[i-1].first) {
      vp.push_back(v[i-1]);
    }
    }

  sort(vp.begin(), vp.end(), ComparePairsBySecond{});
  */

  int count = 1; // how many k-mer support the ec
  std::vector<int> u = index.ecmap[v[0].first];

  for (int i = 1; i < v.size(); i++) {
    if (v[i].first != v[i-1].first) {
      u = index.intersect(v[i].first, u);
      if (u.empty()) {
        return u;
      }
    }
  }

  /*for (auto &x : vp) {
    //tmp = index.intersect(x.first,u);
    u = index.intersect(x.first,u);
    //if (!tmp.empty()) {
     // u = tmp;
      //count++; // increase the count
     // }
  }*/

  // if u is empty do nothing
  /*if (u.empty()) {
    return u;
    }*/

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return u;
}


void MinCollector::loadCounts(ProgramOptions& opt) {
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

double get_mean_frag_len(const MinCollector& mc) {
  auto total_counts = 0;
  double total_mass = 0.0;

  for ( size_t i = 0 ; i < mc.flens.size(); ++i ) {
    total_counts += mc.flens[i];
    total_mass += static_cast<double>(mc.flens[i] * i);
  }

  return total_mass / static_cast<double>(total_counts);
}
