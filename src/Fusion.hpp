#include <sstream>
/** -- fusion functions -- **/

void printTranscripts(const KmerIndex& index, std::stringstream& o, const std::vector<int> u) {
  for (int i = 0; i < u.size(); i++) {
    if (i > 0) {
      o << ",";
    }
    o << index.target_names_[u[i]];
  }
}

std::vector<int> simpleIntersect(const KmerIndex& index, const std::vector<std::pair<KmerEntry,int>>& v) {
  if (v.empty()) {
    return {};
  }
  int ec = index.dbGraph.ecs[v[0].first.contig];
  int lastEC = ec;
  std::vector<int> u = index.ecmap[ec];

  for (int i = 1; i < v.size(); i++) {
    if (v[i].first.contig != v[i-1].first.contig) {
      ec = index.dbGraph.ecs[v[i].first.contig];
      if (ec != lastEC) {
        u = index.intersect(ec, u);
        lastEC = ec;
        if (u.empty()) {
          return u;
        }
      }
    }
  }
  return u;

}

void searchFusion(const KmerIndex &index, const ProgramOptions& opt,
  const MinCollector& tc, MasterProcessor& mp, int ec,
  const std::string &s1, std::vector<std::pair<KmerEntry,int>> &v1,
  const std::string &s2, std::vector<std::pair<KmerEntry,int>> &v2, bool paired) {

  bool partialMap = false;
  if (ec != -1) {
    partialMap = true;
  }
  std::stringstream o;
  // no mapping information
  if (v1.empty() && v2.empty()) {
    o << "UNMAPPED\t\t\t"<< s1 << "\t";
    if (paired) {
      o << s2 << "\t";
    } else {
      o << "\t";
    }
    mp.outputFusion(o);
    return;
  }

  // discordant pairs
  auto u1 = simpleIntersect(index, v1);
  auto u2 = simpleIntersect(index, v2);
  if (!u1.empty() && !u2.empty()) {
    // each pair maps to either end
    return;
  }

  if ((v1.empty() && !u2.empty()) ||  (v2.empty() && !u1.empty())) {
    // read was trimmed
    partialMap = true;
    return;
  }

  // ok so ec == -1 and not both v1 and v2 are empty
  // exactly one of u1 and u2 are empty
  std::vector<std::pair<KmerEntry,int>> vsafe, vsplit;
  // sort v1 and v2 by read position
  auto vsorter =  [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b) {
    return a.second < b.second;
  };

  std::sort(v1.begin(), v1.end(), vsorter);
  std::sort(v2.begin(), v2.end(), vsorter);

  if (!v1.empty() && !v2.empty()) {
    if (u1.empty()) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsafe));
    } else {
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsplit));
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsafe));
    }
  } else if (v1.empty()) {
    if (u2.empty()) {
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsplit));
    } else {
      assert(false);
      return; // can this happen?
    }
  } else if (v2.empty()) {
    if (u1.empty()) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
    } else {
      assert(false);
      return;
    }
  }


  // now we look for a split j s.t. vsplit[0:j] and vsplit[j:] + vsafe both have nonempty ec
  int j = vsplit.size()-1;
  while (!vsplit.empty()) {
    auto x = vsplit.back();
    vsplit.pop_back();
    vsafe.push_back(x);

    auto ut1 = simpleIntersect(index,vsplit);
    auto ut2 = simpleIntersect(index,vsafe);

    if (!ut1.empty() && !ut2.empty()) {
      o << "SPLIT\t";
      printTranscripts(index, o, ut1); o << "\t"; printTranscripts(index, o, ut2);
      o << "\t"  << s1 << "\t";
      // what to put as info?
      if (paired) {
        o << s2 << "\t";
      } else {
        o << "\t";
      }
      o << "splitat=";
      if (u1.empty()) {
        o << "(0,";
      } else {
        o << "(1,";
      }
      o << vsafe.back().second << ")";
      mp.outputFusion(o);
      return;
    } else {
      j--;
    }
  }

  return;

}
