#include <sstream>
#include <set>
/** -- fusion functions -- **/

/**
FUSION output is written to OUT_DIR/fusion.txt where 'OUT_DIR' is the specified
output directory for kallisto. fusion.txt is a tab separated text file with 6 columns
TYPE	SEQ1	SEQ2	INFO	POS1	POS2

TYPE can be one of SPLIT or PAIR, SEQ1 and SEQ2 are the raw read sequences identified
as overlapping the potential fusion.

PAIRs are readpairs where fragment doesn't map, but each read
maps independently, potentially the fusion breakpoint is within the fragment, but not sequenced.

SPLITs are reads where the fusion breakpoint appears within one of the reads

INFO is empty for PAIR but of the form splitat=(a,b)
where a=0,1 indicating whether the first or the second pair is split and b is the 0-based position of the
split within the read.

POS1 and POS2 are comma separated lists of locations of SEQ1 and SEQ2 respectively. Each position is of
the form (tx_name,pos,strand) where pos is 0-based from the start of the transcript tx_name and strand
is either FW or RE depending on whether the read aligns to the forward or reverse of the transcript.
**/

void printTranscripts(const KmerIndex& index, std::stringstream& o, const std::string s,
  const std::vector<std::pair<KmerEntry,int>>& v, const std::vector<int> u) {

  Kmer km;
  KmerEntry val;
  int p;

  if (!v.empty()) {
    val = v[0].first;
    p = v[0].second;
    for (auto &x : v) {
      if (x.second < p) {
        val = x.first;
        p = x.second;
      }
    }
    km = Kmer(s.c_str()+p);
  }


  for (int i = 0; i < u.size(); i++) {
    int tr = u[i];
    if (i > 0) {
      o << ",";
    }
    std::pair<int, bool> xp = index.findPosition(tr, km, val, p);
    o << "(" << index.target_names_[tr] << "," << xp.first << ",";
    if (xp.second) {
      o << "FW)";
    } else {
      o << "RE)";
    }
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

// returns true if the intersection of the union of EC classes for v1 and v2 respectively is empty
bool checkUnionIntersection(const KmerIndex& index, const std::vector<std::pair<KmerEntry,int>>& v1, const std::vector<std::pair<KmerEntry,int>>& v2) {
  if (v1.empty() || v2.empty()) {
    return true;
  }
  
  std::set<int> s1,s2;
  
  auto union_set = [&](const std::vector<std::pair<KmerEntry,int>>& v) {
    std::set<int> s;
    int ec = index.dbGraph.ecs[v[0].first.contig];
    int lastEC = ec;
    //const std::vector<int> &u = index.ecmap[ec];
    s.insert(index.ecmap[ec].begin(), index.ecmap[ec].end());
    for (int i = 1; i < v.size(); i++) {
      if (v[i].first.contig != v[i-1].first.contig) {
        ec = index.dbGraph.ecs[v[i].first.contig];
        if (ec != lastEC) {
          s.insert(index.ecmap[ec].begin(), index.ecmap[ec].end());
          lastEC = ec;
        }
      }
    }
    return s;
  };

  s1 = union_set(v1);
  s2 = union_set(v2);
  if (s2.size() < s1.size()) {
    swap(s1,s2);
  }
  for (auto x : s1) {
    if (s2.count(x) != 0) {
      return false;
    }
  }
  return true;
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
  if (v1.empty() || v2.empty()) {
    return; // consider splitting in case either is empty
  }

  // discordant pairs
  auto u1 = simpleIntersect(index, v1);
  auto u2 = simpleIntersect(index, v2);
  if (!u1.empty() && !u2.empty()) {
    if (checkUnionIntersection(index,v1,v2)) {
      // each pair maps to either end
      // each pair maps to either end
      o << "PAIR\t";
      o << s1 << "\t";
      if (paired) {
      o << s2 << "\t\t";
      } else {
        o << "\t\t";
      }
      printTranscripts(index, o, s1, v1, u1); o << "\t"; printTranscripts(index, o, s2, v2, u2);
      mp.outputFusion(o);
      return;
    }
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
  bool split1 = true;
  if (!v1.empty() && !v2.empty()) {
    if (u1.empty()) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsafe));
    } else {
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsplit));
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsafe));
      split1 = false;
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
      if (checkUnionIntersection(index, vsplit, vsafe)) {
        o << "SPLIT\t";
        o << s1 << "\t";
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
        o << vsafe.back().second << ")\t";
        // fix this
        if (split1) {
          printTranscripts(index, o, s1, vsplit, ut1); o << "\t"; printTranscripts(index, o, s2, vsafe, ut2);
        } else {
          printTranscripts(index, o, s1, vsafe, ut2); o << "\t"; printTranscripts(index, o, s2, vsplit, ut1);
        }
        mp.outputFusion(o);
        return;
      }
    } else {
      j--;
    }
  }

  return;

}
