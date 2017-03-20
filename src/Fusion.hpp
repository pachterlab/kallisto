#include <sstream>
#include <set>

/** -- fusion functions -- **/

/**
FUSION output is written to OUT_DIR/fusion.txt where 'OUT_DIR' is the specified
output directory for kallisto. fusion.txt is a tab separated text file with 10 columns
TYPE  NAME1  SEQ1  KPOS1  NAME2  SEQ2  KPOS2  INFO  POS1  POS2

TYPE can be one of SPLIT or PAIR, NAME1 and NAME2 are the sequence identifiers. 
SEQ1 and SEQ2 are the raw read sequences identifie as overlapping the potential fusion.
KPOS1 and KPOS2 are the positions within the reads where the k-mers match.

PAIRs are readpairs where fragment doesn't map, but each read
maps independently, potentially the fusion breakpoint is within the fragment, but not sequenced.

SPLITs are reads where the fusion breakpoint appears within one of the reads

INFO is NA for PAIR but of the form splitat=(a,b)
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

  // find first mapping k-mer
  if (!v.empty()) {
    p = findFirstMappingKmer(v,val);
    km = Kmer((s.c_str()+p));
  }
  

  for (int i = 0; i < u.size(); i++) {
    int tr = u[i];
    if (i > 0) {
      o << ";";
    }
    std::pair<int, bool> xp = index.findPosition(tr, km, val, p);
    o << "(" << index.target_names_[tr] << "," << xp.first << ",";
    if (xp.second) {
      o << "FW)";
    } else {
      o << "RC)";
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


bool checkMapability(const KmerIndex& index, const std::string &s, const std::vector<std::pair<KmerEntry,int>>& v, std::vector<int> &u) {
  const int maxMismatch = 2;
  const int maxSoftclip = 5;
    
  Kmer km;
  KmerEntry val;
  int p;

  if (!v.empty()) {
    p = findFirstMappingKmer(v,val);
    km = Kmer(s.c_str()+p);
  } else {
    return false;
  }
  
  std::vector<int> vtmp; vtmp.reserve(u.size());
  
  for (auto tr : u) {
    auto trpos = index.findPosition(tr, km, val, p);
    int tpos = (int)trpos.first;
    int sz = (int)s.size();
    bool add = true; 
    if (trpos.second) {
      if (tpos < 1 || tpos + sz - 1 > index.target_seqs_[tr].size()) {
        add = false;
      } else {
        //std::cout << index.target_seqs_[tr].substr(tpos,sz) << std::endl;
        //std::cout << s << std::endl;
        int mis = 0;
        for (int i = 0; i < sz - maxSoftclip; i++) {
          if (index.target_seqs_[tr][tpos-1 + i] != s[i]) {
            ++mis;
            if (mis > maxMismatch) {
              break;
            }
          }
        }
        add = (mis <= maxMismatch);
      }
    }  else {
      if (tpos > index.target_seqs_[tr].size() || tpos - sz < 1) {
        add = false;
      } else {      
        std::string rs = revcomp(s);
        //std::cout << index.target_seqs_[tr].substr(tpos - sz, sz) << std::endl;
        //std::cout << rs << std::endl;
        int mis = 0;
        for (int i = sz-1; i >= maxSoftclip; i--) {
          if (index.target_seqs_[tr][tpos-sz+i] != rs[sz]) {
            ++mis;
            if (mis > maxMismatch) {
              break;
            }
          }
        }
        add = (mis <= maxMismatch);
      }
    }
    
    if (add) {
      vtmp.push_back(tr);
    }
    
  }
 
  
  if (vtmp.empty()) {
    return false;
  }
  
  if (vtmp.size() < u.size()) {
    u = vtmp; // copy
  }
  
  return true;
  
}

// returns true if the intersection of the union of EC classes for v1 and v2 respectively is empty
bool checkUnionIntersection(const KmerIndex& index, const std::string &s1, const std::string &s2, std::pair<int,int> &p1, std::pair<int,int> &p2) { 
//const std::vector<std::pair<KmerEntry,int>>& v1, const std::vector<std::pair<KmerEntry,int>>& v2) {
  
  std::set<int> su1,su2;
  
  auto union_set = [&](const std::string &s, std::pair<int,int> &p) {
    p = {-1,-1};
    std::set<int> su;
    
    KmerIterator kit(s.c_str()), kit_end;
    int lastEC = -1;
    for (int i = 0; kit != kit_end; ++i,++kit) {
      auto search = index.kmap.find(kit->first.rep());
      if (search != index.kmap.end()) {
        if (p.first == -1) {
          p.first = kit->second;
          p.second = p.first +1;
        } else {
          if (p.second + 1 == kit->second) {
            p.second++;
          }
        }
        const KmerEntry &val = search->second;
        int ec = index.dbGraph.ecs[val.contig];
        if (ec != -1 && ec != lastEC) {
          su.insert(index.ecmap[ec].begin(), index.ecmap[ec].end());
          lastEC = ec;
        }
      }
    }
    
    /*
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
    */
    return su;
  };

  su1 = union_set(s1,p1);
  su2 = union_set(s2,p2);
  
  if (su1.empty() || su2.empty()) {
    return false; // TODO, decide on this
  }
  
  if (su2.size() < su1.size()) {
    swap(su1,su2);
  }
  for (auto x : su1) {
    if (su2.count(x) != 0) {
      return false;
    }
  }
  return true;
}


void searchFusion(const KmerIndex &index, const ProgramOptions& opt,
  const MinCollector& tc, MasterProcessor& mp, int ec,
  const std::string &n1, const std::string &s1, std::vector<std::pair<KmerEntry,int>> &v1,
  const std::string &n2, const std::string &s2, std::vector<std::pair<KmerEntry,int>> &v2, bool paired) {

  bool partialMap = false;
  if (ec != -1) {
    partialMap = true;
  }
  std::stringstream o;
  
  // no mapping information
  if (v1.empty() && v2.empty()) {
    return; // consider splitting in case either is empty
  }


  auto u1 = simpleIntersect(index, v1);
  auto u2 = simpleIntersect(index, v2);

  // discordant pairs
  if (!v1.empty() && !v2.empty()) {
    if (!u1.empty() && !u2.empty()) {
      std::pair<int,int> p1,p2;
      if (checkUnionIntersection(index, s1, s2, p1, p2)) {
        // each pair maps to either end
        // each pair maps to either end
        o << "PAIR\t";
        o << n1 << "\t";
        o << s1 << "\t";
        o << p1.first << "," << p1.second << "\t"; 
        if (paired) {
          o << n2 << "\t";
          o << s2 << "\t";
          o << p2.first << "," << p2.second << "\t\tNA\t";
        } else {
          o << "\t\t\tNA\t";
        }
        printTranscripts(index, o, s1, v1, u1); o << "\t"; printTranscripts(index, o, s2, v2, u2);
        mp.outputFusion(o);
        return;
      }
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
      std::pair<int,int> p1,p2;
      std::string st1,st2;
      if (u1.empty()) {
        st1 = s1.substr(0,vsafe.back().second);
        st2 = s2;
      } else {
        st1 = s1;
        st2 = s2.substr(0,vsafe.back().second);
      }
      if (checkUnionIntersection(index, st1, st2, p1, p2)) { // need to check this more carefully
        o << "SPLIT\t";
        o << n1 << "\t" << s1 << "\t" << p1.first << "," << p1.second << "\t";
        // what to put as info?
        if (paired) {
          o << n2 << "\t" << s2 << "\t" << p2.first << "," << p2.second << "\t";
        } else {
          o << "\t\t\t";
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
