/*
#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "MinCollector.h"


#include "common.h"
*/

#include <fstream>

#include "ProcessReads.h"
#include "kseq.h"
#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

void outputPseudoBam(const KmerIndex &index, int ec, const kseq_t *seq1, const std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, const std::vector<std::pair<KmerEntry,int>> &v2, bool paired);
void searchFusion(const KmerIndex &index, const ProgramOptions& opt, const MinCollector& tc, std::ofstream& o, int ec,const kseq_t *seq1,  std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2,  std::vector<std::pair<KmerEntry,int>> &v2, bool paired);
void revseq(char *b1, char *b2, const kseq_t *seq);
void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength);

void printVector(const std::vector<int>& v, std::ostream& o) {
  o << "[";
  int i = 0;
  for (auto x : v) {
    if (i>0) {
      o << ", ";
    }
    o << x;
    i++;
  }
  o << "]";
}

bool isSubset(const std::vector<int>& x, const std::vector<int>& y) {
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      return false;
    } else if (*b < *a) {
      ++b;
    } else {
      ++a;
      ++b;
    }
  }
  return (a==x.end());
}


int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc) {
  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  int tlencount = (opt.fld == 0.0) ? 10000 : 0;
  size_t numreads = 0;
  size_t nummapped = 0;
  size_t testFusion = 0;
  bool paired = !opt.single_end;

  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  std::vector<std::pair<KmerEntry,int>> v1, v2;
  v1.reserve(1000);
  v2.reserve(1000);

  std::ofstream ofusion;
  if (opt.fusion) {
    ofusion.open(opt.output + "/fusion.txt");
    ofusion << "TYPE\tTRANSCRIPTS1\tTRANSCRIPTS2\tNAME1\tSEQ1\tNAME2\tSEQ2\tINFO\n";
  }
  

  int l1 = 0,l2 = 0; // length of read

  if (paired) {
    std::cerr << "[quant] running in paired-end mode" << std::endl;
  } else {
    std::cerr << "[quant] running in single-end mode" << std::endl;
  }

  for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {
    if (paired) {
      std::cerr << "[quant] will process pair " << (i/2 +1) << ": "  << opt.files[i] << std::endl 
                << "                             " << opt.files[i+1] << std::endl;
    } else {
      std::cerr << "[quant] will process file " << i+1 << ": " << opt.files[i] << std::endl;
    }
  }

  // for each file
  std::cerr << "[quant] finding pseudoalignments for the reads ..."; std::cerr.flush();
  int biasCount = 0;
  const int maxBiasCount = 1000000;

  if (opt.pseudobam) {
    index.writePseudoBamHeader(std::cout);
  }
  
  
  for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {
  
    fp1 = gzopen(opt.files[i].c_str(), "r");
    seq1 = kseq_init(fp1);
    if (paired) {
      fp2 = gzopen(opt.files[i+1].c_str(),"r");
      seq2 = kseq_init(fp2);
    }

    

    // for each read
    while (true) {
      l1 = kseq_read(seq1);
      if (paired) {
        l2 = kseq_read(seq2);
      }
      if (l1 <= 0) {
        break;
      }
      if (paired && l2 <= 0) {
        break;
      }

      numreads++;
      v1.clear();
      v2.clear();
      // process read
      index.match(seq1->seq.s, seq1->seq.l, v1);
      if (paired) {
        index.match(seq2->seq.s, seq2->seq.l, v2);
      }

      // collect the target information
      int ec = tc.collect(v1, v2, !paired);

      if (paired && ec != -1 && (v1.empty() || v2.empty()) && tlencount == 0) {
        // inspect the positions
        int fl = (int) tc.get_mean_frag_len();
        if (v2.empty()) {
          int p = -1;
          KmerEntry val;
          Kmer km;
          
          if (!v1.empty()) {
            val = v1[0].first;
            p = v1[0].second;
            for (auto &x : v1) {
              if (x.second < p) {
                val = x.first;
                p = x.second;
              }
            }
            km = Kmer((seq1->seq.s+p));
          }
          if (!v2.empty()) {
            val = v2[0].first;
            p = v2[0].second;
            for (auto &x : v2) {
              if (x.second < p) {
                val = x.first;
                p = x.second;
              }
            }
            km = Kmer((seq2->seq.s+p));
          }

          std::vector<int> u = index.ecmap[ec]; // copy
          std::vector<int> v;
          v.reserve(u.size());

          for (auto tr : u) {
            auto x = index.findPosition(tr, km, val, p);
            if (x.second && x.first + fl <= index.target_lens_[tr]) {
              v.push_back(tr);
            }
            if (!x.second && x.first - fl >= 0) {
              v.push_back(tr);
            }
          }

          if (v.size() < u.size()) {
            // fix the ec
            tc.decreaseCount(ec);
            ec = tc.increaseCount(v);
          }
        }
      }
      
      if (ec != -1) {
        nummapped++;
        // collect sequence specific bias info
        if (opt.bias && biasCount < maxBiasCount) {
          if (tc.countBias(seq1->seq.s, (paired) ? seq2->seq.s : nullptr, v1, v2, paired)) {
            biasCount++;
          }
        }
      }

      // collect fragment length info
      if (tlencount > 0 && paired && 0 <= ec &&  ec < index.num_trans && !v1.empty() && !v2.empty()) {
        // try to map the reads
        int tl = index.mapPair(seq1->seq.s, seq1->seq.l, seq2->seq.s, seq2->seq.l, ec);
        if (0 < tl && tl < tc.flens.size()) {
          tc.flens[tl]++;
          tlencount--;
        }
      }

      if (opt.fusion && (ec == -1 || v1.empty() || v2.empty())) {
        testFusion++;
        searchFusion(index,opt,tc,ofusion,  ec,seq1,v1,seq2,v2,paired);
      }

      // pseudobam
      if (opt.pseudobam) {
        outputPseudoBam(index, ec, seq1, v1, seq2, v2, paired);
      }

      // see if we can do better
      /* if (ec >= index.num_trans) { */
      /*   // non-trivial ec */
      /*   auto maxpair = [](int max_, std::pair<int,int> a) {return (max_ > a.second) ? max_ : a.second;}; */
      /*   int maxPos1 = accumulate(v1.begin(), v1.end(), -1, maxpair); */
      /*   int maxPos2 = -1; */
      /*   if (paired) { */
      /*     maxPos2 = accumulate(v2.begin(), v2.end(), -1, maxpair); */
      /*   } */
      /*   bool better1 = false; */
      /*   bool better2 = false; */
      /*   bool cand = false; */
      /*   if (maxPos1 < seq1->seq.l - index.k) { */
      /*     cand = true; */
      /*     better1 = index.matchEnd(seq1->seq.s, seq1->seq.l, v1, maxPos1); */
      /*   } */
      /*   if (paired) { */
      /*     if (maxPos2 < seq2->seq.l - index.k) { */
      /*       cand = true; */
      /*       better2 = index.matchEnd(seq2->seq.s, seq2->seq.l, v2, maxPos2); */
      /*     } */
      /*   } */

      /*   if (cand) { */
      /*     betterCand++; */
      /*   } */
      /*   if (better1 || better2) { */
      /*     betterCount++; */
      /*     tc.decreaseCount(ec); */
      /*     ec = tc.collect(v1,v2,!paired); */
      /*   } */
      /* } */

      if (opt.verbose && numreads % 100000 == 0 ) {
        std::cerr << "[quant] Processed " << pretty_num(numreads) << std::endl;
      }
    }
    gzclose(fp1);
    if (paired) {
      gzclose(fp2);
    }
  }
    
  kseq_destroy(seq1);
  if (paired) {
    kseq_destroy(seq2);
  }

  if (opt.fusion) {
    ofusion.close();
  }

  std::cerr << " done" << std::endl;

  //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;

  if (opt.bias) {
    std::cerr << "[quant] learning parameters for sequence specific bias" << std::endl;
  }
  
  std::cerr << "[quant] processed " << pretty_num(numreads) << " reads, "
    << pretty_num(nummapped) << " reads pseudoaligned" << std::endl;

  if (opt.fusion) {
    std::cerr << "[quant] classified " << testFusion << " reads for fusion breakpoints" << std::endl;
  }
  
  /*
  for (int i = 0; i < 4096; i++) {
    std::cout << i << " " << tc.bias5[i] << " " << tc.bias3[i] << "\n";
    }*/
  
  // write output to outdir
  if (opt.write_index) {
    std::string outfile = opt.output + "/counts.txt";
    std::ofstream of;
    of.open(outfile.c_str(), std::ios::out);
    tc.write(of);
    of.close();
  }

  return numreads;
}


/** -- fusion functions -- **/


void printTranscripts(const KmerIndex& index, std::ofstream& o, const std::vector<int> u) {
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

void searchFusion(const KmerIndex &index, const ProgramOptions& opt, const MinCollector& tc, std::ofstream& o, int ec,const kseq_t *seq1, std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, std::vector<std::pair<KmerEntry,int>> &v2, bool paired) {

  bool partialMap = false;
  if (ec != -1) {
    partialMap = true;
  }

  // no mapping information
  if (v1.empty() && v2.empty()) {
    o << "UNMAPPED\t\t\t" << seq1->name.s << "\t" << seq1->seq.s << "\t";
    if (paired) {
      o << seq2->name.s << "\t" << seq2->seq.s << "\t\n";
    } else {
      o << "\t\t\n";
    }
    return;
  }

  // discordant pairs
  auto u1 = simpleIntersect(index, v1);
  auto u2 = simpleIntersect(index, v2);
  if (!u1.empty() && !u2.empty()) {
    // each pair maps to either end
    o << "PAIR\t";
    printTranscripts(index, o, u1); o << "\t"; printTranscripts(index, o, u2);
    o << "\t"  << seq1->name.s << "\t" << seq1->seq.s << "\t";
    // what to put as info?
    if (paired) {
      o << seq2->name.s << "\t" << seq2->seq.s << "\t\n";
    } else {
      o << "\t\t\n";
    }
    return;
  }

  if ((v1.empty() && !u2.empty()) ||  (v2.empty() && !u1.empty())) {
    // read was trimmed
    partialMap = true;
  }

  // partial mapping info, e.g. ec != -1 or read was trimmed
  if (partialMap) {
    if (v2.empty()) {
      o << "LEFTMAPPED\t";
      printTranscripts(index, o, index.ecmap[ec]);
      o << "\t\t";
    } else if (v1.empty()) {
      o << "RIGHTMAPPED\t\t";
      printTranscripts(index, o, index.ecmap[ec]);
      o << "\t";
    } else {
      assert(false);
      return;
    }
    o << seq1->name.s << "\t" << seq1->seq.s << "\t";
    if (paired) {
      o << seq2->name.s << "\t" << seq2->seq.s << "\t\n";
    } else {
      o << "\t\t\n";
    }
    return; 
  }
  

  // ok so ec == -1 and not both v1 and v2 are empty
  // exactly one of u1 and u2 are empty
  std::vector<std::pair<KmerEntry,int>> vsafe, vsplit;

  if (!v1.empty() && !v2.empty()) {
    if (u1.empty()) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
      std::copy(v2.rbegin(), v2.rend(), std::back_inserter(vsafe));
    } else {
      std::copy(v2.rbegin(), v2.rend(), std::back_inserter(vsplit));
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsafe));
    }
  } else if (v1.empty()) {
    if (u2.empty()) {
      std::copy(v2.rbegin(), v2.rend(), std::back_inserter(vsplit));
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
      o << "\t"  << seq1->name.s << "\t" << seq1->seq.s << "\t";
      // what to put as info?
      if (paired) {
        o << seq2->name.s << "\t" << seq2->seq.s << "\t";
      } else {
        o << "\t\t";
      }
      o << "splitat=";
      if (u1.empty()) {
        o << "(0,";
      } else {
        o << "(1,";
      }
      o << vsafe.back().second << ")\n";
      return;
    } else {
      j--;
    }
  }

  // can't split the read in two
  o << "CANTSPLIT\t";
  if (!u1.empty()) {
    printTranscripts(index,o,u1);
  }
  o << "\t";
  if (!u2.empty()) {
    printTranscripts(index,o,u2);
  }
  o << "\t" << seq1->name.s << "\t" << seq1->seq.s << "\t";
  if (paired) {
    o << seq2->name.s << "\t" << seq2->seq.s << "\t\n";
  } else {
    o << "\t\t\n";
  }
  return;

}


/** -- pseudobam functions -- **/
void outputPseudoBam(const KmerIndex &index, int ec, const kseq_t *seq1, const std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, const std::vector<std::pair<KmerEntry,int>> &v2, bool paired) {

  static char buf1[32768];
  static char buf2[32768];
  static char cig_[1000];
  char *cig = &cig_[0];
  
  
  if (seq1->name.l > 2 && seq1->name.s[seq1->name.l-2] == '/') {
    seq1->name.s[seq1->name.l-2] = 0; 
  }

  if (paired && seq2->name.l > 2 && seq2->name.s[seq2->name.l-2] == '/') {
    seq2->name.s[seq2->name.l-2] = 0; 
  }
  
  if (ec == -1) {
    // no mapping
    if (paired) {
      printf("%s\t77\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
      printf("%s\t141\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq2->name.s, seq2->seq.s, seq2->qual.s);
      //o << seq1->name.s << "" << seq1->seq.s << "\t" << seq1->qual.s << "\n";
      //o << seq2->name.s << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << seq2->seq.s << "\t" << seq2->qual.s << "\n";
    } else {
      printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
    }
  } else {
    if (paired) {

      int flag1 = 0x01 + 0x40;
      int flag2 = 0x01 + 0x80;
      
      if (v1.empty()) {
        flag1 += 0x04; // read unmapped
        flag2 += 0x08; // mate unmapped
      }

      if (v2.empty()) {
        flag1 += 0x08; // mate unmapped
        flag2 += 0x04; // read unmapped
      }

      if (!v1.empty() && !v2.empty()) {
        flag1 += 0x02; // proper pair
        flag2 += 0x02; // proper pair
      }
      

      int p1 = -1, p2 = -1;
      KmerEntry val1, val2;
      int nmap = index.ecmap[ec].size();
      Kmer km1, km2;
      
      if (!v1.empty()) {
        val1 = v1[0].first;
        p1 = v1[0].second;
        for (auto &x : v1) {
          if (x.second < p1) {
            val1 = x.first;
            p1 = x.second;
          }
        }
        km1 = Kmer((seq1->seq.s+p1));
      }

      if (!v2.empty()) {
        val2 = v2[0].first;
        p2 = v2[0].second;
        for (auto &x : v2) {
          if (x.second < p2) {
            val2 = x.first;
            p2 = x.second;
          }
        }
        km2 = Kmer((seq2->seq.s+p2));
      }

      bool revset = false;
      
      // output pseudoalignments for read 1
      for (auto tr : index.ecmap[ec]) {
        int f1 = flag1;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f1 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], seq1);
              revset = true;
            }
          }
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2 , val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f1 += 0x20; // mate reverse
          }
        }

        int posread = (f1 & 0x10) ? (x1.first - seq1->seq.l + 1) : x1.first;
        int posmate = (f1 & 0x20) ? (x2.first - seq2->seq.l + 1) : x2.first;

        getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, posmate, seq1->seq.l, index.target_lens_[tr]);
        int tlen = x2.first - x1.first;
        if (tlen != 0) {
          tlen += (tlen>0) ? 1 : -1;
        }
        
        printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", seq1->name.s, f1 & 0xFFFF, index.target_names_[tr].c_str(), posread, cig, posmate, tlen, (f1 & 0x10) ? &buf1[0] : seq1->seq.s, (f1 & 0x10) ? &buf2[0] : seq1->qual.s, nmap);
      }

      revset = false;
      // output pseudoalignments for read 2
      for (auto tr : index.ecmap[ec]) {
        int f2 = flag2;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f2 += 0x20; // mate reverse
          }
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2, val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f2 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], seq2);
              revset = true;
            }
            
          }
        }

        int posread = (f2 & 0x10) ? (x2.first - seq2->seq.l + 1) : x2.first;
        int posmate = (f2 & 0x20) ? (x1.first - seq1->seq.l + 1) : x1.first;

        getCIGARandSoftClip(cig, bool(f2 & 0x10), (f2 & 0x04) == 0, posread, posmate, seq2->seq.l, index.target_lens_[tr]);
        int tlen = x1.first - x2.first;
        if (tlen != 0) {
          tlen += (tlen > 0) ? 1 : -1;
        }

        printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", seq2->name.s, f2 & 0xFFFF, index.target_names_[tr].c_str(), posread, cig, posmate, tlen, (f2 & 0x10) ? &buf1[0] : seq2->seq.s,  (f2 & 0x10) ? &buf2[0] : seq2->qual.s, nmap);
      }
      
      
    } else {
      // single end
      int nmap = (int) index.ecmap[ec].size();
      KmerEntry val1 = v1[0].first;
      int p1 = v1[0].second;
      for (auto &x : v1) {
        if (x.second < p1) {
          val1 = x.first;
          p1 = x.second;
        }
      }
      Kmer km1 = Kmer((seq1->seq.s+p1));

      bool revset = false;

      for (auto tr : index.ecmap[ec]) {
        int f1 = 0;
        auto x1 = index.findPosition(tr, km1, val1, p1);

        if (!x1.second) {
          f1 += 0x10;
          if (!revset) {
            revseq(&buf1[0], &buf2[0], seq1);
            revset = true;
          }
        }

        int posread = (f1 & 0x10) ? (x1.first - seq1->seq.l+1) : x1.first;
        int dummy=1;
        getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, dummy, seq1->seq.l, index.target_lens_[tr]);
        
        printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d\n", seq1->name.s, f1 & 0xFFFF, index.target_names_[tr].c_str(), posread, cig, 0, 0, (f1 & 0x10) ? &buf1[0] : seq1->seq.s, (f1 & 0x10) ? &buf2[0] : seq1->qual.s, nmap);
      }
    }
  }
}


void revseq(char *b1, char *b2, const kseq_t *seq) {
  int n = seq->seq.l;
  char* s = seq->seq.s;
  b1[n] = 0;
  for (int i = 0; i < n; i++) {
    switch(s[i]) {
    case 'A': b1[n-1-i] = 'T'; break;
    case 'C': b1[n-1-i] = 'G'; break;
    case 'G': b1[n-1-i] = 'C'; break;
    case 'T': b1[n-1-i] = 'A'; break;
    default:  b1[n-1-i] = 'N';
    }
  }

  n = seq->qual.l;
  s = seq->qual.s;
  b2[n] = 0;
  for (int i = 0; i < n; i++) {
    b2[n-1-i] = s[i];
  }

  
}



void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength) {
  int softclip = 1 - posread;
  int overhang = (posread + length) - targetlength - 1;
  
  if (posread <= 0) {
    posread = 1;
  }

  if (mapped) {
    if (softclip > 0) {
      if (overhang > 0) {
        sprintf(cig, "%dS%dM%dS",softclip, (length-overhang - softclip), overhang);
      } else {
        sprintf(cig, "%dS%dM",softclip,length-softclip);
      }
    } else if (overhang > 0) {
      sprintf(cig, "%dM%dS", length-overhang, overhang);
    } else {
      sprintf(cig, "%dM",length);
    }
  } else {
    sprintf(cig, "*");
  }

  
  if (posmate <= 0) {
    posmate = 1;
  }
}
