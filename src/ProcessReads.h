#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

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


void outputPseudoBam(const KmerIndex &index, int ec, const kseq_t *seq1, const std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, const std::vector<std::pair<KmerEntry,int>> &v2, bool paired);

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


void ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc) {
  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  int tlencount = 10000;
  size_t numreads = 0;
  size_t nummapped = 0;

  bool paired = !opt.single_end;

  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  std::vector<std::pair<KmerEntry,int>> v1, v2;
  v1.reserve(1000);
  v2.reserve(1000);

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
      if (ec != -1) {
        nummapped++;
        // collect sequence specific bias info
        if (opt.bias && biasCount < maxBiasCount) {
          if (tc.countBias(seq1->seq.s, seq2->seq.s, v1, v2, paired)) {
            biasCount++;
          }
        }
      }

      // collect fragment length info
      if (tlencount > 0 && paired && 0 <= ec &&  ec < index.num_trans) {
        //bool allSame = true;
        bool allSame = (index.dbGraph.ecs[v1[0].first.contig] == ec
                        && index.dbGraph.ecs[v2[0].first.contig] == ec)
          && (v1[0].second == 0 && v2[0].second == 0);

        if (allSame) {
          // try to map the reads
          int tl = index.mapPair(seq1->seq.s, seq1->seq.l, seq2->seq.s, seq2->seq.l, ec);
          if (0 < tl && tl < tc.flens.size()) {
            tc.flens[tl]++;
            tlencount--;
          }
        }
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
        std::cerr << "[quant] Processed " << numreads << std::endl;
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

  std::cerr << " done" << std::endl;

  //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;

  if (opt.bias) {
    std::cerr << "[quant] learning parameters for sequence specific bias" << std::endl;
  }
  
  std::cerr << "[quant] processed " << numreads << " reads, " << nummapped << " reads pseudoaligned" << std::endl;

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
}

void outputPseudoBam(const KmerIndex &index, int ec, const kseq_t *seq1, const std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, const std::vector<std::pair<KmerEntry,int>> &v2, bool paired) {
  if (ec == -1) {
    // no mapping
    if (paired) {
      printf("%s\t77\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq1->name.s, seq1->name.s, seq1->qual.s);
      printf("%s\t141\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq2->name.s, seq2->name.s, seq2->qual.s);
      //o << seq1->name.s << "" << seq1->seq.s << "\t" << seq1->qual.s << "\n";
      //o << seq2->name.s << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << seq2->seq.s << "\t" << seq2->qual.s << "\n";
    } else {
      printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", seq1->name.s, seq1->name.s, seq1->qual.s);
    }
  } else {
    
    auto getPos = [&](int tr, bool csense, KmerEntry val, int p) -> std::pair<int, bool> {
      const Contig &c = index.dbGraph.contigs[val.contig];
      int trpos = -1;
      bool trsense = true;
      
      for (auto x : c.transcripts) {
        if (x.trid == tr) {
          trpos = x.pos;
          trsense = x.sense;
          break;
        }
      }

      if (trpos == -1) {
        return {-1,true};
      }

      if (trsense) {
        if (csense) {
          return {trpos + val.getPos() - p + 1, csense}; // 1-based, case I
        } else {
          return {trpos + val.getPos() + index.k + p, csense}; // 1-based, case III
        }
      } else {
        if (csense) {
          return {trpos + (c.length - val.getPos() -1) + index.k + p, !csense};  // 1-based, case IV
        } else {
          return {trpos + (c.length - val.getPos())  - p, !csense}; // 1-based, case II
        }
      }
    };
    
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
      

      int p1 = -1, p2 = -1;
      bool fw1 = true, fw2 = true;
      bool csense1 = true, csense2 = true;
      KmerEntry val1, val2;
      
      if (!v1.empty()) {
        val1 = v1[0].first;
        p1 = v1[0].second;
        for (auto &x : v1) {
          if (x.second < p1) {
            val1 = x.first;
            p1 = x.second;
          }
        }
      
        Kmer km1 = Kmer((seq1->seq.s+p1));
        fw1 = (km1==km1.rep());
        csense1 = (fw1 == val1.isFw()); // is this in the direction of the contig?
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

        Kmer km2 = Kmer((seq2->seq.s+p2));
        fw2 = (km2==km2.rep());
        csense2 = (fw2 == val2.isFw()); // is this in the direction of the contig?
      }
      
      for (auto tr : index.ecmap[ec]) {
        int f1 = flag1;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = getPos(tr, csense1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
        }
        if (p2 != -1) {
          x2 = getPos(tr, csense2, val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
        }

        
        printf("%s\t%d\t%s\t%d\t%dM\t%d\t=\t%d\t%s\t%s\n", seq1->name.s, flag1, index.target_names_[tr].c_str(), x1.first, (int) seq1->seq.l, x2.first, (x2.first-x1.first), seq1->name.s, seq1->qual.s);
        //o << seq1->name.s << "\t?\t" << index.target_names_[tr] << "\t" << x1.first << "\t0\t" << seq1->seq.l << "M\t=\t" << x2.first << "\t" << (x2.first - x1.first) << "\t" << seq1->seq.s << "\t" << seq1->qual.s << "\n";
      }

      for (auto tr : index.ecmap[ec]) {
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = getPos(tr, csense1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
        }
        if (p2 != -1) {
          x2 = getPos(tr, csense2, val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
        }
        // need to figure out flag business
        printf("%s\t%d\t%s\t%d\t%dM\t%d\t=\t%d\t%s\t%s\n", seq2->name.s, flag2, index.target_names_[tr].c_str(), x2.first, (int) seq2->seq.l, x1.first, (x1.first-x2.first), seq2->name.s, seq2->qual.s);
        //o << seq2->name.s << "\t?\t" << index.target_names_[tr] << "\t" << x2.first << "\t0\t" << seq1->seq.l << "M\t=\t" << x1.first << "\t" << (x1.first - x2.first) << "\t" << seq2->seq.s << "\t" << seq2->qual.s << "\n";
      }
      
      
    } else {
      KmerEntry val1 = v1[0].first;
      int p1 = v1[0].second;
      for (auto &x : v1) {
        if (x.second < p1) {
          val1 = x.first;
          p1 = x.second;
        }
      }
      
      Kmer km1 = Kmer((seq1->seq.s+p1));
      bool fw1 = (km1==km1.rep());
      bool csense1 = (fw1 == val1.isFw()); // is this in the direction of the contig?
      
    }
    
    
  }
}


#endif // KALLISTO_PROCESSREADS_H
