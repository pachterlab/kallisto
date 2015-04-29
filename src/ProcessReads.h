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

template<typename Index, typename TranscriptCollector>
void ProcessReads(Index& index, const ProgramOptions& opt, TranscriptCollector& tc) {
  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  int tlencount = 10000;

  bool paired = (opt.files.size() == 2);

  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2;
  std::vector<std::pair<int,int>> v1, v2;
  v1.reserve(1000);
  v2.reserve(1000);

  int l1,l2; // length of read
  size_t nreads = 0;


  // for each file

  fp1 = gzopen(opt.files[0].c_str(), "r");
  seq1 = kseq_init(fp1);
  if (paired) {
    fp2 = gzopen(opt.files[1].c_str(),"r");
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

    nreads++;
    v1.clear();
    v2.clear();
    // process read
    index.match(seq1->seq.s, seq1->seq.l, v1);
    if (paired) {
      index.match(seq2->seq.s, seq2->seq.l, v2);
    }

    // collect the transcript information
    int ec = tc.collect(v1, v2, !paired);
    if (paired && 0 <= ec &&  ec < index.num_trans && tlencount > 0) {
      //bool allSame = true;
      bool allSame = (v1[0].first == ec && v2[0].first == ec) && (v1[0].second == 0 && v2[0].second == 0);

      if (allSame) {
        // try to map the reads
        int tl = index.mapPair(seq1->seq.s, seq1->seq.l, seq2->seq.s, seq2->seq.l, ec);
        if (0 < tl && tl < tc.flens.size()) {
          tc.flens[tl]++;
          tlencount--;
        }
      }
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
    
    if (opt.verbose && nreads % 100000 == 0 ) {
      std::cerr << "[quant] Processed " << nreads << std::endl;
    }
  }
  gzclose(fp1);
  if (paired) {
    gzclose(fp2);
  }

  kseq_destroy(seq1);
  if (paired) {
    kseq_destroy(seq2);
  }

  //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;
  
  // write output to outdir
  std::string outfile = opt.output + "/counts.txt"; 
  std::ofstream of;
  of.open(outfile.c_str(), std::ios::out);
  tc.write(of);
  of.close();
}


#endif // KALLISTO_PROCESSREADS_H
