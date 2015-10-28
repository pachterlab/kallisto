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

#include "ProcessReads.h"
#include "kseq.h"
#include "PseudoBam.h"


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

  int limit = 1048576;
  char *buf = new char[limit];
  std::vector<std::pair<const char*, int>> seqs;
  seqs.reserve(limit/50);

  SequenceReader SR(opt);

  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  //int tlencount = (opt.fld == 0.0) ? 10000 : 0;
  size_t numreads = 0;
  size_t nummapped = 0;

  bool paired = !opt.single_end;

  /*
  std::vector<std::pair<KmerEntry,int>> v1, v2;
  v1.reserve(1000);
  v2.reserve(1000);
  */


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

  if (opt.pseudobam) {
    index.writePseudoBamHeader(std::cout);
  }

  MasterProcessor MP(index, opt, tc);
  MP.processReads();
  numreads = MP.numreads;
  nummapped = MP.nummapped;

  std::cerr << " done" << std::endl;

  //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;

  if (opt.bias) {
    std::cerr << "[quant] learning parameters for sequence specific bias" << std::endl;
  }

  std::cerr << "[quant] processed " << pretty_num(numreads) << " reads, "
    << pretty_num(nummapped) << " reads pseudoaligned" << std::endl;

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


/** -- read processors -- **/

void MasterProcessor::processReads() {
  // start worker threads
  std::vector<std::thread> workers;
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(ReadProcessor(index,opt,tc,*this)));
  }
  // let the workers do their thing

  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }

  // now handle the modification of the mincollector
  for (auto &t : newECcount) {
    if (t.second <= 0) {
      continue;
    }
    int ec = tc.increaseCount(t.first); // modifies the ecmap

    if (ec != -1 && t.second > 1) {
      tc.counts[ec] += (t.second-1);
    }
  }

}

void MasterProcessor::update(const std::vector<int> &c, const std::vector<std::vector<int> > &newEcs,
                            int n, std::vector<int>& flens, std::vector<int>& bias) {
  // acquire the writer lock
  std::lock_guard<std::mutex> lock(this->writer_lock);

  for (int i = 0; i < c.size(); i++) {
    tc.counts[i] += c[i]; // add up ec counts
    nummapped += c[i];
  }

  for(auto &u : newEcs) {
    ++newECcount[u];
  }
  nummapped += newEcs.size();

  if (!flens.empty()) {
    int local_tlencount = 0;
    for (int i = 0; i < flens.size(); i++) {
      tc.flens[i] += flens[i];
      local_tlencount += flens[i];
    }
    tlencount += local_tlencount;
  }

  if (!bias.empty()) {
    int local_biasCount = 0;
    for (int i = 0; i < bias.size(); i++) {
      tc.bias5[i] += bias[i];
      local_biasCount += bias[i];
    }
    biasCount += local_biasCount;
  }

  numreads += n;
  // releases the lock
}

ReadProcessor::ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp) :
 paired(!opt.single_end), tc(tc), index(index), mp(mp) {
   // initialize buffer
   bufsize = 1ULL<<23;
   buffer = new char[bufsize];

   seqs.reserve(bufsize/50);
   newEcs.reserve(1000);
   counts.reserve((int) (tc.counts.size() * 1.25));
   clear();
}

ReadProcessor::~ReadProcessor() {
  if (buffer) {
      /*delete[] buffer;
    buffer = nullptr;*/
  }
}

void ReadProcessor::operator()() {
  while (true) {
    // grab the reader lock
    {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR.empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals,mp.opt.pseudobam);
      }
    } // release the reader lock

    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    mp.update(counts, newEcs, paired ? seqs.size()/2 : seqs.size(), flens, bias5);
    clear();
  }
}

void ReadProcessor::processBuffer() {
  // set up thread variables
  std::vector<std::pair<KmerEntry,int>> v1, v2;
  std::vector<int> vtmp;
  std::vector<int> u;

  u.reserve(1000);
  v1.reserve(1000);
  v2.reserve(1000);
  vtmp.reserve(1000);

  const char* s1 = 0;
  const char* s2 = 0;
  int l1,l2;

  bool findFragmentLength = (mp.opt.fld == 0) && (mp.tlencount < 10000);

  int flengoal = 0;
  flens.clear();
  if (findFragmentLength) {
    flengoal = (10000 - mp.tlencount);
    if (flengoal <= 0) {
      findFragmentLength = false;
      flengoal = 0;
    } else {
      flens.resize(tc.flens.size(), 0);
    }
  }

  int maxBiasCount = 0;
  bool findBias = mp.opt.bias && (mp.biasCount < mp.maxBiasCount);


  int biasgoal  = 0;
  bias5.clear();
  if (findBias) {
    biasgoal = (mp.maxBiasCount - mp.biasCount);
    if (biasgoal <= 0) {
      findBias = false;
    } else {
      bias5.resize(tc.bias5.size(),0);
    }
  }


  // actually process the sequences
  for (int i = 0; i < seqs.size(); i++) {
    s1 = seqs[i].first;
    l1 = seqs[i].second;
    if (paired) {
      i++;
      s2 = seqs[i].first;
      l2 = seqs[i].second;
    }

    numreads++;
    v1.clear();
    v2.clear();
    u.clear();

    // process read
    index.match(s1,l1, v1);
    if (paired) {
      index.match(s2,l2, v2);
    }

    // collect the target information
    int ec = -1;
    int r = tc.intersectKmers(v1, v2, !paired,u);
    if (u.empty()) {
      continue;
    } else {
      ec = tc.findEC(u);
    }

    /* --  possibly modify the pseudoalignment  -- */

    // If we have paired end reads where one end maps, check if some transcsripts
    // are not compatible with the mean fragment length
    if (paired && !u.empty() && (v1.empty() || v2.empty()) && tc.has_mean_fl) {
      vtmp.clear();
      // inspect the positions
      int fl = (int) tc.get_mean_frag_len();
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
        km = Kmer((s1+p));
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
        km = Kmer((s2+p));
      }

      // for each transcript in the pseudoalignment
      for (auto tr : u) {
        auto x = index.findPosition(tr, km, val, p);
        // if the fragment is within bounds for this transcript, keep it
        if (x.second && x.first + fl <= index.target_lens_[tr]) {
          vtmp.push_back(tr);
        }
        if (!x.second && x.first - fl >= 0) {
          vtmp.push_back(tr);
        }
      }

      if (vtmp.size() < u.size()) {
        u = vtmp; // copy
      }
    }

    // count the pseudoalignment
    if (ec == -1 || ec >= counts.size()) {
      // something we haven't seen before
      newEcs.push_back(u);
    } else {
      // add to count vector
      ++counts[ec];
    }



    /* -- collect extra information -- */
    // collect bias info
    if (findBias && !u.empty() && biasgoal > 0) {
      // collect sequence specific bias info
      if (tc.countBias(s1, (paired) ? s2 : nullptr, v1, v2, paired, bias5)) {
        biasgoal--;
      }
    }

    // collect fragment length info
    if (findFragmentLength && flengoal > 0 && paired && 0 <= ec &&  ec < index.num_trans && !v1.empty() && !v2.empty()) {
      // try to map the reads
      int tl = index.mapPair(s1, l1, s2, l2, ec);
      if (0 < tl && tl < flens.size()) {
        flens[tl]++;
        flengoal--;
      }
    }

    // pseudobam
    if (mp.opt.pseudobam) {
      if (paired) {
        outputPseudoBam(index, u,
          s1, names[i-1].first, quals[i-1].first, l1, names[i-1].second, v1,
          s2, names[i].first, quals[i].first, l2, names[i].second, v2,
          paired);
      } else {
        outputPseudoBam(index, u,
          s1, names[i].first, quals[i].first, l1, names[i].second, v1,
          nullptr, nullptr, nullptr, 0, 0, v2,
          paired);
      }
    }

    /*
    if (opt.verbose && numreads % 100000 == 0 ) {
      std::cerr << "[quant] Processed " << pretty_num(numreads) << std::endl;
    }*/
  }

}

void ReadProcessor::clear() {
  numreads=0;
  memset(buffer,0,bufsize);
  newEcs.clear();
  counts.clear();
  counts.resize(tc.counts.size(),0);
}





/** -- sequence reader -- **/
SequenceReader::~SequenceReader() {
  if (fp1) {
    gzclose(fp1);
  }
  if (paired && fp2) {
    gzclose(fp2);
  }

  kseq_destroy(seq1);
  if (paired) {
    kseq_destroy(seq2);
  }
}



// returns true if there is more left to read from the files
bool SequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals, bool full) {
  seqs.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
  int bufpos = 0;
  int pad = (paired) ? 2 : 1;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current file
        if(fp1) {
          gzclose(fp1);
        }
        if (paired && fp2) {
          gzclose(fp2);
        }
        // open the next one
        fp1 = gzopen(files[current_file].c_str(),"r");
        seq1 = kseq_init(fp1);
        l1 = kseq_read(seq1);
        state = true;
        if (paired) {
          current_file++;
          fp2 = gzopen(files[current_file].c_str(),"r");
          seq2 = kseq_init(fp2);
          l2 = kseq_read(seq2);
        }
      }
    }
    // the file is open and we have read into seq1 and seq2

    if (l1 > 0 && (!paired || l2 > 0)) {
      int bufadd = l1 + l2 + pad;
      // fits into the buffer
      if (full) {
        nl1 = seq1->name.l;
        if (paired) {
          nl2 = seq2->name.l;
        }
        bufadd += (l1+l2) + pad + (nl1+nl2)+ pad;
      }
      if (bufpos+bufadd< limit) {
        char *p1 = buf+bufpos;
        memcpy(p1, seq1->seq.s, l1+1);
        bufpos += l1+1;
        seqs.emplace_back(p1,l1);
        if (full) {
          p1 = buf+bufpos;
          memcpy(p1, seq1->qual.s,l1+1);
          bufpos += l1+1;
          quals.emplace_back(p1,l1);
          p1 = buf+bufpos;
          memcpy(p1, seq1->name.s,nl1+1);
          bufpos += nl1+1;
          names.emplace_back(p1,nl1);
        }

        if (paired) {
          char *p2 = buf+bufpos;
          memcpy(p2, seq2->seq.s,l2+1);
          bufpos += l2+1;
          seqs.emplace_back(p2,l2);
          if (full) {
            p2 = buf+bufpos;
            memcpy(p2,seq2->qual.s,l2+1);
            bufpos += l2+1;
            quals.emplace_back(p2,l2);
            p2 = buf + bufpos;
            memcpy(p2,seq2->name.s,nl2+1);
            bufpos += nl2+1;
            names.emplace_back(p2,nl2);
          }
        }
      } else {
        return true; // read it next time
      }

      // read for the next one
      l1 = kseq_read(seq1);
      if (paired) {
        l2 = kseq_read(seq2);
      }
    } else {
      current_file++; // move to next file
      state = false; // haven't opened file yet
    }
  }
}

bool SequenceReader::empty() {
  return (!state && current_file >= files.size());
}
