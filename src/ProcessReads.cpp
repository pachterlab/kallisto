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


void outputPseudoBam(const KmerIndex &index, int ec, const kseq_t *seq1, const std::vector<std::pair<KmerEntry,int>>& v1, const kseq_t *seq2, const std::vector<std::pair<KmerEntry,int>> &v2, bool paired);
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
  /*
  int biasCount = 0;
  const int maxBiasCount = 1000000;

  if (opt.pseudobam) {
    //index.writePseudoBamHeader(std::cout);
  }
  */

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
        mp.SR.fetchSequences(buffer, bufsize, seqs);
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
    if (paired && !u.empty() && (v1.empty() || v2.empty()) && flengoal == 0) {
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

    /*
    // Deal with strand specific reads
    if (opt.strand_specific && ec != -1 && !v1.empty()) {

      Kmer km;
      KmerEntry val = v1[0].first;
      int p = v1[0].second;
      for (auto &x : v1) {
        if (x.second < p) {
          val = x.first;
          p = x.second;
        }
      }
      km = Kmer((s1+p));

      std::vector<int> u = index.ecmap[ec];
      std::vector<int> v;
      v.reserve(u.size());

      for (auto tr : u)  {
        bool strand = ((km == km.rep()) == val.isFw());
        if (!strand) {
          v.push_back(tr);
        }
      }
      if (v.size() < u.size()) {
        tc.decreaseCount(ec);
        ec = tc.increaseCount(v);
      }
    }
    */

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
      //outputPseudoBam(index, ec, seq1, v1, seq2, v2, paired);
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
bool SequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs) {
  seqs.clear();
  int bufpos = 0;
  char * p1 = 0;
  char * p2 = 0;
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

      // fits into the buffer
      if (bufpos+l1+l2+pad < limit) {
        p1 = buf+bufpos;
        memcpy(p1, seq1->seq.s, l1+1);
        bufpos += l1+1;
        if (paired) {
          p2 = buf+bufpos;
          memcpy(p2, seq2->seq.s,l2+1);
          bufpos += l2+1;
        }
        seqs.emplace_back(p1,l1);
        if (paired) {
          seqs.emplace_back(p2,l2);
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

bool SequenceReader::fetchSequencesFull(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs, std::vector<std::pair<const char *, int> > &names, std::vector<std::pair<const char *, int> > &quals) {
  return false;
}




/** --- pseudobam functions -- **/

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
