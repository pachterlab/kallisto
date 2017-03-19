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
#include "PseudoBam.h"
#include "Fusion.hpp"

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


int findFirstMappingKmer(const std::vector<std::pair<KmerEntry,int>> &v,KmerEntry &val) {
  int p = -1;
  if (!v.empty()) {
    val = v[0].first;
    p = v[0].second;
    for (auto &x : v) {
      if (x.second < p) {
        val = x.first;
        p = x.second;
      }
    }
  }
  return p;
}

int ProcessBatchReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, std::vector<std::vector<int>> &batchCounts) {
  int limit = 1048576; 
  std::vector<std::pair<const char*, int>> seqs;
  seqs.reserve(limit/50);


  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  //int tlencount = (opt.fld == 0.0) ? 10000 : 0;
  size_t numreads = 0;
  size_t nummapped = 0;

  bool paired = !opt.single_end;
  
  if (paired) {
    std::cerr << "[quant] running in paired-end mode" << std::endl;
  } else {
    std::cerr << "[quant] running in single-end mode" << std::endl;
  }

  for (const auto& fs : opt.batch_files) { 
    for (int i = 0; i < fs.size(); i += (paired) ? 2 : 1) {
      if (paired) {
        std::cerr << "[quant] will process pair " << (i/2 +1) << ": "  << fs[i] << std::endl
                  << "                             " << fs[i+1] << std::endl;
      } else {
        std::cerr << "[quant] will process file " << i+1 << ": " << fs[i] << std::endl;
      }
    }
  }
  
  std::cerr << "[quant] finding pseudoalignments for all files ..."; std::cerr.flush();
  

  MasterProcessor MP(index, opt, tc);
  MP.processReads();
  numreads = MP.numreads;
  nummapped = MP.nummapped;
  batchCounts = std::move(MP.batchCounts);

  std::cerr << " done" << std::endl;

  if (opt.bias) {
    std::cerr << "[quant] learning parameters for sequence specific bias" << std::endl;
  }

  std::cerr << "[quant] processed " << pretty_num(numreads) << " reads, "
    << pretty_num(nummapped) << " reads pseudoaligned";
  if (nummapped == 0) {
    std::cerr << "[~warn] no reads pseudoaligned." << std::endl;
  }
  if (!opt.umi) {
    std::cerr << std::endl;
  } else {
    std::cerr << ", " << pretty_num(MP.num_umi) << " unique UMIs mapped" << std::endl;
  }

  return numreads;
  

}

int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc) {

  int limit = 1048576;
  std::vector<std::pair<const char*, int>> seqs;
  seqs.reserve(limit/50);

  //SequenceReader SR(opt);

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
  if (nummapped == 0) {
    std::cerr << "[~warn] no reads pseudoaligned." << std::endl;
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



/** -- read processors -- **/

void MasterProcessor::processReads() {
  // start worker threads
  if (!opt.batch_mode) {
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
  } else {
    std::vector<std::thread> workers;
    int num_ids = opt.batch_ids.size();
    int id =0;
    while (id < num_ids) {
      // TODO: put in thread pool
      workers.clear();
      int nt = std::min(opt.threads, (num_ids - id));
      for (int i = 0; i < nt; i++,id++) {
        workers.emplace_back(std::thread(ReadProcessor(index, opt, tc, *this, id)));
      }
      
      for (int i = 0; i < nt; i++) {
        workers[i].join();
      }
      
      if (opt.umi) {
        // process the regular EC umi now
        for (int i = 0; i < nt; i++) {
          int l_id = id - nt + i;
          auto &umis = batchUmis[l_id];
          std::sort(umis.begin(), umis.end());
          size_t sz = umis.size();
          nummapped += sz;
          if (sz > 0) {
            ++batchCounts[l_id][umis[0].first];
          }
          for (int j = 1; j < sz; j++) {
            if (umis[j-1] != umis[j]) {
              ++batchCounts[l_id][umis[j].first];
            }
          }
          umis.clear();
        }
      }
    }
    
    int num_newEcs = 0;
    if (!opt.umi) {      
      // for each cell
      for (int id = 0; id < num_ids; id++) {
        // for each new ec
        for (auto &t : newBatchECcount[id]) {
          // add the ec and count number of new ecs
          if (t.second <= 0) {
            continue;          
          }
          int ecsize = index.ecmap.size();
          int ec = tc.increaseCount(t.first);
          if (ec != -1 && ecsize < index.ecmap.size()) {          
            num_newEcs++; 
          }
        }
      }
      // for each cell
      for (int id = 0; id < num_ids; id++) {
        auto& c = batchCounts[id];
        c.resize(c.size() + num_newEcs,0);
        // for each new ec
        for (auto &t : newBatchECcount[id]) {
          // count the ec
          if (t.second <= 0) {
            continue;
          }
          int ec = tc.findEC(t.first);
          assert(ec != -1);
          ++c[ec];
        }
      }
    } else {
      // UMI case
      // for each cell
      for (int id = 0; id < num_ids; id++) {
        // for each new ec
        for (auto &t : newBatchECumis[id]) {
          // add the new ec
          int ecsize = index.ecmap.size();
          int ec = tc.increaseCount(t.first);
          if (ec != -1 && ecsize < index.ecmap.size()) {
            num_newEcs++;  
          }
        }
      }
      // for each cell
      for (int id = 0; id < num_ids; id++) {
        auto& c = batchCounts[id];
        c.resize(c.size() + num_newEcs,0);
        std::vector<std::pair<int, std::string>> umis;
        umis.reserve(newBatchECumis[id].size());
        // for each new ec
        for (auto &t : newBatchECumis[id]) {
          // record the ec,umi          
          int ec = tc.findEC(t.first);
          umis.push_back({ec, std::move(t.second)});
        }
        // find unique umis per ec
        std::sort(umis.begin(), umis.end());
        size_t sz = umis.size();
        if (sz > 0) {
          ++batchCounts[id][umis[0].first];
        }
        for (int j = 1; j < sz; j++) {
          if (umis[j-1] != umis[j]) {
            ++batchCounts[id][umis[j].first];
          }
        }
        for (auto x : c) {
          num_umi += x;
        }        
      }
    }
  }
}

void MasterProcessor::update(const std::vector<int>& c, const std::vector<std::vector<int> > &newEcs, 
                            std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, 
                            int n, std::vector<int>& flens, std::vector<int> &bias, int id) {
  // acquire the writer lock
  std::lock_guard<std::mutex> lock(this->writer_lock);

  if (!opt.batch_mode) {
    for (int i = 0; i < c.size(); i++) {
      tc.counts[i] += c[i]; // add up ec counts
      nummapped += c[i];
    }
  } else {
    if (!opt.umi) {
      for (int i = 0; i < c.size(); i++) {
        batchCounts[id][i] += c[i];
        nummapped += c[i];
      }
    } else {
      for (auto &t : ec_umi) {
        batchUmis[id].push_back(std::move(t));
      }
    }    
  }

  if (!opt.batch_mode) {
    for(auto &u : newEcs) {
      ++newECcount[u];
    }
  } else {
    if (!opt.umi) {
      for(auto &u : newEcs) {
        ++newBatchECcount[id][u];
      }
    } else {
      for (auto &u : new_ec_umi) {
        newBatchECumis[id].push_back(std::move(u));
      }
    }
  }
  if (!opt.umi) {
    nummapped += newEcs.size();
  } else {
    nummapped += new_ec_umi.size();
  }
  
  

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

void MasterProcessor::outputFusion(const std::stringstream &o) {
  std::string os = o.str();
  if (!os.empty()) {
    std::lock_guard<std::mutex> lock(this->writer_lock);
    ofusion << os << "\n";
  }
}


ReadProcessor::ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int _id) :
 paired(!opt.single_end), tc(tc), index(index), mp(mp), id(_id) {
   // initialize buffer
   bufsize = 1ULL<<23;
   buffer = new char[bufsize];

   if (opt.batch_mode) {
     assert(id != -1);
     batchSR.files = opt.batch_files[id];
     if (opt.umi) {
       batchSR.umi_files = {opt.umi_files[id]};
     }
     batchSR.paired = !opt.single_end;
   }

   seqs.reserve(bufsize/50);
   if (opt.umi) {
    umis.reserve(bufsize/50);
   }
   newEcs.reserve(1000);
   counts.reserve((int) (tc.counts.size() * 1.25));
   clear();
}

ReadProcessor::ReadProcessor(ReadProcessor && o) :
  paired(o.paired),
  tc(o.tc),
  index(o.index),
  mp(o.mp),
  id(o.id),
  bufsize(o.bufsize),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  umis(std::move(o.umis)),
  newEcs(std::move(o.newEcs)),
  flens(std::move(o.flens)),
  bias5(std::move(o.bias5)),
  batchSR(std::move(o.batchSR)),
  counts(std::move(o.counts)) {
    buffer = o.buffer;
    o.buffer = nullptr;
    o.bufsize = 0;
}

ReadProcessor::~ReadProcessor() {
  if (buffer != nullptr) {
      delete[] buffer;
      buffer = nullptr;
  }
}

void ReadProcessor::operator()() {
  while (true) {
    // grab the reader lock
    if (mp.opt.batch_mode) {
      if (batchSR.empty()) {
        return;
      } else {
        batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, false);
      }
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR.empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, mp.opt.pseudobam || mp.opt.fusion);
      }
      // release the reader lock
    }

    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    mp.update(counts, newEcs, ec_umi, new_ec_umi, paired ? seqs.size()/2 : seqs.size(), flens, bias5, id);
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
      if (mp.opt.fusion && !(v1.empty() || v2.empty())) {
        searchFusion(index,mp.opt,tc,mp,ec,names[i-1].first,s1,v1,names[i].first,s2,v2,paired);
      }
    } else {
      ec = tc.findEC(u);
    }


    /* --  possibly modify the pseudoalignment  -- */

    // If we have paired end reads where one end maps or single end reads, check if some transcsripts
    // are not compatible with the mean fragment length
    if (!mp.opt.umi && !u.empty() && (!paired || v1.empty() || v2.empty()) && tc.has_mean_fl) {
      vtmp.clear();
      // inspect the positions
      int fl = (int) tc.get_mean_frag_len();
      int p = -1;
      KmerEntry val;
      Kmer km;

      if (!v1.empty()) {
        p = findFirstMappingKmer(v1,val);
        km = Kmer((s1+p));
      }
      if (!v2.empty()) {
        p = findFirstMappingKmer(v2,val);
        km = Kmer((s2+p));
      }

      // for each transcript in the pseudoalignment
      for (auto tr : u) {
        auto x = index.findPosition(tr, km, val, p);
        // if the fragment is within bounds for this transcript, keep it
        if (x.second && x.first + fl <= index.target_lens_[tr]) {
          vtmp.push_back(tr);
        } else {
          //pass
        }
        if (!x.second && x.first - fl >= 0) {
          vtmp.push_back(tr);
        } else {
          //pass
        }
      }

      if (vtmp.size() < u.size()) {
        u = vtmp; // copy
      }
    }
    
    if (mp.opt.strand_specific && !u.empty()) {
      int p = -1;
      Kmer km;
      KmerEntry val;
      if (!v1.empty()) {
        vtmp.clear();
        bool firstStrand = (mp.opt.strand == ProgramOptions::StrandType::FR); // FR have first read mapping forward
        p = findFirstMappingKmer(v1,val);
        km = Kmer((s1+p));
        bool strand = (val.isFw() == (km == km.rep())); // k-mer maps to fw strand?
        // might need to optimize this
        const auto &c = index.dbGraph.contigs[val.contig];
        for (auto tr : u) {
          for (auto ctx : c.transcripts) {
            if (tr == ctx.trid) {
              if ((strand == ctx.sense) == firstStrand) {
                // swap out 
                vtmp.push_back(tr);
              } 
              break;
            }
          }          
        }
        if (vtmp.size() < u.size()) {
          u = vtmp; // copy
        }
      }
      
      if (!v2.empty()) {
        vtmp.clear();
        bool secondStrand = (mp.opt.strand == ProgramOptions::StrandType::RF);
        p = findFirstMappingKmer(v2,val);
        km = Kmer((s2+p));
        bool strand = (val.isFw() == (km == km.rep())); // k-mer maps to fw strand?
        // might need to optimize this
        const auto &c = index.dbGraph.contigs[val.contig];
        for (auto tr : u) {
          for (auto ctx : c.transcripts) {
            if (tr == ctx.trid) {
              if ((strand == ctx.sense) == secondStrand) {
                // swap out 
                vtmp.push_back(tr);
              } 
              break;
            }
          }          
        }
        if (vtmp.size() < u.size()) {
          u = vtmp; // copy
        }
      }
    }

    // find the ec
    if (!u.empty()) {
      ec = tc.findEC(u);

      if (!mp.opt.umi) {
        // count the pseudoalignment
        if (ec == -1 || ec >= counts.size()) {
          // something we haven't seen before
          newEcs.push_back(u);
        } else {
          // add to count vector
          ++counts[ec];
        }
      } else {       
        if (ec == -1 || ec >= counts.size()) {
          new_ec_umi.emplace_back(u, std::move(umis[i]));          
        } else {
          ec_umi.emplace_back(ec, std::move(umis[i]));
        }
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
  ec_umi.clear();
  new_ec_umi.clear();
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
  
  // check if umi stream is open, then close
}



// returns true if there is more left to read from the files
bool SequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals,
  std::vector<std::string> &umis, 
  bool full) {
    
  std::string line;
  std::string umi;
  
    
  seqs.clear();
  umis.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
   
  bool usingUMIfiles = !umi_files.empty();
  int umis_read = 0;
  
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
        // close current umi file
        if (usingUMIfiles) {
          // read up the rest of the files          
          f_umi->close();
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
        if (usingUMIfiles) {
          // open new umi file
          f_umi->open(umi_files[current_file]);          
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
        
        if (usingUMIfiles) {
          std::stringstream ss;
          std::getline(*f_umi, line);
          ss.str(line);
          ss >> umi;
          umis.emplace_back(std::move(umi));
        }
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

SequenceReader::SequenceReader(SequenceReader&& o) :
  fp1(o.fp1),
  fp2(o.fp2),
  seq1(o.seq1),
  seq2(o.seq2),
  l1(o.l1),
  l2(o.l2),
  nl1(o.nl1),
  nl2(o.nl2),
  paired(o.paired),
  files(std::move(o.files)),
  umi_files(std::move(o.umi_files)),
  f_umi(std::move(o.f_umi)),
  current_file(o.current_file),
  state(o.state) {
  o.fp1 = nullptr;
  o.fp2 = nullptr;
  o.seq1 = nullptr;
  o.seq2 = nullptr;
  o.state = false;
  
}