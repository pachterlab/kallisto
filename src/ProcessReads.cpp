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
#include <limits>

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

// constants
const int auxlen = 7;

//methods

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

int ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt) {

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

  /*if (opt.pseudobam) {
    bam_hdr_t *t = createPseudoBamHeader(index);
    index.writePseudoBamHeader(std::cout);
  }*/

  
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
    MP.tc.write(of);
    of.close();
  }


  return numreads;
}



/** -- read processors -- **/

void MasterProcessor::processReads() {
  // open bamfile and fetch header
  if (opt.pseudobam) {
    std::string bamfn = opt.output + "/pseudoalignments.bam";
    bamh = createPseudoBamHeader(index);
    bamfp = sam_open(bamfn.c_str(), "wb");
    int r = sam_hdr_write(bamfp, bamh);
  }

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

  if (opt.pseudobam) {
    pseudobatchf_out.close();
  }
}

void MasterProcessor::processAln() {
  assert(opt.pseudobam);
  pseudobatchf_in.open(opt.output + "/pseudoaln.bin", std::ios::in | std::ios::binary);
  SR.reset();
  if (!opt.batch_mode) {
    std::vector<std::thread> workers;
    for (int i = 0; i < opt.threads; i++) {
      workers.emplace_back(std::thread(AlnProcessor(index,opt,*this)));
    }
    
    // let the workers do their thing
    for (int i = 0; i < opt.threads; i++) {
      workers[i].join(); //wait for them to finish
    }

  } else {
    assert(false); // not implemented yet
  }


  pseudobatchf_in.close();
}

void MasterProcessor::update(const std::vector<int>& c, const std::vector<std::vector<int> > &newEcs, 
                            std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, 
                            int n, std::vector<int>& flens, std::vector<int> &bias, const PseudoAlignmentBatch& pseudobatch, int id) {
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

  if (opt.pseudobam) {
    // copy the pseudo alignment information, either write to disk or queue up
    pseudobatch_stragglers.push_back(std::move(pseudobatch));
    while (true) {
      if (pseudobatch_stragglers.empty()) {
        break;
      }
      // find the smallest batch id
      auto min_it = std::min_element(pseudobatch_stragglers.begin(), pseudobatch_stragglers.end(), 
      [](const PseudoAlignmentBatch &p1, const PseudoAlignmentBatch &p2) -> bool {
        return p1.batch_id < p2.batch_id;
      });      
      if ((last_pseudobatch_id + 1) != min_it->batch_id) {
        break;
      }
      // if it is in sequence, write it out
      writePseudoAlignmentBatch(pseudobatchf_out, *min_it);
      // remove from processing
      pseudobatch_stragglers.erase(min_it); 
      last_pseudobatch_id += 1;
    }
  }

  numreads += n;
  // releases the lock
}

void MasterProcessor::writePseudoBam(const std::vector<bam1_t> &bv) {
  std::lock_guard<std::mutex> lock(this->writer_lock);

  for (const auto &b : bv) {
    int r = sam_write1(bamfp, bamh, &b);    
  }
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
   bufsize = mp.bufsize;
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
    int readbatch_id;
    // grab the reader lock
    if (mp.opt.batch_mode) {
      if (batchSR.empty()) {
        return;
      } else {
        batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, readbatch_id, false );
      }
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR.empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, readbatch_id, mp.opt.pseudobam || mp.opt.fusion);
      }
      // release the reader lock
    }
    pseudobatch.aln.clear();
    pseudobatch.batch_id = readbatch_id;
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    mp.update(counts, newEcs, ec_umi, new_ec_umi, paired ? seqs.size()/2 : seqs.size(), flens, bias5, pseudobatch, id);
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
      PseudoAlignmentInfo info;
      info.id = (paired) ? (i/2) : i; // read id
      info.paired = paired;
      if (!u.empty()) {
        info.r1empty = v1.empty();
        info.r2empty = v2.empty();
        KmerEntry val;
        info.k1pos = (!info.r1empty) ? findFirstMappingKmer(v1,val) : -1;
        info.k2pos = (!info.r2empty) ? findFirstMappingKmer(v2,val) : -1;
        
        
        if (ec != -1) {
          info.ec_id = ec;
        } else {
          info.ec_id = -1;
          info.u = u; // copy
        }
      }
      pseudobatch.aln.push_back(std::move(info));
      /*
      if (paired) {
        outputPseudoBam(index, u,
          s1, names[i-1].first, quals[i-1].first, l1, names[i-1].second, v1,
          s2, names[i].first, quals[i].first, l2, names[i].second, v2,
          paired, mp.bamh, mp.bamfp);
      } else {
        outputPseudoBam(index, u,
          s1, names[i].first, quals[i].first, l1, names[i].second, v1,
          nullptr, nullptr, nullptr, 0, 0, v2,
          paired, mp.bamh, mp.bamfp);
      }
      */
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


AlnProcessor::AlnProcessor(const KmerIndex& index, const ProgramOptions& opt, MasterProcessor& mp, int _id) :
 paired(!opt.single_end), index(index), mp(mp), id(_id) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];
   bambufsize = 1<<20;
   bambuffer = new char[bambufsize]; // refactor this?


   if (opt.batch_mode) {
     /* need to check this later */
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
   
   clear();
   
}


AlnProcessor::AlnProcessor(AlnProcessor && o) :
  paired(o.paired),
  index(o.index),
  mp(o.mp),
  id(o.id),
  bufsize(o.bufsize),
  bambufsize(o.bambufsize),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  umis(std::move(o.umis)),
  batchSR(std::move(o.batchSR)) {
    buffer = o.buffer;
    o.buffer = nullptr;
    o.bufsize = 0;
    bambuffer = o.bambuffer;
    o.bambuffer = nullptr;
    o.bambufsize = 0;
}

AlnProcessor::~AlnProcessor() {
  if (buffer != nullptr) {
    delete[] buffer;
    buffer = nullptr;
  }
  if (bambuffer != nullptr) {
    delete[] bambuffer;
    bambuffer = nullptr;
  }
}

void AlnProcessor::clear() {
  numreads=0;
  memset(buffer,0,bufsize);
  memset(bambuffer, 0, bambufsize);
  pseudobatch.aln.clear();
  pseudobatch.batch_id = -1;
}


void AlnProcessor::operator()() {
  while (true) {
    clear();
    int readbatch_id;
    // grab the reader lock
    if (mp.opt.batch_mode) {
      if (batchSR.empty()) {
        return;
      } else {
        batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, readbatch_id, true );
      }
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR.empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals, umis, readbatch_id, true);
        readPseudoAlignmentBatch(mp.pseudobatchf_in, pseudobatch);
        assert(pseudobatch.batch_id == readbatch_id);
        assert(pseudobatch.aln.size() == ((paired) ? seqs.size()/2 : seqs.size())); // sanity checks
      }
      // release the reader lock
    }
    // process our sequences
    processBuffer();


  }
}


void AlnProcessor::processBuffer() {

  /* something simple where we can construct the bam records */
  std::vector<bam1_t> bv;

  /*
  int num_pairs = pseudobatch.aln.size();
  int num_aln = 0;
  int bbsize = 0;
  for (int i = 0; i < n; i++) {
    int naln = 1, bsz = 0;
    const auto &x = pseudobatch.aln[i];
    if (x.ec_id != -1) {
      naln = index.ecmap[x.ec_id].size();
    } else if (!x.r1empty || !x.r2empty ) {
      naln = x.u.size();
    } 
    if (naln == 0) {
      naln = 1;
    }

    bsz = names[2*i].size() + names[2*i+1].size() + 8 + (seqs[2*i].size() + seqs[2*i+1].size() )/2 + 10;
    bbsize += bsz * naln;
  }

  std::cout << names[2*i].first << "\t" << x.id << "\t" << x.ec_id <<  "\t" << x.k1pos << "\t" << x.k2pos << "\n";

  bv.reserve(num_aln);
  if (bbsize * 2 < bambufsize) {
    char* t = new char[2*bbsize];
    memcpy(t, bambuffer, bambufsize);
    delete[] bambuffer;
    bambuffer = t;
    bambufsize = bbsize;
  } // let's hope this is enough for now
  
  char* buf = bambuffer;
  */
  int n = pseudobatch.aln.size();
  std::vector<int> u;
  

  Kmer km1,km2;
  KmerEntry val1, val2;

  for (int i = 0; i < n; i++) {
    bam1_t b1,b2, b1c, b2c;
    int si1 = (paired) ? 2*i : i;
    int si2 = (paired) ? 2*i +1 : -1;
    int rlen1 = seqs[si1].second;
    int rlen2;
    if (paired) {
      rlen2 = seqs[si2].second;
    }
    
    // fill in the bam core info
    b1.core.tid = -1;
    b1.core.pos = -1;
    b1.core.bin = 4680; // magic bin for unmapped reads
    b1.core.qual = 0;
    b1.core.flag = BAM_FUNMAP;
    b1.core.mtid = -1;
    b1.core.mpos = -1;
    b1.core.isize = 0;

    if (paired) {
      // fix this flag for b1
      b1.core.flag = BAM_FPAIRED | BAM_FREAD1 | BAM_FUNMAP | BAM_FMUNMAP;
      // set for b2
      b2.core.tid = -1;
      b2.core.pos = -1;
      b2.core.bin = 4680; // magic bin for unmapped reads
      b2.core.qual = 0;
      b2.core.flag = BAM_FPAIRED | BAM_FREAD2 | BAM_FUNMAP | BAM_FMUNMAP;      
      b2.core.mtid = -1;
      b2.core.mpos = -1;
      b2.core.isize = 0;
    }

    PseudoAlignmentInfo &pi = pseudobatch.aln[i];

    // fill in the data info
    fillBamRecord(b1, nullptr, seqs[si1].first, names[si1].first,  quals[si1].first, seqs[si1].second, names[si1].second, pi.r1empty);
    if (paired) {
      fillBamRecord(b2, nullptr, seqs[si2].first, names[si2].first,  quals[si2].first, seqs[si2].second, names[si2].second, pi.r2empty);
    }

    if (pi.r1empty && pi.r2empty) {
      bv.push_back(b1);
      if (paired) {
        bv.push_back(b2);
      }
    } else {
      u.clear();
      int ec = pi.ec_id;
      if (ec != -1) {
        u = index.ecmap[ec]; // copy, but meh
      } else {
        u = pi.u;
      }
      if (u.empty()) {
        // shouldn't happen
        bv.push_back(b1);
        if (paired) {
          bv.push_back(b2);
        }        
      } else {
        // set flags
        if (!pi.r1empty) {
          b1.core.flag &= ~BAM_FUNMAP;
          b2.core.flag &= ~BAM_FMUNMAP;
        }
        if (!pi.r2empty) {
          b1.core.flag &= ~BAM_FMUNMAP;
          b2.core.flag &= ~BAM_FUNMAP;
        }
        if (!pi.r1empty && !pi.r2empty) {
          b1.core.flag |= BAM_FPROPER_PAIR;
          b2.core.flag |= BAM_FPROPER_PAIR;
        }
        
        int32_t nmap = u.size();
        
        // set aux
        assert(b1.l_data + auxlen <= b1.m_data);
        assert(b2.l_data + auxlen <= b2.m_data);

        b1.data[b1.l_data] = 'N';
        b1.data[b1.l_data+1] = 'H';
        b1.data[b1.l_data+2] = 'i';
        memcpy(b1.data + b1.l_data + 3, &nmap, 4);
        b1.l_data += auxlen;

        b2.data[b2.l_data] = 'N';
        b2.data[b2.l_data+1] = 'H';
        b2.data[b2.l_data+2] = 'i';
        memcpy(b2.data + b2.l_data + 3, &nmap, 4);
        b2.l_data += auxlen;

        auto strandednessInfo = [&](Kmer km, KmerEntry& val, const std::vector<int> &u) -> std::pair<bool,bool> {
          bool reptrue = (km == km.rep());
          auto search = index.kmap.find(km.rep());
          if (search == index.kmap.end()) {
            return {false,reptrue};
          } else {
            val = search->second;
            if (val.contig == -1) {
              return {false,reptrue};
            } else {
              const Contig &c = index.dbGraph.contigs[val.contig];
              if (c.transcripts.empty()) {
                return {false,reptrue};
              }
              bool trsense = c.transcripts[0].sense;
              for (const auto & x : c.transcripts) {
                if (x.sense != trsense) {
                  if (std::find(u.begin(), u.end(), x.trid) != u.end()) {
                    return {false,reptrue};
                  }
                }
              }
              return {true, trsense == (reptrue == val.isFw())};
            }
          }
        };
        std::pair<bool,bool> strInfo1 = {true,true}, strInfo2 = {true,true};
        
        if (!pi.r1empty) {
          km1 = Kmer(seqs[si1].first + pi.k1pos);
          strInfo1 = strandednessInfo(km1, val1, u);

        }
        if (paired && !pi.r2empty) {
          km2 = Kmer(seqs[si2].first + pi.k2pos);
          strInfo2 = strandednessInfo(km2, val2, u);
        }
        
        
        // first read is all on reverse strand
        if (strInfo1.first && !strInfo1.second) {
          // reverse complement
          reverseComplementSeqInData(b1);
        }
        
        // we have second read and it is all on reverse strand
        if (paired && strInfo2.first && !strInfo2.second) {
          reverseComplementSeqInData(b2);
        }


        bool firstTr = true;
        for (auto t : u) {
          std::pair<int,bool> pos1, pos2;
          
          if (!pi.r1empty) {
            pos1 = index.findPosition(t, km1, val1, pi.k1pos);
          } else {
            pos1 = {std::numeric_limits<int>::min(), true};
          }
          

          if (paired) {
            if (!pi.r2empty) {
              pos2 = index.findPosition(t, km2, val2, pi.k2pos);
            } else {
              pos2 = {std::numeric_limits<int>::min(), true}; // use true so we don't reverse complement it
            }
          }


          // need to check the p value
          // check probabilities, discard 0.0

          // move data part of bamcore.
          b1c = b1;
          b1c.data = new uint8_t[b1c.m_data];
          memcpy(b1c.data, b1.data, b1c.m_data*sizeof(uint8_t));
          if (!strInfo1.first && !pos1.second) {
            reverseComplementSeqInData(b1c);
          }
          if (paired) {
            b2c = b2;    
            b2c.data = new uint8_t[b2c.m_data];
            memcpy(b2c.data, b2.data, b2c.m_data*sizeof(uint8_t));
            if (!strInfo2.first && !pos2.second) {
              reverseComplementSeqInData(b2c);
            }
          }
          
          
          // set strandedness
          if (paired) {
            // todo handle non mapping reads
            if (!pi.r1empty && !pos1.second) {
              b1c.core.flag |= BAM_FREVERSE;
              b2c.core.flag |= BAM_FMREVERSE;
            }
            if (!pi.r2empty && !pos2.second) {
              b1c.core.flag |= BAM_FMREVERSE;
              b2c.core.flag |= BAM_FREVERSE;
            }            
          } else {
            if (!pi.r1empty && !pos1.second) {
              b1c.core.flag |= BAM_FREVERSE;
            }
          }

          if (!firstTr) { // replace with highest prob
            b1c.core.flag |= BAM_FSECONDARY;
            b2c.core.flag |= BAM_FSECONDARY;
          }
          

          int softclip1, softclip2, overhang1, overhang2;
          int targetlength = index.target_lens_[t];
          // set core info
          b1c.core.tid = t;
          if (!pi.r1empty) {
            b1c.core.pos = (pos1.second) ? pos1.first-1 : pos1.first - rlen1;
            softclip1 = -b1c.core.pos;
            overhang1 = b1c.core.pos + rlen1 - targetlength;
            if (b1c.core.pos < 0) {
              b1c.core.pos = 0;
            }
            b1c.core.bin = hts_reg2bin(b1.core.pos, b1.core.pos + seqs[si1].second-1, 14, 5);
            b1c.core.qual = 255;
          }

          if (paired) {
            b2c.core.tid = t;
            if (!pi.r2empty) {
              b2c.core.pos = (pos2.second) ? pos2.first-1 : pos2.first - rlen2;
              softclip2 = -b2c.core.pos;
              overhang2 = b2c.core.pos + rlen2 - targetlength;
              if (b2c.core.pos < 0) {
                b2c.core.pos = 0;
              }

              b2c.core.bin = hts_reg2bin(b2.core.pos, b2.core.pos + seqs[si2].second, 14, 5);
              b2c.core.qual = 255;

              if (pi.r1empty) {
                b1c.core.pos = b2c.core.pos;
                b1c.core.bin = b2c.core.bin;
                b1c.core.qual = 0;
              }
            } else {
              b2c.core.pos = b1c.core.pos;
              b2c.core.bin = b1c.core.bin;
              b2c.core.qual = 0;
            }

            // fill in mate info
            b1c.core.mtid = t;
            b1c.core.mpos = b2c.core.pos;
            b2c.core.mtid = t;
            b2c.core.mpos = b1c.core.pos;

            if (!pi.r1empty && !pi.r2empty) {
              int tlen = pos2.first - pos1.first;
              if (tlen != 0) {
                tlen += (tlen>0) ? 1 : -1;
              }
              b1c.core.isize = tlen;
              b2c.core.isize = -tlen;
            }
          }

          if (!pi.r1empty) {
            if (softclip1 > 0 || overhang1 > 0) {
              fixCigarString(b1c, rlen1, softclip1, overhang1);
            }
          }
          if (!pi.r2empty) {
            if (softclip2 > 0 || overhang2 > 0) {
              fixCigarString(b2c, rlen2, softclip2, overhang2);
            }
          }

          if (!pi.r1empty || firstTr) {
            bv.push_back(b1c);
          }
          if (!pi.r2empty || firstTr) {
            bv.push_back(b2c);
          }

          firstTr = false;
        }        
      }
    }    
  }

  mp.writePseudoBam(bv);

  // clean up our mess
  for (auto &b : bv) {
    delete[] b.data;
  }
  bv.clear();
}

void fixCigarString(bam1_t &b, int rlen, int softclip, int overhang) {
  int ncig = 1;
  if (softclip > 0) {
    if (overhang > 0) {
      ncig = 3;
    } else {
      ncig = 2;
    }
  } else {
    ncig = 2;
  }
  if (ncig == 1) {
    return; // no-op
  }
  
  uint8_t *cigafter = b.data + b.core.l_qname + b.core.n_cigar * sizeof(uint32_t);
  int copylen = b.l_data - (b.core.l_qname + b.core.n_cigar * sizeof(uint32_t));
  uint8_t *dst = cigafter+(ncig-b.core.n_cigar)*sizeof(uint32_t);
  memmove(dst, cigafter, copylen);
  uint32_t *cig = (uint32_t*) (b.data + b.core.l_qname);
  b.l_data += (ncig - b.core.n_cigar)*sizeof(uint32_t);
  b.core.n_cigar = ncig;

  if (softclip > 0) {  
    if (overhang > 0) {       
      *cig = BAM_CSOFT_CLIP | (((uint32_t) softclip) << BAM_CIGAR_SHIFT);
      ++cig;
      *cig = BAM_CMATCH | (((uint32_t) (rlen - overhang - softclip)) << BAM_CIGAR_SHIFT);
      ++cig;
      *cig = BAM_CSOFT_CLIP | (((uint32_t) overhang) << BAM_CIGAR_SHIFT);      
    } else {
      *cig = BAM_CSOFT_CLIP | (((uint32_t) softclip) << BAM_CIGAR_SHIFT);
      ++cig;
      *cig = BAM_CMATCH | (((uint32_t) (rlen-softclip)) << BAM_CIGAR_SHIFT);
    }
  } else {
    *cig = BAM_CMATCH | (((uint32_t) (rlen-overhang)) << BAM_CIGAR_SHIFT);
    ++cig;
    *cig = BAM_CSOFT_CLIP | (((uint32_t) overhang) << BAM_CIGAR_SHIFT);
  }
}

void reverseComplementSeqInData(bam1_t &b) {

  const int bitrev[] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
  int slen = b.core.l_qseq;
  int lseq = (b.core.l_qseq+1)>>1;  
  uint8_t *s = b.data + (b.core.n_cigar << 2) + b.core.l_qname;
  
  for (int i = 0; i < (slen>>1); ++i) {
    int i1 = i>>1, i2=((slen-i-1)>>1);
    // let's just do the obvious thing and optimize later if needed
    int ca = s[i1] >> ((~i&1)<<2) & 0xf;
    int cb = s[i2] >> ((~(slen-1-i)&1)<<2) & 0xf;
    
    s[i1] &= 0xF << ((i&1)<<2);
    s[i1] |= bitrev[cb] << ((~i&1)<<2);
    s[i2] &= 0xF << (((slen-i-1)&1)<<2);
    s[i2] |= bitrev[ca] << ((~(slen-1-i)&1)<<2);
  }
  
  if ((slen & 1) == 1) {
    int ca = s[lseq>>1] >> 4 & 0xf;
    s[lseq>>1] &= 0xF;
    s[lseq>>1] |= bitrev[ca] << 4;
  }

  uint8_t *q = b.data + (b.core.n_cigar << 2) + b.core.l_qname + ((b.core.l_qseq+1)>>1);  
  for (int i = 0; i < (slen>>1); ++i) {
    uint8_t t = q[i];
    q[i] = q[slen-1-i];
    q[slen-1-i] = t;
  }
}

int fillBamRecord(bam1_t &b, uint8_t* buf, const char *seq, const char *name, const char *qual, int slen, int nlen, bool unmapped) {
 
  b.core.l_extranul =  (3 - (nlen % 4));  
  b.core.l_qname = nlen + b.core.l_extranul + 1;
  b.core.l_qseq = slen;

  if (buf == nullptr) {
    // allocate memory for buffer
    int blen = b.core.l_qname + 16 + ((b.core.l_qseq+2)>>1) + b.core.l_qseq + auxlen; // 16 for cigar, auxlen for NI:i:int
    b.data = new uint8_t[blen];
    b.m_data = blen;
    b.l_data = 0;
    buf = b.data;
  }

  memcpy(buf, name, nlen);
  int p = b.core.l_qname;
  for (int i = nlen; i < p; i++) {
    buf[i] = '\0';
  }

  // 99% of the time we just report a *, match. 1% there is some softclipping
  
  if (!unmapped) {
    b.core.n_cigar = 1;  
    uint32_t* cig = (uint32_t*) (buf+p);
    *cig = BAM_CMATCH | (((uint32_t) slen) << BAM_CIGAR_SHIFT);
    p += sizeof(uint32_t);
  } else {
    b.core.n_cigar = 0;    
  }
  
  
  // copy the sequence
  int lseq = (slen+1)>>1;
  uint8_t *seqp = (uint8_t *) (buf+p);
  memset(seqp,0,lseq);
  for (int i = 0; i < slen; ++i) {
    seqp[i>>1] |= seq_nt16_table[(int)seq[i]] << ((~i&1)<<2);
  }
  p += lseq;
  // copy qual
  for (int i = 0; i < slen; i++) {
    buf[p+i] = qual[i] - 33;
  }
  p += slen;

  b.l_data = p;
  return p;
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
  if (f_umi && f_umi->is_open()) {
    f_umi->close();
  }
}


void SequenceReader::reset() {
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

  if (f_umi && f_umi->is_open()) {
    f_umi->close();    
  }

  f_umi->clear();

  fp1 = 0;
  fp2 = 0;
  seq1 = 0;
  seq2 = 0;
  l1 = 0; 
  l2 = 0;
  nl1 = 0;
  nl2 = 0;
  current_file = 0;
  state = false;
  readbatch_id = -1;
}


// returns true if there is more left to read from the files
bool SequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals,
  std::vector<std::string> &umis, int& read_id,
  bool full) {
    
  std::string line;
  std::string umi;
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
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