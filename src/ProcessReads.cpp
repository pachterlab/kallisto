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

#include <iomanip>

#include "ProcessReads.h"
#include "kseq.h"
#include "PseudoBam.h"
#include "Fusion.hpp"
#include "BUSData.h"
#include "BUSTools.h"
#include <htslib/kstring.h>


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
const int default_trans_auxlen = 14; // for NI:i:int and ZW:f:0.0
const int default_genome_auxlen = 7; // for ZW:f:0.0

//methods

int64_t ProcessBatchReads(MasterProcessor& MP, const ProgramOptions& opt) {
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
  
  MP.processReads();
  numreads = MP.numreads;
  nummapped = MP.nummapped;

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

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt) {

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
  if (opt.verbose) {
    std::cerr << std::endl;
  }

  /*if (opt.pseudobam) {
    bam_hdr_t *t = createPseudoBamHeader(index);
    index.writePseudoBamHeader(std::cout);
  }*/

  
  MP.processReads();
  numreads = MP.numreads;
  nummapped = MP.nummapped;
  if (opt.verbose) {
    std::cerr << std::endl << "[quant] done " << std::endl;
  } else {
    std::cerr << " done" << std::endl;
  }

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


int64_t ProcessBUSReads(MasterProcessor& MP, const  ProgramOptions& opt) {

  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  //int tlencount = (opt.fld == 0.0) ? 10000 : 0;
  size_t numreads = 0;
  size_t nummapped = 0;
  bool paired = !opt.single_end;

  for (int i = 0, si=1; i < opt.files.size(); si++) {
    auto& busopt = opt.busOptions;
    std::cerr << "[quant] will process sample " << si<< ": ";
    for (int j = 0; j < busopt.nfiles; j++,i++) {
      //              "[quant] will process sample 1: "
      if (j>0) {
        std::cerr << "                               ";
      }
      std::cerr << opt.files[i] << std::endl;
    }
  }

  // for each file
  std::cerr << "[quant] finding pseudoalignments for the reads ..."; std::cerr.flush();

  if (opt.genomebam) {
    /*
    // open bam files for writing
    MP.bamh = createPseudoBamHeaderGenome(MP.model);
    MP.bamfps = new htsFile*[MP.numSortFiles];
    for (int i = 0; i < MP.numSortFiles; i++) {
      MP.bamfps[i] = sam_open((opt.output + "/tmp." + std::to_string(i) + ".bam").c_str(), "wb1");
      int r = sam_hdr_write(MP.bamfps[i], MP.bamh);
    }

    // assign breakpoints to chromosomes
    MP.breakpoints.clear();
    MP.breakpoints.assign(MP.numSortFiles -1 , (((uint64_t)-1)<<32));
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> chrWeights(MP.model.chr.size());
    for (const auto& t : MP.model.genes) {
      if (t.chr != -1 && t.stop > 0) {
        chrWeights[t.chr].push_back({t.stop, 1});
      }      
    }

    double sum = 0;
    for (const auto &chrw : chrWeights) {
      for (const auto &p : chrw) {
        sum += p.second;
      }
    }
    double bpLimit = sum / (MP.numSortFiles-1);

    for (auto& chrw : chrWeights) {
      // sort each by stop point
      std::sort(chrw.begin(), chrw.end(), [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t> b) { return a.first < b.first;});
    }

    double bp = 0.0;
    int ifile = 0;
    for (int i = 0; i < chrWeights.size(); i++) {
      const auto &chrw = chrWeights[i];
      for (const auto &x : chrw) {
        bp += x.second;
        if (bp > bpLimit) {
          uint64_t pos = ((uint64_t) i) << 32 | ((uint64_t) x.first+1) << 1;
          MP.breakpoints[ifile] = pos;
          ifile++;
          while (bp > bpLimit) {
            bp -= bpLimit;
          }
        }
      }
    }
    */
  }




  
  MP.processReads();
  numreads = MP.numreads;
  nummapped = MP.nummapped;
  if (opt.verbose) {
    std::cerr << std::endl << "[quant] done " << std::endl;
  } else {
    std::cerr << " done" << std::endl;
  }

  //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;


  std::cerr << "[quant] processed " << pretty_num(numreads) << " reads, "
    << pretty_num(nummapped) << " reads pseudoaligned" << std::endl;
  if (nummapped == 0) {
    std::cerr << "[~warn] no reads pseudoaligned." << std::endl;
  }
  
  return numreads;
}


/** -- read processors -- **/

void MasterProcessor::processReads() {
  

  // start worker threads
  if (!opt.batch_mode && !opt.bus_mode) {
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
  } else if (opt.bus_mode) {
    std::vector<std::thread> workers;
    for (int i = 0; i < opt.threads; i++) {
      workers.emplace_back(std::thread(BUSProcessor(index,opt,tc,*this)));
    }
    
    // let the workers do their thing
    for (int i = 0; i < opt.threads; i++) {
      workers[i].join(); //wait for them to finish
    }

    // now handle the modification of the mincollector
    for (int i = 0; i < bus_ecmap.size(); i++) {
      auto &u = bus_ecmap[i];
      int ec = index.ecmapinv.size();
      auto it = bus_ecmapinv.find(u);
      if (it->second != ec) {
        std::cout << "Error" << std::endl;
        exit(1);
      }      
      index.ecmapinv.insert({u,ec});
      index.ecmap.push_back(u);
    }

  
  } else if (opt.batch_mode) {    
    std::vector<std::thread> workers;
    int num_ids = opt.batch_ids.size();
    int id =0;
    tmp_bc.assign(opt.threads, {});
    while (id < num_ids) {
      // TODO: put in thread pool
      workers.clear();
      int nt = std::min(opt.threads, (num_ids - id));
      int first_id = id;
      for (int i = 0; i < nt; i++,id++) {
        tmp_bc[i].assign(tc.counts.size(), 0);     
        workers.emplace_back(std::thread(ReadProcessor(index, opt, tc, *this, id,i)));
      }
      
      for (int i = 0; i < nt; i++) {
        workers[i].join();
      }
      
      if (!opt.umi) {
        for (int i = 0; i < nt; i++) {
          auto &bc = batchCounts[first_id+i];
          auto &tmp_c = tmp_bc[i];
          for (int j = 0; j < tmp_c.size(); j++) {
            if (tmp_c[j] > 0) {
              bc.push_back({j, tmp_c[j]});
            }
          }
        }
      } else {
        std::vector<int32_t> tmp_counts;
        tmp_counts.resize(tc.counts.size(), 0);
        // process the regular EC umi now
        for (int i = 0; i < nt; i++) {
          tmp_counts.assign(tmp_counts.size(), 0);
          int l_id = id - nt + i;
          auto &umis = batchUmis[l_id];
          std::sort(umis.begin(), umis.end());
          size_t sz = umis.size();
          nummapped += sz;
          if (sz > 0) {
            ++tmp_counts[umis[0].first];
          }
          for (int j = 1; j < sz; j++) {
            if (umis[j-1] != umis[j]) {
              ++tmp_counts[umis[j].first];
            }
          }
          auto &bc = batchCounts[l_id];
          bc.clear();
          for (int j = 0; j < tmp_counts.size(); j++) {
            if (tmp_counts[j] > 0) {
              bc.push_back({j, tmp_counts[j]});
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
      std::vector<int> tmp_counts;
      tmp_counts.resize(tc.counts.size(), 0);
      for (int id = 0; id < num_ids; id++) {
        tmp_counts.assign(tmp_counts.size(), 0);
        // for each new ec
        for (auto &t : newBatchECcount[id]) {
          // count the ec
          if (t.second <= 0) {
            continue;
          }
          int ec = tc.findEC(t.first);
          assert(ec != -1);
          tmp_counts[ec] += t.second;
        }
        auto& bc = batchCounts[id];
        for (int j = 0; j < tmp_counts.size(); j++) {
          if (tmp_counts[j] > 0) {
            bc.push_back({j,tmp_counts[j]});
          }
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
      
      std::vector<int> tmp_counts;
      tmp_counts.resize(tc.counts.size(), 0);
      // for each cell
      for (int id = 0; id < num_ids; id++) {
        tmp_counts.assign(tmp_counts.size(), 0);
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
          ++tmp_counts[umis[0].first];
        }
        for (int j = 1; j < sz; j++) {
          if (umis[j-1] != umis[j]) {
            ++tmp_counts[umis[j].first];
          }
        }
        
        auto& bc = batchCounts[id];
        for (int j = 0; j < tmp_counts.size(); j++) {
          if (tmp_counts[j] > 0) {
            bc.push_back({j,tmp_counts[j]});
          }
        }
        for (auto x : bc) {
          num_umi += x.second;
        }
      }
    }
  }

  if (opt.pseudobam) {
    pseudobatchf_out.close();
  }
  if (opt.bus_mode) {
    busf_out.close();
  }
}

void MasterProcessor::processAln(const EMAlgorithm& em, bool useEM = true) {
  // open bamfile and fetch header
  std::string bamfn = opt.output + "/pseudoalignments.bam";
  if (opt.pseudobam) {
    if (opt.genomebam) {
      bamh = createPseudoBamHeaderGenome(model);

      bamfps = new htsFile*[numSortFiles];
      for (int i = 0; i < numSortFiles; i++) {
        bamfps[i] = sam_open((opt.output + "/tmp." + std::to_string(i) + ".bam").c_str(), "wb1");
        int r = sam_hdr_write(bamfps[i], bamh);
      }

    } else {
      bamh = createPseudoBamHeaderTrans(index);
      bamfp = sam_open(bamfn.c_str(), "wb");
      int r = sam_hdr_write(bamfp, bamh);    
    }
    
    if (opt.threads > 1 && !opt.genomebam) {
      // makes no sens to use threads on unsorted bams
      hts_set_threads(bamfp, opt.threads);
    }
  }
  std::cerr << "[  bam] writing pseudoalignments to BAM format .. "; std::cerr.flush();

  // figure out where to place breakpoints!
  breakpoints.clear();
  breakpoints.assign(numSortFiles -1 , (((uint64_t)-1)<<32));  
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> chrWeights(model.chr.size());
  for (const auto& t : model.transcripts) {
    if (t.id >= 0 && t.id < index.num_trans) {
      // valid id
      if (t.chr != -1 && t.stop > 0) {
        chrWeights[t.chr].push_back({t.stop, ((int32_t) (em.alpha_[t.id]+1))});
      }
    }
  }

  double sum = 0;
  for (const auto &chrw : chrWeights) {
    for (const auto &p : chrw) {
      sum += p.second;
    }
  }
  
  double bpLimit = sum / (numSortFiles-1);

  // place breakpoints on evenly spaced   
  for (auto& chrw : chrWeights) {
    // sort each by stop point
    std::sort(chrw.begin(), chrw.end(), [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t> b) { return a.first < b.first;});
  }

  double bp = 0.0;
  int ifile = 0;
  for (int i = 0; i < chrWeights.size(); i++) {
    const auto &chrw = chrWeights[i];
    for (const auto &x : chrw) {
      bp += x.second;
      if (bp > bpLimit) {
        uint64_t pos = ((uint64_t) i) << 32 | ((uint64_t) x.first+1) << 1;
        breakpoints[ifile] = pos;
        ifile++;
        while (bp > bpLimit) {
          bp -= bpLimit;
        }
      }
    }
  }

  // debug and show breakpoints
  /*
  std::cout << "bplimit = " << bpLimit << ", sum = " << sum << std::endl;
  std::cout << "breakpoints (" << breakpoints.size() << ") = {" << std::endl;
  for (auto &x : breakpoints) {
    int32_t tid = (int32_t) (x >> 32);
    int32_t pos = (int32_t) ((x >> 1) & (0xFFFF)) - 1;
    std::string chr = "*";
    if (tid != -1) {
      chr = model.chr[tid].name;
    }
    std::cout << "  " << chr << " (" << model.chr[tid].len << ") : " << pos << std::endl;
  }

  std::cout << "}" << std::endl;
  */

  assert(opt.pseudobam);
  pseudobatchf_in.open(opt.output + "/pseudoaln.bin", std::ios::in | std::ios::binary);
  SR->reset();

  std::vector<std::thread> workers;
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(AlnProcessor(index,opt,*this, em, model, useEM)));
  }
  
  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }

  pseudobatchf_in.close();
  remove((opt.output + "/pseudoaln.bin").c_str());
  std::cerr << "done" << std::endl;

  if (opt.genomebam) {
    std::cerr << "[  bam] sorting BAM files .. "; std::cerr.flush();
    for (int i = 0; i < numSortFiles; i++) {
      sam_close(bamfps[i]);
      bamfps[i] = nullptr;
    }

    // at this point we don't need the index
    index.clear();
    bamfp = sam_open(bamfn.c_str(), "wb9");
    int r = sam_hdr_write(bamfp, bamh);
    if (opt.threads > 1) {
      // makes no sense to use threads on unsorted bams
      hts_set_threads(bamfp, opt.threads);
    }

    int ret;
    std::vector<bam1_t> bv;
    std::vector<std::pair<uint64_t, uint64_t>> bb;
    bam1_t b;  
    int tid;
    BGZF* ofp = bamfp->fp.bgzf;
    for (int i = 0; i < numSortFiles; i++) {
      bv.clear();
      bb.clear();
      std::string tmpFileName = opt.output + "/tmp." + std::to_string(i) + ".bam";
      bamfps[i] =  sam_open(tmpFileName.c_str(), "rb");

      bam_hdr_t *tmp_hdr = sam_hdr_read(bamfps[i]); // unused results

      // init
      memset(&b, 0, sizeof(b));
      
      if (i < numSortFiles-1) {
        // sort the vector
        while ((ret = bam_read1(bamfps[i]->fp.bgzf, &b)) >= 0) {        
          bv.push_back(b);  
          memset(&b, 0, sizeof(b));        
        }

        sam_close(bamfps[i]);
        bamfps[i] = nullptr;
        remove(tmpFileName.c_str());

        uint64_t n = bv.size();
        for (uint64_t j = 0; j < n; j++) {
          b = bv[j];          
          uint64_t pos = ((uint64_t) b.core.tid) << 32 | ((uint64_t) b.core.pos+1) << 1 | (b.core.flag & BAM_FREVERSE) >>4;
          bb.push_back({pos,j});
        }
        
        std::sort(bb.begin(), bb.end(), [](std::pair<uint64_t, uint64_t> a, std::pair<uint64_t, uint64_t> b) {return a.first < b.first;});

        for (auto &x : bb) {
          b = bv[x.second];           
          ret = bam_write1(ofp, &b);               
          free(bv[x.second].data);
          bv[x.second].l_data = 0;
          bv[x.second].m_data = 0;
          
        }
      } else {
        // for unsorted files, just copy directly
        memset(&b, 0, sizeof(b));
        while ((ret = bam_read1(bamfps[i]->fp.bgzf, &b)) >= 0) {        
          ret = bam_write1(ofp, &b);
        }  
        sam_close(bamfps[i]);
        bamfps[i] = nullptr;
        remove(tmpFileName.c_str());
      }
    }

    sam_close(bamfp);
    bamfp = nullptr;    
    std::cerr << "done" << std::endl;
    // if we are multithreaded we need to construct the index last
    
    std::cerr << "[  bam] indexing BAM file .. "; std::cerr.flush();

    ret = sam_index_build3(bamfn.c_str(), (bamfn+".bai").c_str(), 0, opt.threads);
    if (ret != 0) {
      std::cerr << " invalid return code when indexing file " << ret << " .. ";
    }
    std::cerr << "done" << std::endl;
    
  }
}

void MasterProcessor::update(const std::vector<int>& c, const std::vector<std::vector<int> > &newEcs, 
                            std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, 
                            int n, std::vector<int>& flens, std::vector<int> &bias, const PseudoAlignmentBatch& pseudobatch, std::vector<BUSData> &bv, std::vector<std::pair<BUSData, std::vector<int32_t>>> newBP, int *bc_len, int *umi_len,  int id, int local_id) {
  // acquire the writer lock
  std::lock_guard<std::mutex> lock(this->writer_lock);

  if (!opt.batch_mode) {
    for (int i = 0; i < c.size(); i++) {
      tc.counts[i] += c[i]; // add up ec counts
      nummapped += c[i];
    }
  } else {
    if (!opt.umi) {
      auto& bc = tmp_bc[local_id];
      for (int i = 0; i < c.size(); i++) {
        bc[i] += c[i];
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
    if (opt.batch_mode) {
      auto &bflen = batchFlens[id];
      auto &tcount = tlencounts[id];
      for (int i = 0; i < flens.size(); i++) {
        bflen[i] += flens[i];
        tcount += flens[i];
      }
    } else {
      int local_tlencount = 0;
      for (int i = 0; i < flens.size(); i++) {
        tc.flens[i] += flens[i];
        local_tlencount += flens[i];
      }
      tlencount += local_tlencount;
    }
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

  if (opt.bus_mode) {
    int bus_bc_sum = 0;
    int bus_umi_sum = 0;
    for (int i = 0; i <= 32; i++) {
      bus_bc_sum += bus_bc_len[i];
      bus_umi_sum += bus_umi_len[i];
    }

    if (bus_bc_sum < 10000 or bus_umi_sum < 10000) {
      for (int i = 0; i < 32; i++) {
        bus_bc_len[i] += bc_len[i];
        bus_umi_len[i] += umi_len[i];
      }
    }

    // add new equiv classes to extra format
    int offset = index.ecmapinv.size();
    for (auto &bp : newBP) {
      auto& u = bp.second;
      int ec = -1;
      auto it = bus_ecmapinv.find(u);
      if (it != bus_ecmapinv.end()) {
        ec = it->second;
      } else {
        ec = offset + bus_ecmapinv.size();
        bus_ecmapinv.insert({u,ec});
        bus_ecmap.push_back(u);
      }
      auto &b = bp.first;
      b.ec = ec;
      bv.push_back(b);
    }

    //copy bus mode information, write to disk or queue up
    writeBUSData(busf_out, bv); 
    /*for (auto &bp : newBP) {
      newB.push_back(std::move(bp));
    } */
  }

  numreads += n;
  // releases the lock
}

void MasterProcessor::writePseudoBam(const std::vector<bam1_t> &bv) {
  std::lock_guard<std::mutex> lock(this->writer_lock);
  // locking is handled by htslib
  //kstring_t str = { 0, 0, NULL };
  for (const auto &b : bv) {
    /*
    std::cout << "name: " <<  b.data << ", ldata:" << (int)b.l_data << " ,mdata: " << (int) b.m_data  
      <<  ", lqname " << (int) b.core.l_qname <<", lqseq " <<  (int) b.core.l_qseq << ", enull " << (int) b.core.l_extranul << std::endl;
    for (int i = 0; i < b.l_data; ++i) {      
      std::cout << (int) b.data[i] << " ";
    }
    std::cout << std::endl;    
    std::cout << "tid: " << b.core.tid << ", pos: " << (int) b.core.pos
              << ", flag: " << b.core.flag << ", ncigar: " << b.core.n_cigar << ", mtid: " <<(int) b.core.mtid << ", mpos: " << (int) b.core.mpos << ", isize" << (int) b.core.isize <<  std::endl;
    kstring_t str = { 0, 0, NULL };
    sam_format1(bamh, &b, &str);
    kputc('\n', &str);
    std::cout << str.s << std::endl;
    */
    int r = sam_write1(bamfp, bamh, &b);    
  } 
}

// bvv needs to be sorted according to MasterProcessor::bucketSplits. 
void MasterProcessor::writeSortedPseudobam(const std::vector<std::vector<bam1_t>> &bvv) { 
  assert(bvv.size() == numSortFiles);
  
  for (int i = 0; i < numSortFiles; i++) {
    std::lock_guard<std::mutex> lock(this->writer_lock);
    for (const auto &b : bvv[i]) {
      int r = sam_write1(bamfps[i], bamh, &b);
    }
  }
}

void MasterProcessor::outputFusion(const std::stringstream &o) {
  std::string os = o.str();
  if (!os.empty()) {
    std::lock_guard<std::mutex> lock(this->writer_lock);
    ofusion << os << "\n";
  }
}


ReadProcessor::ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int _id, int _local_id) :
 paired(!opt.single_end), tc(tc), index(index), mp(mp), id(_id), local_id(_local_id) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];

   if (opt.batch_mode) {
     assert(id != -1);
     batchSR.files = opt.batch_files[id];
     batchSR.nfiles = opt.batch_files[id].size();
     batchSR.reserveNfiles(opt.batch_files[id].size());
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
  local_id(o.local_id),
  bufsize(o.bufsize),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)),
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
        batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, flags, umis, readbatch_id, mp.opt.pseudobam );
      }
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, umis, readbatch_id, mp.opt.pseudobam || mp.opt.fusion);
      }
      // release the reader lock
    }
    pseudobatch.aln.clear();
    pseudobatch.batch_id = readbatch_id;
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    std::vector<BUSData> tmp_v{};
    mp.update(counts, newEcs, ec_umi, new_ec_umi, paired ? seqs.size()/2 : seqs.size(), flens, bias5, pseudobatch, tmp_v, std::vector<std::pair<BUSData, std::vector<int32_t>>>{}, nullptr, nullptr, id, local_id);
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
  if (mp.opt.batch_mode) {
    findFragmentLength = (mp.opt.fld == 0) && (mp.tlencounts[id] < 10000);
  }

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
    if (!mp.opt.single_overhang && !mp.opt.umi && !u.empty() && (!paired || v1.empty() || v2.empty()) && tc.has_mean_fl) {
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
    
    if (mp.opt.verbose && numreads > 0 && numreads % 1000000 == 0 ) {   
      int nmap = mp.nummapped;
      for (int i = 0; i < counts.size(); i++) {
        nmap += counts[i];
      }
      nmap += newEcs.size();

      std::cerr << '\r' << (numreads/1000000) << "M reads processed (" 
        << std::fixed << std::setw( 3 ) << std::setprecision( 1 ) << ((100.0*nmap)/double(numreads))
        << "% pseudoaligned)"; std::cerr.flush();
    }
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



BUSProcessor::BUSProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int _id, int _local_id) :
 paired(!opt.single_end), bam(opt.bam), num(opt.num), tc(tc), index(index), mp(mp), id(_id), local_id(_local_id), numreads(0) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];  
   seqs.reserve(bufsize/50);
   newEcs.reserve(1000);
   bv.reserve(1000);
   counts.reserve((int) (tc.counts.size() * 1.25));
   memset(&bc_len[0],0,33);
   memset(&umi_len[0],0,33);

   clear();
}

BUSProcessor::BUSProcessor(BUSProcessor && o) :
  paired(o.paired),
  bam(o.bam),
  num(o.num),
  tc(o.tc),
  index(o.index),
  mp(o.mp),
  id(o.id),
  local_id(o.local_id),
  bufsize(o.bufsize),
  numreads(o.numreads),
  seqs(std::move(o.seqs)), 
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)),
  newEcs(std::move(o.newEcs)),
  flens(std::move(o.flens)),
  bias5(std::move(o.bias5)),
  bv(std::move(o.bv)),  
  counts(std::move(o.counts)) {
    memcpy(&bc_len[0], &o.bc_len[0], 33);
    memcpy(&umi_len[0], &o.umi_len[0], 33);
    buffer = o.buffer;
    o.buffer = nullptr;
    o.bufsize = 0;
}

BUSProcessor::~BUSProcessor() {
  if (buffer != nullptr) {
      delete[] buffer;
      buffer = nullptr;
  }
}

void BUSProcessor::operator()() {
  while (true) {
    int readbatch_id;
    // grab the reader lock
    {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        std::vector<std::string> umis;
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, umis, readbatch_id, mp.opt.pseudobam);
      }
      // release the reader lock
    }
    // do the same for BUS ?!?
    
    pseudobatch.aln.clear();
    pseudobatch.batch_id = readbatch_id;
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    std::vector<std::pair<int, std::string>> ec_umi;
    std::vector<std::pair<std::vector<int>, std::string>> new_ec_umi;
    mp.update(counts, newEcs, ec_umi, new_ec_umi, seqs.size() / mp.opt.busOptions.nfiles , flens, bias5, pseudobatch, bv, newB, &bc_len[0], &umi_len[0], id, local_id);
    clear();
  }
}

void BUSProcessor::processBuffer() {
  // set up thread variables  
  std::vector<std::pair<KmerEntry,int>> v,v2;
  std::vector<int> vtmp;
  std::vector<int> u;
  
  u.reserve(1000);
  v.reserve(1000);
  v2.reserve(1000);
  vtmp.reserve(1000);

  memset(&bc_len[0], 0, 33);
  memset(&umi_len[0], 0, 33);

  const BUSOptions& busopt = mp.opt.busOptions;
  
  bool findFragmentLength = mp.tlencount < 10000 && busopt.paired && !mp.opt.tagsequence.empty();
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

 
  char buffer[100];
  memset(buffer, 0, 100);
  char *umi = &(buffer[50]);
  char *bc  = &(buffer[0]);
  //char *seqbuffer[1000];
  std::string seqbuffer;
  seqbuffer.reserve(1000);
  
  // actually process the sequence
  
  int incf, jmax;
  if (bam) {
    incf = 1;
    jmax = 2;

  } else {
    incf = busopt.nfiles-1;
    jmax = busopt.nfiles;
  }

  std::vector<const char*> s(jmax, nullptr);
  std::vector<int> l(jmax,0);

  bool singleSeq = busopt.seq.size() ==1 ;
  const BUSOptionSubstr seqopt = busopt.seq.front();
  bool check_tag_sequence = !mp.opt.tagsequence.empty();
  int taglen = 0;
  uint64_t tag_binary = 0;
  int hamming_dist = 1;
  if (check_tag_sequence) {
    taglen = mp.opt.tagsequence.length();
    uint32_t f = 0;
    tag_binary = stringToBinary(mp.opt.tagsequence, f);
  }


  //int incf = (bam) ? 1 : busopt.nfiles-1;
  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax /*(bam) ? 2 : busopt.nfiles*/; j++) {
      s[j] = seqs[i+j].first;
      l[j] = seqs[i+j].second;      
    }
    i += incf;
    
    // find where the sequence is
    const char *seq = nullptr;
    const char *seq2 = nullptr;
    
    bool ignore_umi = false;
    
    // copy the umi
    int umilen = (busopt.umi.start == busopt.umi.stop) ? l[busopt.umi.fileno] - busopt.umi.start : busopt.umi.stop - busopt.umi.start;
    if (l[busopt.umi.fileno] < busopt.umi.start + umilen) {
      continue; // too short
    }
    if (check_tag_sequence) {
      umilen += taglen; // Expand to include tag sequence
      memcpy(umi, s[busopt.umi.fileno] + busopt.umi.start - taglen, umilen);
    } else {
      memcpy(umi, s[busopt.umi.fileno] + busopt.umi.start, umilen);
    }
    umi[umilen] = 0;
    uint64_t umi_binary = -1;
    if (umilen >= 0 && umilen <= 32) {
      if (check_tag_sequence) {
        uint32_t f = 0;
        umi_binary = stringToBinary(umi, umilen, f);
        if (hamming(tag_binary, umi_binary >> 2*(umilen-taglen), taglen) <= hamming_dist) { // if tag present
          umi_binary = umi_binary & ~(~0ULL << (2*(umilen - taglen)));
          umilen = umilen - taglen;
          umi_len[umilen]++;
        } else {
          ignore_umi = true; // tag not present, it's not a UMI-read
          umi_binary = -1;
        }
      } else {
        umi_len[umilen]++;
      }
    }
    
    size_t seqlen = 0;
    size_t seqlen2 = 0;
    if (singleSeq) {
      int seqstart = (ignore_umi ? busopt.umi.start - taglen : seqopt.start);
      seq = s[seqopt.fileno] + seqstart;
      seqlen = (seqopt.stop == 0) ? l[seqopt.fileno]-seqstart : seqopt.stop - seqstart;
    } else if (busopt.paired) {
      const auto &sopt1 = busopt.seq[0];
      const auto &sopt2 = busopt.seq[1];
      int seqstart1 = (ignore_umi && busopt.umi.fileno == sopt1.fileno ? busopt.umi.start - taglen : sopt1.start);
      int seqstart2 = (ignore_umi && busopt.umi.fileno == sopt2.fileno ? busopt.umi.start - taglen : sopt2.start);
      int cplen1 = (sopt1.stop == 0) ? l[sopt1.fileno]-seqstart1 : sopt1.stop - seqstart1;
      int cplen2 = (sopt2.stop == 0) ? l[sopt2.fileno]-seqstart2 : sopt2.stop - seqstart2;
      seq = s[sopt1.fileno] + seqstart1;
      seq2 = s[sopt2.fileno] + seqstart2;
      seqlen = cplen1;
      seqlen2 = cplen2;
    } else {
      seqbuffer.clear();
      for (int j = 0; j < busopt.seq.size(); j++) {
        const auto &sopt = busopt.seq[j];
        int seqstart = (ignore_umi && busopt.umi.fileno == sopt.fileno ? busopt.umi.start - taglen : sopt.start);
        int cplen = (sopt.stop == 0) ? l[sopt.fileno]-seqstart : sopt.stop - seqstart;
        seqbuffer.append(s[sopt.fileno] + seqstart, cplen);
        seqbuffer.push_back('N');
      }
      seqbuffer.push_back(0);
      seq = seqbuffer.c_str();
      seqlen = seqbuffer.size();
    }

    
//    auto &bcc = busopt.bc[0];
    int blen = 0;
    bool bad_bc = false;
    for (auto &bcc : busopt.bc) {
      int bclen = (bcc.start == bcc.stop) ? l[bcc.fileno] - bcc.start : bcc.stop - bcc.start;
      if (l[bcc.fileno] < bcc.start + bclen) {
        bad_bc = true;
        break;
      }
      memcpy(bc+blen, s[bcc.fileno] + bcc.start, bclen);
      blen += bclen;
    }
    if (bad_bc) {
      continue;
    }
    bc[blen] = 0;

    if (blen >= 0 && blen <= 32) {
      bc_len[blen]++;
    }
    /* debugging
    std::cout << "seq " << seq << std::endl;
    std::cout << "bc  " << bc << std::endl;
    std::cout << "umi " << umi << std::endl << std::endl;
    */

    numreads++;
    v.clear();
    u.clear();

    // process 2nd read
    index.match(seq,seqlen, v);
    if (busopt.paired) {
      v2.clear();
      index.match(seq2,seqlen2, v2);
    }

    // collect the target information
    int ec = -1;
    int r = tc.intersectKmers(v, v2, false,u);
    if (!u.empty()) {      
      ec = tc.findEC(u);
    }

    if (!ignore_umi && mp.opt.tag_strand_specific && !u.empty()) {
      int p = -1;
      Kmer km;
      KmerEntry val;
      if (!v.empty()) {
        vtmp.clear();
        bool firstStrand = (mp.opt.tag_strand == ProgramOptions::StrandType::FR); // FR have first read mapping forward
        p = findFirstMappingKmer(v,val);
        km = Kmer((seq+p));
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
        bool secondStrand = (mp.opt.tag_strand == ProgramOptions::StrandType::RF);
        p = findFirstMappingKmer(v2,val);
        km = Kmer((seq2+p));
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
      ec = -1;
      if (!u.empty()) { // update the ec
        ec = tc.findEC(u);
      }
    }

    // find the ec
    if (!u.empty()) {
      BUSData b;      
      uint32_t f = 0;
      b.flags = 0;
      b.barcode = stringToBinary(bc, blen, f);
      b.flags |= f;
      b.UMI = check_tag_sequence ? umi_binary : stringToBinary(umi, umilen, f);
      b.flags |= (f) << 8;
      b.count = 1;
      //std::cout << std::string(s1,10)  << "\t" << b.barcode << "\t" << std::string(s1+10,16) << "\t" << b.UMI << "\n";
      if (num) {
        b.flags = (uint32_t) flags[i / jmax];
      }
      
      //ec = tc.findEC(u);

      if (busopt.paired && ignore_umi) {
        if (findFragmentLength && flengoal > 0 && 0 <= ec && ec < index.num_trans && !v.empty() && !v2.empty()) {
          // try to map the reads
          int tl = index.mapPair(seq, seqlen, seq2, seqlen2, ec);
          if (0 < tl && tl < flens.size()) {
            flens[tl]++;
            flengoal--;
          }
        }
      }
      
      // count the pseudoalignment
      if (ec == -1 || ec >= counts.size()) {
        // something we haven't seen before
        newEcs.push_back(u);
        newB.push_back({b,u});
      } else {
        // add to count vector
        ++counts[ec];
        // push back BUS record
        b.ec = ec;
        bv.push_back(b);
      }      
    }

    if (mp.opt.pseudobam) {
      PseudoAlignmentInfo info;
      info.id = i / jmax;
      info.paired = busopt.paired;
      uint32_t f = 0;
      info.barcode = stringToBinary(bc, blen, f);
      info.UMI = check_tag_sequence ? umi_binary : stringToBinary(umi, umilen, f);
      if (!u.empty()) {
        info.r1empty = v.empty();
        KmerEntry val;
        info.k1pos = findFirstMappingKmer(v,val);
        info.k2pos = -1;
        if (busopt.paired) {
          info.r2empty = v2.empty();
          info.k2pos = findFirstMappingKmer(v2,val);
        }
                
        if (ec != -1) {
          info.ec_id = ec;
        } else {
          info.ec_id = -1;
          info.u = u; // copy
        }
      }
      pseudobatch.aln.push_back(std::move(info));
    }

    if (mp.opt.verbose && numreads > 0 && numreads % 1000000 == 0 ) {   
        int nmap = mp.nummapped;
        for (int i = 0; i < counts.size(); i++) {
          nmap += counts[i];
        }
        nmap += newEcs.size();

        std::cerr << '\r' << (numreads/1000000) << "M reads processed (" 
          << std::fixed << std::setw( 3 ) << std::setprecision( 1 ) << ((100.0*nmap)/double(numreads))
          << "% pseudoaligned)"; std::cerr.flush();
      }
  }
}

void BUSProcessor::clear() {
  numreads=0;
  memset(buffer,0,bufsize);
  newEcs.clear();
  counts.clear();  
  counts.resize(tc.counts.size(),0);
  bv.clear();
  newB.clear();
}


AlnProcessor::AlnProcessor(const KmerIndex& index, const ProgramOptions& opt, MasterProcessor& mp, const EMAlgorithm& em, const Transcriptome &model, bool useEM, int _id) :
 paired(!opt.single_end), index(index), mp(mp), em(em), model(model), useEM(useEM), id(_id) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];
   bambufsize = 1<<20;
   bambuffer = new char[bambufsize]; // refactor this?
   
   if (opt.batch_mode) {
     /* need to check this later */
     assert(id != -1);
     batchSR.files = opt.batch_files[id];
     batchSR.nfiles = opt.batch_files[id].size();
     batchSR.reserveNfiles(opt.batch_files[id].size());
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
  em(o.em),
  mp(o.mp),
  model(o.model),
  id(o.id),
  useEM(o.useEM),
  bufsize(o.bufsize),
  bambufsize(o.bambufsize),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)),
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
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (batchSR.empty()) {
        return;
      } else {
        batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, flags, umis, readbatch_id, true );
        readPseudoAlignmentBatch(mp.pseudobatchf_in, pseudobatch);
        assert(pseudobatch.batch_id == readbatch_id);
        assert(pseudobatch.aln.size() == ((paired) ? seqs.size()/2 : seqs.size())); // sanity checks
      }
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences        
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, umis, readbatch_id, true);
        readPseudoAlignmentBatch(mp.pseudobatchf_in, pseudobatch);
        assert(pseudobatch.batch_id == readbatch_id);
        if (mp.opt.bus_mode) {
          assert(pseudobatch.aln.size() == seqs.size()/mp.opt.busOptions.nfiles); // sanity checks
        } else {
          assert(pseudobatch.aln.size() == ((paired) ? seqs.size()/2 : seqs.size())); // sanity checks
        }
      }
      // release the reader lock
    }
    // process our sequences
    if (mp.opt.genomebam) {
      processBufferGenome();

      // append files in order. clean up temporary stuff
    } else {
      processBufferTrans();
    }


  }
}


void AlnProcessor::processBufferTrans() {
  /* something simple where we can construct the bam records */
  std::vector<bam1_t> bv;


  int n = pseudobatch.aln.size();
  std::vector<int> u;
  std::vector<std::pair<int,double>> ua;
  u.reserve(1000);
  ua.reserve(1000);

  int trans_auxlen = default_trans_auxlen;
  if (!useEM) {
    trans_auxlen -= 7;
  }
  

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
    fillBamRecord(b1, nullptr, seqs[si1].first, names[si1].first,  quals[si1].first, seqs[si1].second, names[si1].second, pi.r1empty, trans_auxlen);
    if (paired) {
      fillBamRecord(b2, nullptr, seqs[si2].first, names[si2].first,  quals[si2].first, seqs[si2].second, names[si2].second, pi.r2empty, trans_auxlen);
    }

    if (pi.r1empty && pi.r2empty) {      
      bv.push_back(b1);
      if (paired) {      
        bv.push_back(b2);
      }
    } else {
      u.clear();
      ua.clear();
      int ec = pi.ec_id;
      if (ec != -1) {
        u = index.ecmap[ec]; // copy, but meh
      } else {
        u = pi.u;
        auto it = index.ecmapinv.find(u);
        if (it != index.ecmapinv.end()) {
          ec = it->second;
        }

        if (ec == -1) {
          assert(false && "Problem with ecmapinv");
        }
      }

      
      // modify u and compute norm
      int32_t nmap = u.size();
      int bestProbTr = -1;
      double bestProb = 0.0;

      if (useEM) {

        double denom = 0.0;
        const auto& wv = em.weight_map_[ec];
        
        for (int i = 0; i < nmap; ++i) {
          denom += em.alpha_[u[i]] * wv[i];
        }
        
        if (denom < TOLERANCE) {
          u.clear();
          ua.clear();
        } else {
          // compute the update step
          for (int i = 0; i < nmap; ++i) {
            if (em.alpha_[u[i]] > 0.0) {
              double prob = em.alpha_[u[i]] * wv[i] / denom;
              ua.push_back({u[i],prob});
              if (bestProb < prob) {
                bestProb = prob;
                bestProbTr = u[i];
              }
            }
          }
        }
      } else {
        for (int i = 0; i < nmap; i++) {
          ua.push_back({u[i], 0.0}); // probability is never used
        }
        bestProbTr = u[0];
      }

      nmap = ua.size();

      if (ua.empty()) {
        // shouldn't happen
        bv.push_back(b1);
        if (paired) {
          bv.push_back(b2);
        }        
      } else {
        // set flags
        if (!pi.r1empty) {
          b1.core.flag &= ~BAM_FUNMAP;
          if (paired) {
            b2.core.flag &= ~BAM_FMUNMAP;
          }
        }
        if (paired) {
          if (!pi.r2empty) {
            b1.core.flag &= ~BAM_FMUNMAP;
            b2.core.flag &= ~BAM_FUNMAP;
          }
          if (!pi.r1empty && !pi.r2empty) {
            b1.core.flag |= BAM_FPROPER_PAIR;
            b2.core.flag |= BAM_FPROPER_PAIR;
          }              
        }
        
        // set aux
        float zero = 0.0;
        assert(b1.l_data + trans_auxlen <= b1.m_data);      

        b1.data[b1.l_data] = 'N';
        b1.data[b1.l_data+1] = 'H';
        b1.data[b1.l_data+2] = 'i';
        memcpy(b1.data + b1.l_data + 3, &nmap, 4);
        b1.l_data += 7;

        if (useEM) {
          b1.data[b1.l_data] = 'Z';
          b1.data[b1.l_data+1] = 'W';
          b1.data[b1.l_data+2] = 'f';
          memcpy(b1.data + b1.l_data + 3, &zero, 4);
          b1.l_data += 7;
        }

        if (paired) {
          assert(b2.l_data + trans_auxlen <= b2.m_data);
          
          b2.data[b2.l_data] = 'N';
          b2.data[b2.l_data+1] = 'H';
          b2.data[b2.l_data+2] = 'i';
          memcpy(b2.data + b2.l_data + 3, &nmap, 4);
          b2.l_data += 7;

          if (useEM) {
            b2.data[b2.l_data] = 'Z';
            b2.data[b2.l_data+1] = 'W';
            b2.data[b2.l_data+2] = 'f';
            memcpy(b2.data + b2.l_data + 3, &zero, 4);
            b2.l_data += 7;
          }
        }

        auto strandednessInfo = [&](Kmer km, KmerEntry& val, const std::vector<std::pair<int,double>> &ua) -> std::pair<bool,bool> {
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
                  for (const auto &y : ua) {
                    if (y.first == x.trid) {
                      return {false,reptrue}; 
                    }
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
          strInfo1 = strandednessInfo(km1, val1, ua);

        }
        if (paired && !pi.r2empty) {
          km2 = Kmer(seqs[si2].first + pi.k2pos);
          strInfo2 = strandednessInfo(km2, val2, ua);
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


        bool bestTr = false;
        for (auto tp : ua) {        
          auto t = tp.first;
          float prob = (float) tp.second; // explicit cast
          bestTr = (t == bestProbTr);
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

          if (!bestTr) {
            b1c.core.flag |= BAM_FSECONDARY;
            if (paired) {
              b2c.core.flag |= BAM_FSECONDARY;
            }
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
            b1c.core.bin = hts_reg2bin(b1c.core.pos, b1c.core.pos + seqs[si1].second-1, 14, 5);
            b1c.core.qual = 255;
            if (useEM) {
              memcpy(b1c.data + b1c.l_data - 4, &prob, 4); // set ZW tag           
            }
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

              b2c.core.bin = hts_reg2bin(b2c.core.pos, b2c.core.pos + seqs[si2].second, 14, 5);
              b2c.core.qual = 255;
              if (useEM) {
                memcpy(b2c.data + b2c.l_data - 4, &prob, 4);
              }

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
              fixCigarStringTrans(b1c, rlen1, softclip1, overhang1);
            }
          }
          if (paired && !pi.r2empty) {
            if (softclip2 > 0 || overhang2 > 0) {
              fixCigarStringTrans(b2c, rlen2, softclip2, overhang2);
            }
          }

          if (!pi.r1empty || bestTr) {
            bv.push_back(b1c);
          }
          if (paired && (!pi.r2empty || bestTr)) {
            bv.push_back(b2c);
          }
        } 

        delete[] b1.data;
        if (paired) {
          delete[] b2.data;
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


void AlnProcessor::processBufferGenome() {
  /* something simple where we can construct the bam records */
  std::vector<bam1_t> bv;
  int idnum = 0;

  int n = pseudobatch.aln.size();
  std::vector<int> u;
  std::vector<std::pair<int,double>> ua;
  u.reserve(1000);
  ua.reserve(1000);

  int bclen = 0;
  int umilen = 0;
  
  int genome_auxlen = default_genome_auxlen;
  if (!useEM) {
    genome_auxlen -= 7;
  }
  if (mp.opt.bus_mode) {
    //todo replace with what is written to busfile
    bclen = mp.opt.busOptions.getBCLength();
    if (bclen == 0) {
      bclen = 32;
    }
    umilen = mp.opt.busOptions.getUMILength();
    if (umilen == 0) {
      umilen = 32;
    }

    genome_auxlen += (bclen + umilen) + 8;
  }
  std::string bc, umi;
  

  Kmer km1,km2;
  KmerEntry val1, val2;
  if  (mp.opt.bus_mode) {
    paired = false;
    if (mp.opt.busOptions.paired) {
      paired = true;
    }
  }

  for (int i = 0; i < n; i++) {
    bam1_t b1,b2, b1c, b2c;
    int si1 = (paired) ? 2*i : i;
    int si2 = (paired) ? 2*i +1 : -1;
    // For BUS file processing with paired-end reads:
    if (mp.opt.bus_mode) {
      const BUSOptions& busopt = mp.opt.busOptions;
      si1 = i*busopt.nfiles + busopt.seq[0].fileno;
      if (paired) {
        si2 = i*busopt.nfiles + busopt.seq[1].fileno;
      }
    }
    int rlen1 = seqs[si1].second;
    int rlen2;
    int seq1_offset = 0;
    int seq2_offset = 0;
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
    
    if (!mp.opt.bus_mode) {
      // fill in the data info
      fillBamRecord(b1, nullptr, seqs[si1].first, names[si1].first,  quals[si1].first, seqs[si1].second, names[si1].second, pi.r1empty, genome_auxlen);
      //b1.id = idnum++;
      if (paired) {
        fillBamRecord(b2, nullptr, seqs[si2].first, names[si2].first,  quals[si2].first, seqs[si2].second, names[si2].second, pi.r2empty, genome_auxlen);
        // b2.id = idnum++;
      }
    }
    else {
      bool tag_present = !mp.opt.tagsequence.empty() && pi.UMI != -1; // tag+UMI present in the sequence
      bool umi_first_file = mp.opt.busOptions.umi.fileno == mp.opt.busOptions.seq[0].fileno; // UMI-containing file is the first file (file number 0)
      int tagseqstart;
      if (tag_present) {
        tagseqstart = umi_first_file ? mp.opt.busOptions.seq[0].start : mp.opt.busOptions.seq[1].start; // Position where the actual sequence begins when tag+UMI present
      }

      // fill in the data info
      if (tag_present && umi_first_file) {
        rlen1 = seqs[si1].second-tagseqstart;
        seq1_offset = tagseqstart;
        fillBamRecord(b1, nullptr, seqs[si1].first+tagseqstart, names[si1].first,  quals[si1].first+tagseqstart, seqs[si1].second-tagseqstart, names[si1].second, pi.r1empty, genome_auxlen);
      } else {
        fillBamRecord(b1, nullptr, seqs[si1].first, names[si1].first,  quals[si1].first, seqs[si1].second, names[si1].second, pi.r1empty, genome_auxlen);
      }
      if (paired) {
        if (tag_present && !umi_first_file) {
          rlen2 = seqs[si2].second-tagseqstart;
          seq2_offset = tagseqstart;
          fillBamRecord(b2, nullptr, seqs[si2].first+tagseqstart, names[si2].first,  quals[si2].first+tagseqstart, seqs[si2].second-tagseqstart, names[si2].second, pi.r2empty, genome_auxlen);
        } else {
          fillBamRecord(b2, nullptr, seqs[si2].first, names[si2].first,  quals[si2].first, seqs[si2].second, names[si2].second, pi.r2empty, genome_auxlen);
        }
      }

      b1.data[b1.l_data] = 'C';
      b1.data[b1.l_data+1] = 'R';
      b1.data[b1.l_data+2] = 'Z';
      bc = binaryToString(pi.barcode, bclen);
      memcpy(b1.data + b1.l_data + 3, bc.c_str(), bclen+1);
      b1.l_data += bclen + 4;
      
      b1.data[b1.l_data] = 'U';
      b1.data[b1.l_data+1] = 'R';
      b1.data[b1.l_data+2] = 'Z';
      if (!mp.opt.tagsequence.empty() && pi.UMI == -1) { // Non-UMI read
        umi = std::string(umilen, 'N');
      } else {
        umi = binaryToString(pi.UMI, umilen);
      }
      memcpy(b1.data + b1.l_data + 3, umi.c_str(), umilen+2);
      b1.l_data += umilen + 4;
      if (paired) {
        b2.data[b2.l_data] = 'C';
        b2.data[b2.l_data+1] = 'R';
        b2.data[b2.l_data+2] = 'Z';
        memcpy(b2.data + b2.l_data + 3, bc.c_str(), bclen+1);
        b2.l_data += bclen + 4;
        b2.data[b2.l_data] = 'U';
        b2.data[b2.l_data+1] = 'R';
        b2.data[b2.l_data+2] = 'Z';
        memcpy(b2.data + b2.l_data + 3, umi.c_str(), umilen+2);
        b2.l_data += umilen + 4;
      }
    }

    if (pi.r1empty && pi.r2empty) {
      bv.push_back(b1);
      if (paired) {
        bv.push_back(b2);
      }
    } else {
      u.clear();
      ua.clear();
      int ec = pi.ec_id;
      if (ec != -1) {
        u = index.ecmap[ec]; // copy, but meh
      } else {
        u = pi.u;
        auto it = index.ecmapinv.find(u);
        if (it != index.ecmapinv.end()) {
          ec = it->second;
        }

        if (ec == -1) {
          assert(false && "Problem with ecmapinv");
        }
      }

      
      // modify u and compute norm
      int32_t nmap = u.size();
      if (useEM) {
        double denom = 0.0;
        const auto& wv = em.weight_map_[ec];
        
        for (int i = 0; i < nmap; ++i) {
          denom += em.alpha_[u[i]] * wv[i];
        }
        
        if (denom < TOLERANCE) {
          u.clear();
          ua.clear();
        } else {
           // compute the update step
          for (int i = 0; i < nmap; ++i) {
            if (em.alpha_[u[i]] > 0.0) {
              double prob = em.alpha_[u[i]] * wv[i] / denom;
              ua.push_back({u[i],prob});
            }
          }
        }
      } else {
        for (int i = 0; i < nmap; i++) {
          ua.push_back({u[i], 0.0}); // never used
        }
      }

      nmap = ua.size();

      if (ua.empty()) {
        // shouldn't happen
        bv.push_back(b1);
        if (paired) {
          bv.push_back(b2);
        }        
      } else {
        // set flags
        if (!pi.r1empty) {
          b1.core.flag &= ~BAM_FUNMAP;
          if (paired) {
            b2.core.flag &= ~BAM_FMUNMAP;
          }
        }
        if (paired) {
          if (!pi.r2empty) {
            b1.core.flag &= ~BAM_FMUNMAP;
            b2.core.flag &= ~BAM_FUNMAP;
          }
          if (!pi.r1empty && !pi.r2empty) {
            b1.core.flag |= BAM_FPROPER_PAIR;
            b2.core.flag |= BAM_FPROPER_PAIR;
          }              
        }
          
        // set aux
        float zero = 0.0;
        //assert(b1.l_data + genome_auxlen <= b1.m_data);
        
        if (useEM) {
          b1.data[b1.l_data] = 'Z';
          b1.data[b1.l_data+1] = 'W';
          b1.data[b1.l_data+2] = 'f';
          memcpy(b1.data + b1.l_data + 3, &zero, 4);
          b1.l_data += 7;

          if (paired) {
            //assert(b2.l_data + genome_auxlen <= b2.m_data);
            b2.data[b2.l_data] = 'Z';
            b2.data[b2.l_data+1] = 'W';
            b2.data[b2.l_data+2] = 'f';
            memcpy(b2.data + b2.l_data + 3, &zero, 4);
            b2.l_data += 7;
          }
        }


        // everything maps to the same strand on all transcriptomes
        auto strandednessInfo = [&](Kmer km, KmerEntry& val, const std::vector<std::pair<int,double>> &ua) -> std::pair<bool,bool> {          
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
              bool chrsense = model.transcripts[c.transcripts[0].trid].strand == c.transcripts[0].sense;
              for (const auto & x : c.transcripts) {
                if ((model.transcripts[x.trid].strand == x.sense) != chrsense) {
                  for (const auto &y : ua) {
                    if (y.first == x.trid) {
                      return {false,reptrue}; 
                    }
                  }
                }
              }
              return {true, chrsense == (reptrue == val.isFw())};
            }
          }
        };
        std::pair<bool,bool> strInfo1 = {true,true}, strInfo2 = {true,true};
        
        if (!pi.r1empty) {
          km1 = Kmer(seqs[si1].first + seq1_offset + pi.k1pos);
          strInfo1 = strandednessInfo(km1, val1, ua);
        }
        if (paired && !pi.r2empty) {
          km2 = Kmer(seqs[si2].first + seq2_offset + pi.k2pos);
          strInfo2 = strandednessInfo(km2, val2, ua);
        }
        
        
        // first read is all on reverse strand
        if (strInfo1.first && !strInfo1.second) {
          reverseComplementSeqInData(b1);          
        }
        
        // we have second read and it is all on reverse strand
        if (paired && strInfo2.first && !strInfo2.second) {
          reverseComplementSeqInData(b2);
        }

        std::unordered_map<std::pair<TranscriptAlignment, TranscriptAlignment>, double> alnmap;
        TranscriptAlignment tra1, tra2;
        for (auto tp : ua) {
          auto t = tp.first;
          double prob = tp.second;
          std::pair<int,bool> pos1, pos2;

          if (!pi.r1empty) {
            pos1 = index.findPosition(t, km1, val1, pi.k1pos);
            int trpos = (pos1.second) ? pos1.first-1 : pos1.first - rlen1;
            if (!model.translateTrPosition(t, trpos, rlen1, pos1.second, tra1)) {
              continue;
            }
          }
          
          if (paired) {
            if (!pi.r2empty) {
              pos2 = index.findPosition(t, km2, val2, pi.k2pos);
              int trpos = (pos2.second) ? pos2.first-1 : pos2.first - rlen2;
              if (!model.translateTrPosition(t, trpos, rlen2, pos2.second, tra2)) {
                continue;
              }
            }
          }

          alnmap[{tra1,tra2}] += prob;

        }

        if (alnmap.size() == 0) {
          bv.push_back(b1);
          if (paired) {
            bv.push_back(b2);
          }
          continue;
        }
    
        double bestprob = 0.0;
        std::pair<TranscriptAlignment, TranscriptAlignment> bestTra;
        if (alnmap.size() == 1) {
          // common case
          bestTra.first = std::move(tra1);
          bestTra.second = std::move(tra2);
          bestprob = 1.0;
        } else {
          for (auto &x : alnmap) {
            if (x.second > bestprob) {
              bestprob = x.second;
              bestTra = x.first; // copy
            }
          }
          if (useEM) {
            bestTra = alnmap.begin()->first; // arbitrary
          }
        }

        for (auto &x : alnmap) {
          auto &tra = x.first;
          float prob = (float) x.second; // explicit cast
          
          bool bestTr = (bestprob == 1.0) || (tra == bestTra);

          b1c = b1;
          // b1c.id = idnum++;
          b1c.data = new uint8_t[b1c.m_data];
          memcpy(b1c.data, b1.data, b1c.m_data*sizeof(uint8_t));          
          if (!strInfo1.first && !tra.first.strand) {
            reverseComplementSeqInData(b1c);
          }

          if (paired) {
            b2c = b2;
            // b2c.id = idnum++;
            b2c.data = new uint8_t[b2c.m_data];
            memcpy(b2c.data, b2.data, b2c.m_data*sizeof(uint8_t));
            if (!strInfo2.first && !tra.second.strand) {
              reverseComplementSeqInData(b2c);
            }
          }

          if (paired) {
            // todo handle non mapping reads
            if (!pi.r1empty && !tra.first.strand) {
              b1c.core.flag |= BAM_FREVERSE;
              b2c.core.flag |= BAM_FMREVERSE;
            }
            if (!pi.r2empty && !tra.second.strand) {
              b1c.core.flag |= BAM_FMREVERSE; // 
              b2c.core.flag |= BAM_FREVERSE;
            }            
          } else {
            if (!pi.r1empty && !tra.first.strand) {
              b1c.core.flag |= BAM_FREVERSE;
            }
          }

          if (!bestTr) {
            b1c.core.flag |= BAM_FSECONDARY;
            if (paired) {
              b2c.core.flag |= BAM_FSECONDARY;
            }
          }
        
          b1c.core.tid = tra.first.chr;
          if (!pi.r1empty) {
            b1c.core.pos = tra.first.chrpos;
            b1c.core.bin = hts_reg2bin(b1c.core.pos, b1c.core.pos + seqs[si1].second-seq1_offset-1, 14, 5);
            b1c.core.qual = 255;
            if (useEM) {
              memcpy(b1c.data + b1c.l_data - 4, &prob, 4); // set ZW tag
            }
          }

          if (paired) {
            b2c.core.tid = tra.second.chr;
            if (!pi.r2empty) {
              b2c.core.pos = tra.second.chrpos;
              b2c.core.bin = hts_reg2bin(b2c.core.pos, b2c.core.pos + seqs[si2].second-seq2_offset, 14, 5);
              b2c.core.qual = 255;
              if (useEM) {
                memcpy(b2c.data + b2c.l_data - 4, &prob, 4);
              }

              if (pi.r1empty) {
                b1c.core.tid = b2c.core.tid;
                b1c.core.pos = b2c.core.pos;
                b1c.core.bin = b2c.core.bin;
                b1c.core.qual = 0;
              }
            } else {
              b2c.core.tid = b1c.core.tid;
              b2c.core.pos = b1c.core.pos;
              b2c.core.bin = b2c.core.bin;
              b2c.core.qual = 0;
            }
            b1c.core.mtid = b2c.core.tid;
            b1c.core.mpos = b2c.core.pos;
            b2c.core.mtid = b1c.core.tid;
            b2c.core.mpos = b1c.core.pos;
  
          }

          
          if (!pi.r1empty && !pi.r2empty) {
            int tlen = bam_endpos(&b2c) - b1c.core.pos;
            b1c.core.isize = tlen;
            b2c.core.isize = -tlen;
          }          
        

          if (!pi.r1empty && tra.first.cigar.size() > 1) {
            fixCigarStringGenome(b1c, tra.first);
          }
          if (paired && !pi.r2empty && tra.second.cigar.size() > 1) {
            fixCigarStringGenome(b2c, tra.second);
          }

        
          if (!pi.r1empty || bestTr) {
            bv.push_back(b1c);
          }
          if (paired && (!pi.r2empty || bestTr)) {
            bv.push_back(b2c);
          }
        }  
      
        delete[] b1.data;
        // std::cerr << "deleting" << b1.id << std::endl;
        if (paired) {
          delete[] b2.data;
          // std::cerr << "deleting" << b2.id << std::endl;
        }
      }
    }    
  }


  // partition the vectors. 
  int k = mp.numSortFiles;
  int logk = 1;
  while (k > (1<<logk)) {
    logk++;
  }
  const auto& bp = mp.breakpoints;
  std::vector<std::vector<bam1_t>> bvv(k);
  
  for (auto &b : bv) {
    uint64_t pos = ((uint64_t) b.core.tid) << 32 | ((uint64_t) b.core.pos+1) << 1 | (b.core.flag & BAM_FREVERSE) >>4;
    int i=0,j=k;
    
    if (b.core.tid == -1) {
      i = k-1;
    } else {
      // invariant bp[i-1].pos <= b.pos < bp[j].pos
      for (int r = 0; r <= logk; r++) {
        int m = (i+j)/2;
        /*assert(i <= j);
        assert(i == 0 || bp[i-1] <= pos);
        assert(j == k || bp[j] > pos);         
        assert(i <= m);
        assert(m <= j);*/
        if (pos < bp[m]) {
          j = m;
        } else {
          i = m+1;
        }
      }
      //assert(i==j);
    }
    bvv[i].push_back(b);
  }


  // output to disk
  mp.writeSortedPseudobam(bvv);

  // clean up our mess
  for (auto &b : bv) {
    delete[] b.data;
    // std::cerr << "deleting" << b.id  << std::endl;
  }  
  bv.clear();  
}

void fixCigarStringGenome(bam1_t &b, const TranscriptAlignment& tra) {
  int ncig = tra.cigar.size();
  if (ncig == 1) {
    return;
  }

  int extraspace = b.m_data - b.l_data;
  extraspace -= ((int) (ncig - b.core.n_cigar)*sizeof(uint32_t));

  if (extraspace < 0) {
    //allocate more space
    int n = b.m_data - extraspace;
    uint8_t* bf = new uint8_t[n];
    memcpy(bf, b.data, b.m_data);
    delete[] b.data;
    // std::cerr << "temporarily deleting " << b.id  << std::endl;
    b.data = bf;
    b.m_data = n;
  }


  uint8_t *cigafter = b.data + b.core.l_qname + b.core.n_cigar * sizeof(uint32_t);
  int copylen = b.l_data - (b.core.l_qname + b.core.n_cigar * sizeof(uint32_t));
  uint8_t *dst = cigafter+((int) (ncig-b.core.n_cigar))*sizeof(uint32_t);
  memmove(dst, cigafter, copylen);
  uint32_t *cig = (uint32_t*) (b.data + b.core.l_qname);
  b.l_data += ((int)(ncig - b.core.n_cigar))*sizeof(uint32_t);
  b.core.n_cigar = ncig;

  for (int i = 0; i < ncig; i++) {
    *cig = tra.cigar[i];
    ++cig;
  }
}

void fixCigarStringTrans(bam1_t &b, int rlen, int softclip, int overhang) {
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
  uint8_t *dst = cigafter+((int) (ncig-b.core.n_cigar))*sizeof(uint32_t);
  memmove(dst, cigafter, copylen);
  uint32_t *cig = (uint32_t*) (b.data + b.core.l_qname);
  b.l_data += ((int) (ncig - b.core.n_cigar))*sizeof(uint32_t);
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
    int i = slen >> 1;
    int i1 = i >> 1;
    int ca = s[i1] >> ((~(i)&1)<<2) & 0xf;
    s[i1] &= 0xF << ((i&1)<<2);
    s[i1] |= bitrev[ca] << ((~i&1)<<2);
  }

  uint8_t *q = b.data + (b.core.n_cigar << 2) + b.core.l_qname + ((b.core.l_qseq+1)>>1);  
  for (int i = 0; i < (slen>>1); ++i) {
    uint8_t t = q[i];
    q[i] = q[slen-1-i];
    q[slen-1-i] = t;
  }
}

int fillBamRecord(bam1_t &b, uint8_t* buf, const char *seq, const char *name, const char *qual, int slen, int nlen, bool unmapped, int auxlen) {
 
  b.core.l_extranul =  (3 - (nlen % 4));  
  b.core.l_qname = nlen + b.core.l_extranul + 1;
  b.core.l_qseq = slen;

  if (buf == nullptr) {
    // allocate memory for buffer
    int blen = b.core.l_qname + 16 + ((b.core.l_qseq+2)>>1) + b.core.l_qseq + auxlen; // 16 for cigar, auxlen for auxiliary fields
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


/** -- sequence readers -- **/

void SequenceReader::reset() {
  state = false;
  readbatch_id = -1;
}

FastqSequenceReader::~FastqSequenceReader() {
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
  }

  for (auto &s : seq) {
    kseq_destroy(s);
  }

  
  // check if umi stream is open, then close
  if (f_umi && f_umi->is_open()) {
    f_umi->close();
  }
}


bool FastqSequenceReader::empty() {
  return (!state && current_file >= files.size());
}

void FastqSequenceReader::reset() {
  SequenceReader::reset();
   
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
    f = nullptr;
  }

  if (f_umi && f_umi->is_open()) {
    f_umi->close();    
  }

  f_umi->clear();

  for (auto &ll : l) {
    ll = 0;
  }
  for (auto &nll : nl) {
    nll = 0;
  }
  
  current_file = 0;
  for (auto &s : seq) {
    kseq_destroy(s);
    s = nullptr;
  }
}

void FastqSequenceReader::reserveNfiles(int n) {
  fp.resize(nfiles);
  l.resize(nfiles, 0);
  nl.resize(nfiles, 0);
  seq.resize(nfiles, nullptr);
}

// returns true if there is more left to read from the files
bool FastqSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals,
  std::vector<uint32_t>& flags,
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
  flags.clear();
   
  bool usingUMIfiles = !umi_files.empty();
  int umis_read = 0;
  
  int bufpos = 0;
  int pad = nfiles; //(paired) ? 2 : 1;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current files
        for (auto &f : fp) {
          if (f) {
            gzclose(f);
          }
        }
        // close current umi file
        if (usingUMIfiles) {
          // read up the rest of the files          
          f_umi->close();
        }
        
        // open the next one
        for (int i = 0; i < nfiles; i++) {
          fp[i] = gzopen(files[current_file+i].c_str(), "r");
          seq[i] = kseq_init(fp[i]);
          l[i] = kseq_read(seq[i]);
          
        }
        if (usingUMIfiles) {
          // open new umi file
          f_umi->open(umi_files[current_file]);  
          current_file++;        
        }
        current_file+=nfiles;
        state = true; 
      }
    }
    // the file is open and we have read into seq1 and seq2
    bool all_l = true;
    int bufadd = nfiles;
    for (int i = 0; i < nfiles; i++) {
      all_l = all_l && l[i] >= 0;
      bufadd += l[i]; // includes seq
    }
    if (all_l) {      
      // fits into the buffer
      if (full) {
        for (int i = 0; i < nfiles; i++) {
          nl[i] = seq[i]->name.l;
          bufadd += l[i] + nl[i]; // includes name and qual
        }
        bufadd += 2*pad;
      }

      if (bufpos+bufadd< limit) {

        for (int i = 0; i < nfiles; i++) {
          char *pi = buf + bufpos;
          memcpy(pi, seq[i]->seq.s, l[i]+1);
          bufpos += l[i]+1;
          seqs.emplace_back(pi,l[i]);

          if (full) {
            pi = buf + bufpos;
            memcpy(pi, seq[i]->qual.s,l[i]+1);
            bufpos += l[i]+1;
            quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          }
        }

        if (usingUMIfiles) {
          std::stringstream ss;
          std::getline(*f_umi, line);
          ss.str(line);
          ss >> umi;
          umis.emplace_back(std::move(umi));
        }

        numreads++;
        flags.push_back(numreads-1);
      } else {
        return true; // read it next time
      }

      // read for the next one
      for (int i = 0; i < nfiles; i++) {
        l[i] = kseq_read(seq[i]);
      }        
    } else {
      state = false; // haven't opened file yet
    }
  }
}

// move constructor

FastqSequenceReader::FastqSequenceReader(FastqSequenceReader&& o) :
  nfiles(o.nfiles),
  numreads(o.numreads),
  fp(std::move(o.fp)),
  l(std::move(o.l)),
  nl(std::move(o.nl)),
  paired(o.paired),
  files(std::move(o.files)),
  umi_files(std::move(o.umi_files)),
  f_umi(std::move(o.f_umi)),
  current_file(o.current_file),
  seq(std::move(o.seq)) {

  o.fp.resize(nfiles);
  o.l.resize(nfiles, 0);
  o.nl.resize(nfiles, 0);
  o.seq.resize(nfiles, nullptr);
  o.state = false;
}

const std::string BamSequenceReader::seq_enc = "=ACMGRSVTWYHKDBN";

BamSequenceReader::~BamSequenceReader() {
  if (fp) {
    bgzf_close(fp);
  }
  if (head) {
    bam_hdr_destroy(head);
  }
  if (rec) {
    bam_destroy1(rec);
  }
}

bool BamSequenceReader::empty() {
  return !state;
}

void BamSequenceReader::reset() {
  SequenceReader::reset();
}

void BamSequenceReader::reserveNfiles(int n) {
}

// returns true if there is more left to read from the files
bool BamSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals,
  std::vector<uint32_t>& flags,
  std::vector<std::string> &umis, int& read_id,
  bool full) {
  
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
  seqs.clear();
  flags.clear();
  
  int bufpos = 0;
  while (true) {
    if (!state) { // should we open a file
        return false;
    }

    // the file is open and we have read into seq1 and seq2
    if (err >= 0) {
      if (!(rec->core.flag & BAM_FSECONDARY)) { // record is primary alignment
        int bufadd = l_seq + l_bc + l_umi + 2;
        // fits into the buffer
        if (bufpos+bufadd < limit) {
          char *pi;

          pi = buf + bufpos;
          memcpy(pi, bc, l_bc);
          memcpy(pi + l_bc, umi, l_umi + 1);
          seqs.emplace_back(pi, l_bc + l_umi);
          bufpos += l_bc + l_umi + 1;

          pi = buf + bufpos;
          int len = (l_seq + 1) / 2;
          if (l_seq % 2) --len;
          int j = 0;
          for (int i = 0; i < len; ++i, ++eseq) {
            buf[bufpos++] = seq_enc[*eseq >> 4];
            buf[bufpos++] = seq_enc[*eseq & 0x0F];
          }
          if (l_seq % 2) {
            buf[bufpos++] = seq_enc[*eseq >> 4];
          }
          buf[bufpos++] = '\0';
          seqs.emplace_back(pi, l_seq);

        } else {
          return true; // read it next time
        }
      }
      
      // read for the next one
      err = bam_read1(fp, rec);
      if (err < 0) {
        return true;
      }
      eseq = bam_get_seq(rec);
      l_seq = rec->core.l_qseq;

      bc = bam_aux2Z(bam_aux_get(rec, "CR"));
      l_bc = 0;
      for (char *pbc = bc; *pbc != '\0'; ++pbc) {
        ++l_bc;
      }

      umi = bam_aux2Z(bam_aux_get(rec, "UR"));
      l_umi = 0;
      for (char *pumi = umi; *pumi != '\0'; ++pumi) {
        ++l_umi;
      }

    } else {
      state = false; // haven't opened file yet
    }
  }
}

