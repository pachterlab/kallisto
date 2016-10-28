#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>

#include "common.h"
#include "KmerIndex.h"
#include "EMAlgorithm.h"
#include "H5Writer.h"


void LoadEcFromFile(KmerIndex &index, const std::string &fn, const ProgramOptions& opt) {
  int N = index.num_trans;
  int num_old_ecs = index.ecmap.size();

  // open the file, parse and fill the batch_files values
  std::vector<std::vector<int>> new_ecmap;
  std::ifstream in(fn);
  std::string line;
  int ec, old_ec = -1;
  std::string trx;
  while (std::getline(in,line)) {
    if (line.size() == 0) {
      continue;
    }
    std::stringstream ss(line);
    ss >> ec >> trx;
    if (ec != old_ec +1) {
      std::cerr << "Error: file " << fn << " is not in correct format" << std::endl;
      exit(1);
    }
    old_ec = ec;
    
    std::stringstream trxs(trx);
    std::string tr;
    std::vector<int> v;
    while (getline(trxs, tr, ',')) {
      v.push_back(std::atoi(tr.c_str()));
    }
    new_ecmap.push_back(std::move(v));
  }


  bool badEC = false;
  for (int i = 0; i < new_ecmap.size(); i++) {
    if (i < num_old_ecs) {
      const auto &u = index.ecmap[i];
      const auto &v = new_ecmap[i];
      if (u != v) {
        badEC = true; assert(false);
        break;
      }
    } else {
      const auto &v = new_ecmap[i];

      auto it = index.ecmapinv.find(v);
      if (it != index.ecmapinv.end()) {
        badEC = true; assert(false);
        break;
      }

      for (int x : v) {
        if (x < 0 || x >= N) {
          badEC = true;assert(false);
          break;
        }
      }
      if (badEC) {
        break;
      }
      // ok add to ecmap and ecmapinv
      index.ecmap.push_back(v); // copy vector 
      index.ecmapinv.insert({v, i});
    }

  }
  if (badEC) {
    std::cerr << "Error: malformed ec file compared to index; ec file: " << opt.input_directory + "/matrix.ec, index: " << opt.index << std::endl;
    exit(1);
  }
}



void PseudoQuant(KmerIndex &index, const ProgramOptions& opt, const std::string &call, const std::string &start_time) {
  // load ec matrix
  LoadEcFromFile(index, opt.input_directory + "/matrix.ec", opt);
  
  // read batch cell ids;
  std::vector<std::string> batch_ids;  
  { 
    std::unordered_set<std::string> batch_set; 
    std::ifstream in(opt.input_directory + "/matrix.cells");
    std::string line;
    while (std::getline(in,line)) {
      if (batch_set.count(line) != 0) {
        std::cerr << "Error: batch id " << line << " is repeated" << std::endl;
        exit(1);
      }
      batch_set.insert(line);
      batch_ids.push_back(line);
    }
    batch_set.clear();
  }


  // load count data
  // TODO, do this in a streaming fashion without loading all the data into memory._pos

  int cellcount = batch_ids.size();
  int numec = index.ecmap.size();

  std::vector<std::vector<int>> bigCount; // bigCount[cellid][eccid] = count
  bigCount.resize(cellcount,{});
  for (int i = 0; i < cellcount; i++) {
    auto &v = bigCount[i];
    v.resize(numec,0);
  }


  {
    std::ifstream in(opt.input_directory + "/matrix.tsv");
    std::string line;
    int i,j,c;
    while (std::getline(in, line)) {
      std::stringstream ss(line);
      ss >> j >> i >> c;
      if (i < 0 || i >= cellcount || j < 0 || j >= numec || c < 0 || bigCount[i][j] > 0) {
        std::cerr << "Error: count file " << opt.input_directory + "/matrix.tsv refers to cell/ec that doesn't exist or is listed twice, " 
        << "cellid: " << i << ", ec: " << j << ", count: " << c << std::endl;
        exit(1);
      }
      bigCount[i][j] = c;
    }
  }


  // precompute fragment length distribution
  std::vector<double> fl_means;
  std::vector<int> fld;
  std::vector<int> preBias(4096,1);
  { 
    auto mean_fl = opt.fld;
    auto sd_fl = opt.sd;
    MinCollector collection(index, opt);
    collection.init_mean_fl_trunc( mean_fl, sd_fl );
    fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);
    fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);
  }


  
  // for each cell, we can make this parallel ... later
  for (int i = 0; i < cellcount; i++) {
    MinCollector collection(index, opt);
    collection.counts = bigCount[i]; // copy counts

    EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
    em.run(10000, 50, true, false); // no bias

    H5Writer writer;

        
    // setting num_processed to 0 because pseudoquant
    writer.init(opt.output + "/abundance_" + batch_ids[i] + ".h5", 0, 0, fld, preBias, em.post_bias_, 6,
        index.INDEX_VERSION, call, start_time);
    writer.write_main(em, index.target_names_, index.target_lens_);
    plaintext_writer(opt.output + "/abundance_" + batch_ids[i] + ".tsv", em.target_names_,
        em.alpha_, em.eff_lens_, index.target_lens_);
  }

  

}

