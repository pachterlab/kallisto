#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <stdint.h>
#include <regex>

#include "KmerIndex.h"
#include "common.h"

bool MergeBatchDirectories(const ProgramOptions &opt, int& num_targets, int64_t& num_processed, 
                           int64_t& num_pseudoaligned, int64_t&  num_unique, int& index_version) {
  int ndirs = opt.files.size();
  std::vector<int> numCells;
  std::unordered_set<std::string> cellSet;
  std::vector<std::string> cells;

  // stats to collect
  num_targets = -1;
  num_processed = 0;
  num_pseudoaligned = 0;
  num_unique = 0;
  index_version = -1;

  // get stat info from run_info.json
  for (const auto &fn : opt.files) {
    std::ifstream jfile(fn + "/run_info.json");
    std::ostringstream ss;
    ss << jfile.rdbuf();
    std::string s = ss.str();
    jfile.close();

    auto match_int = [&](std::string name, int64_t& val ) -> bool {
      std::string re_s = "\"" + name + "\": (\\d+),";
      std::regex re(re_s);
      std::smatch m;
      if (std::regex_search(s, m, re)) {
        if (m.empty() || m.size() <= 1) {
          return false;
        }
        val = std::stoi(m[1]);
        return true;
      } else {
        return false;
      }
    };

    int64_t tmp = -1;
    if (match_int("n_processed", tmp)) {
      num_processed += tmp;
    }
    if (match_int("n_pseudoaligned", tmp)) {
      num_pseudoaligned += tmp;
    }
    if (match_int("n_unique", tmp)) {
      num_unique += tmp;
    }
    if (index_version == -1 && match_int("index_version",tmp)) {
      index_version = (int) tmp;
    }
    if (num_targets == -1 && match_int("n_targets", tmp)) {
      num_targets = (int) tmp;
    }
          
  }  

  // merge cell lists
  for (const auto &fn : opt.files) {
    std::ifstream cfile((fn + "/matrix.cells"));
    
    std::string line;
    std::string id,f1,f2;
    int num = 0;
    while (std::getline(cfile,line)) {
      if (line.size() == 0) {
        continue;
      }
      std::stringstream ss(line);
      ss >> id;

      if (cellSet.find(id) == cellSet.end()) {
        cells.push_back(id);
        cellSet.insert(id);
        num++;
      } else {
        std::cerr << "Error: cell id \"" << id << "\" repeated in file " << fn << "/matrix.cells" << std::endl;
        return false;
      }
    }
    numCells.push_back(num);
  }

  std::ofstream cfile((opt.output + "/matrix.cells"));
  for (const auto &id : cells) {
    cfile << id << "\n";
  }
  cfile.close();

  // merge ec lists
  KmerIndex index(opt);
  std::vector<std::vector<int>> ecTrans;
  std::vector<int> numEcs;

  for (int i = 0; i < ndirs; i++) {
    std::vector<int> ctrans;
    std::ifstream ecFile(opt.files[i] + "/matrix.ec");

    std::string line,t;    
    std::vector<int> c;
    while (std::getline(ecFile,line)) {   
      c.clear();
      int ec = -1;
      if (line.size() == 0) {
        continue;
      }
      std::stringstream ss(line);
      ss >> ec;
      while (std::getline(ss, t, ',')) {
        c.push_back(std::stoi(t));
      }

      assert(ctrans.size() == ec);

      int index_ec = -1;
      auto search = index.ecmapinv.find(c);
      if (search != index.ecmapinv.end()) {
        index_ec = search->second;
      } else {        
         index_ec = index.ecmap.size();
         index.ecmap.push_back(c);
         index.ecmapinv.insert({c,index_ec});         
      }
      ctrans.push_back(index_ec);            
    }
    numEcs.push_back(ctrans.size());
    ecTrans.push_back(std::move(ctrans));
  }

  //debug
  /*
  for (auto &ct : ecTrans) {
    std::cout << "-- " << std::endl;
    for (int i = 0; i < ct.size(); i++) {
      std::cout << i << "\t" << ct[i] << std::endl;
    }
  }
  */

  std::ofstream ecFile(opt.output + "/matrix.ec");
  for (int j = 0; j < index.ecmap.size(); j++) {
    ecFile << j << "\t";
    // output the rest of the class
    const auto &v = index.ecmap[j];
    bool first = true;
    for (auto x : v) {
      if (!first) {
        ecFile << ",";
      } else {
        first = false;
      }
      ecFile << x;
    }
    ecFile << "\n";
  }
  ecFile.close();

  // get stats
  std::vector<int64_t> numNonzero;
  for (int i = 0; i < ndirs; i++) {
    std::ifstream mFile(opt.files[i] + "/matrix.tcc.mtx");
    std::string line;
    std::getline(mFile, line);
    assert(!line.empty() && line[0] =='%');
    int ncell,nec;
    int64_t nz;
    mFile >> ncell >> nec >> nz;
    mFile.close();

    if (ncell != numCells[i]) {
      std::cerr << "Error: number of cells in mtx and cell file do not match in directory  " << opt.files[i] << std::endl;
      return false; 
    }

    if (nec != numEcs[i]) {
      std::cerr << "Error: number of equivalence classes in mtx and cell file do not match in directory  " << opt.files[i] << std::endl;
      return false; 
    }
    numNonzero.push_back(nz);
  }


  // merge into single file
  std::ofstream omFile(opt.output + "/matrix.tcc.mtx");
  omFile << "%%MatrixMarket matrix coordinate real general\n";
  int sNumCells = std::accumulate(numCells.begin(), numCells.end(), 0);
  int sNumEcs = index.ecmap.size();
  int64_t sNz = std::accumulate(numNonzero.begin(), numNonzero.end(), 0);
  omFile << sNumCells << "\t" << sNumEcs << "\t" << sNz << "\n";
  int64_t nz = 0;
  int cellOffset = 0;
  for (int i = 0; i < ndirs; i++) {
    // merge all files
    std::ifstream mFile(opt.files[i] + "/matrix.tcc.mtx");
    std::string line;
    std::getline(mFile, line); // mtx format
    std::getline(mFile, line); // specification
    if (i > 0) {
      cellOffset += numCells[i-1];
    }

    int local_cellid, local_ec, count;
    while (mFile >> local_cellid >> local_ec >> count) {
      local_ec -= 1;
      int ec = ecTrans[i][local_ec];
      int cellid = cellOffset + local_cellid;
      nz++;
      omFile << cellid << "\t" << (ec+1) << "\t" << count << "\n";
    }
  }

  if (nz != sNz) {
    std::cerr << "Warning: number of matrix entries inconsistent with header" << std::endl;
  }
  omFile.close();

  

  return true;
}



