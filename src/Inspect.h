#ifndef KALLISTO_INSPECTINDEX_H
#define KALLISTO_INSPECTINDEX_H

#include <iostream>

#include "KmerIndex.h"

using namespace std;

void printVector(const vector<int>& v) {
  cout << "[";
  int i = 0;
  for (auto x : v) {
    if (i>0) {
      cout << ", ";
    }
    cout << x;
    i++;
  }
  cout << "]";
}

void printHisto(const unordered_map<int,int>& m, const string& header) {
  cout << header << "\n";
  int mn = std::numeric_limits<int>::max();
  int mx = 0;

  for (auto kv : m) {
    mn = min(mn,kv.first);
    mx = max(mx,kv.first);
  }

  for (int i = mn; i <= mx; i++) {
    auto search = m.find(i);
    if (search == m.end()) {
      cout << i << "\t0\n";
    } else {
      cout << i << "\t" << search->second << "\n";
    }
  }
}

void InspectIndex(const KmerIndex& index) {
  cout << "#[inspect] Index version number = " << index.INDEX_VERSION << endl;
  cout << "#[inspect] k = " << index.k << endl;;
  cout << "#[inspect] number of transcripts = " << index.num_trans << endl;

  cout << "#[inspect] number of equivalence classes = " << index.ecmap.size() << endl;

  bool verify_ecmap;

  if (index.ecmap.size() != index.ecmapinv.size()) {
    cout << "Error: sizes do not match. ecmap.size = " << index.ecmap.size()
         << ", ecmapinv.size = " << index.ecmapinv.size() << endl;
    exit(1);
  }

  unordered_map<int,int> echisto;

  //for (auto& ecv : index.ecmap) {
  for (int ec = 0; ec < index.ecmap.size(); ec++) {
    const vector<int>& v = index.ecmap[ec];
    ++echisto[v.size()];

    if (v.empty()) {
      cout << "Error: ec = " << ec  << " is empty!" << endl;
      exit(1);
    }
    for (int i = 0; i < v.size(); i++) {
      if (v[i] < 0 || v[i] >= index.num_trans) {
        cout << "Error: ec = " << ec  << " has invalid transcript id " << v[i] << endl;
        exit(1);
      }

      if (i > 0 && v[i] == v[i-1]) {
        cout << "Error: ec = " << ec  << " has repeated transcript id " << v[i] << endl;
        exit(1);
      }

      if (i > 0 && v[i] < v[i-1]) {
        cout << "Error: ec = " << ec  << " is not sorted!" << endl;
        exit(1);
      }
    }

    auto search = index.ecmapinv.find(v);
    if (search == index.ecmapinv.end()) {
      cout << "Error: could not find inverse for " << ec << endl;
      exit(1);
    } else {
      if (search->second != ec) {
        cout << "Error: inverse incorrect for ecmap -> ecmapinv,  ecv.first = "
             << ec <<  ", ecmapinv[ecv.second] = " << search->second << endl;
        exit(1);
      }
    }
  }

  for (auto& eiv : index.ecmapinv) {
    //auto search = index.ecmap.find(eiv.second);
    //if (search == index.ecmap.end()) {
    if (eiv.second < 0 || eiv.second >= index.ecmap.size()) {
      cout << "Error: could not find inverse for ";
      printVector(eiv.first);
      cout << ", ecid = " << eiv.second << endl;
      exit(1);
    } else {
      auto &v = index.ecmap[eiv.second];
      if (v != eiv.first) {
        cout << "Error: inverse incorrect for ecmapinv -> ecmap,  eiv.first = ";
        printVector(eiv.first);
        cout <<  ", ecmap[eiv.second] = ";
        printVector(v);
        cout << endl;
        exit(1);
      }
    }
  }

  cout << "#[inspect] Number of k-mers in index = " << index.kmap.size() << endl;
  unordered_map<int,int> kmhisto;
  unordered_map<int,int> fjumphisto;
  unordered_map<int,int> bjumphisto;

  for (auto& kv : index.kmap) {
    ++fjumphisto[kv.second.fdist];
    ++bjumphisto[kv.second.bdist];
    //auto search = index.ecmap.find(kv.second.id);
//    if (search == index.ecmap.end()) {
    if (kv.second.id < 0 || kv.second.id >= index.ecmap.size()) {
      cerr << "Kmer " << kv.first.toString() << " mapped to ec " << kv.second.id << ", which is not in the index" << endl;
      exit(1);
    } else {
      ++kmhisto[index.ecmap[kv.second.id].size()];
    }
  }


  printHisto(echisto, "#EC.size\tNum.transcripts");
  cout << endl << endl;

  printHisto(kmhisto, "#EC.size\tNum.kmers");

  cout << endl << endl;
  printHisto (fjumphisto, "#Jump.fw\tNum.kmers");

  cout << endl << endl;
  printHisto (bjumphisto, "#Jump.bw\tNum.kmers");

}

#endif // KALLISTO_INSPECTINDEX_H
