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
  int k = index.k;
  cout << "#[inspect] Index version number = " << index.INDEX_VERSION << endl;
  cout << "#[inspect] k = " << index.k << endl;;
  cout << "#[inspect] number of targets = " << index.num_trans << endl;

  cout << "#[inspect] number of equivalence classes = " << index.ecmap.size() << endl;


  if (index.ecmap.size() != index.ecmapinv.size()) {
    cout << "Error: sizes do not match. ecmap.size = " << index.ecmap.size()
         << ", ecmapinv.size = " << index.ecmapinv.size() << endl;
    exit(1);
  }

  if (index.dbGraph.ecs.size() != index.dbGraph.contigs.size()) {
    cout << "Error: sizes do not match. ecs.size = " << index.dbGraph.ecs.size()
         << ", contigs.size = " << index.dbGraph.contigs.size() << endl;
    exit(1);
  }

  cout << "#[inspect] number of contigs = " << index.dbGraph.contigs.size() << endl;
  

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
        cout << "Error: ec = " << ec  << " has invalid target id " << v[i] << endl;
        exit(1);
      }

      if (i > 0 && v[i] == v[i-1]) {
        cout << "Error: ec = " << ec  << " has repeated target id " << v[i] << endl;
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

  for (auto& kv : index.kmap) {
    int id = kv.second.contig;
    int pos = kv.second.getPos();
    int fw = kv.second.isFw();

    if (id < 0 || id >= index.dbGraph.contigs.size()) {
      cerr << "Kmer " << kv.first.toString() << " mapped to contig " << id << ", which is not in the de Bruijn Graph" << endl;
      exit(1);
    } else {
      ++kmhisto[index.ecmap[index.dbGraph.ecs[id]].size()];
    }

    const Contig& c = index.dbGraph.contigs[id];
    const char* s = c.seq.c_str();
    Kmer x = Kmer(s+pos);
    Kmer xr = x.rep();

    bool bad = (fw != (x==xr)) || (xr != kv.first);
    if (bad) {
      cerr << "Kmer " << kv.first.toString() << " mapped to contig " << id << ", pos = " << pos << ", on " << (fw ? "forward" : "reverse") << " strand" << endl;
      cerr << "seq = " << c.seq << endl;
      cerr << "x  = " << x.toString() << endl;
      cerr << "xr = " << xr.toString() << endl;
      exit(1);
    }
  }

  for (int i = 0; i < index.dbGraph.contigs.size(); i++) {
    const Contig& c = index.dbGraph.contigs[i];

    if (c.seq.size() != c.length + k-1) {
      cerr << "Length and string dont match " << endl << "seq = " << c.seq << " (length = " << c.seq.size() << "), c.length = " << c.length << endl;
      exit(1);
    }


    const char *s = c.seq.c_str();
    KmerIterator kit(s), kit_end;
    for (; kit != kit_end; ++kit) {
      Kmer x = kit->first;
      Kmer xr = x.rep();
      auto search = index.kmap.find(xr);
      if (search == index.kmap.end()) {
        cerr << "could not find kmer " << x.toString() << " in map " << endl << "seq = " << c.seq << ", pos = " << kit->second << endl;
        exit(1);
      }

      KmerEntry val = search->second;
      if (val.contig != i || val.getPos() != kit->second || val.isFw() != (x==xr)) {
        cerr << "mismatch " << x.toString() << " in map " << endl << "id = " << i << ", seq = " << c.seq << ", pos = " << kit->second << endl;
        cerr << "val = " << val.contig << ", pos = (" << val.getPos() << ", " << (val.isFw() ? "forward" :  "reverse") << ")" << endl;
        exit(1);
      }
    }
    
  }

  printHisto(echisto, "#EC.size\tNum.targets");
  cout << endl << endl;

  printHisto(kmhisto, "#EC.size\tNum.kmers");



}

#endif // KALLISTO_INSPECTINDEX_H
