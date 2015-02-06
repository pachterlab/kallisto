#ifndef KALLISTO_INSPECTINDEX_H
#define KALLISTO_INSPECTINDEX_H

#include <iostream>

#include "KmerIndex.h"

using namespace std;

void InspectIndex(const KmerIndex& index) {
  cout << "#[inspect] Index version number = " << index.INDEX_VERSION << endl;
  cout << "#[inspect] k = " << index.k << endl;;
  cout << "#[inspect] number of transcripts = " << index.num_trans << endl;

  cout << "#[inspect] number of equivalence classes = " << index.ecmap.size() << endl;

  // bool verify_ecmap;

  /*
  	cout << endl << endl << "#ecid\tnumtrans\ttranslist\n";
  for (auto &ecv : index.ecmap) {
  	cout << ecv.first << "\t" << ecv.second.size() << "\t";
  	int j = 0;
  	for (auto &tid : ecv.second) {
  		if (j>0) {
  			cout << ",";
  		}
  		cout << tid;
  		j++;
  	}
  	cout << "\n";
  }
  */

  cout << "#[inspect] Number of k-mers in index = " << index.kmap.size() << endl;
  int num_unique = 0, num_common = 0;
  int limit = 30;

  for (auto& kv : index.kmap) {
    auto search = index.ecmap.find(kv.second);
    if (search == index.ecmap.end()) {
      cerr << "Kmer " << kv.first.toString() << " mapped to ec " << kv.second << ", which is not in the index" << endl;
      exit(1);
    } else {
      if (search->second.size() == 1) {
        num_unique++;
      } else if (search->second.size() >= limit) {
        num_common++;
      }
    }
  }
  cout << "#[inspect] Found " << num_unique << " k-mers that map to a unique transcript" << endl;
  cout << "#[inspect] Found " << num_common << " k-mers that map to more than " << limit << " transcripts" << endl;


}

#endif // KALLISTO_INSPECTINDEX_H
