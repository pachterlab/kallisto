#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H

#include <vector>
#include <unordered_map>
#include <map>
#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"


struct KmerIndex
{
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0) {
		//LoadTranscripts(opt.transfasta);
	}
	~KmerIndex() {}

	void match(const char *s, int l, std::vector<int> & v) const {}

	void BuildTranscripts(const std::string fasta) {
		// TODO: add code to check if binary file exists and load it directly
		int l;
		std::cerr << "Loading fasta file " << fasta
							<< std::endl;
		gzFile fp = gzopen(fasta.c_str(),"r");
		kseq_t *seq = kseq_init(fp);
		int transid = 0;
		std::unordered_map<Kmer, int, KmerHash> kmcount; // temporary
		// maps kmers to set of transcript ids that contain them
		std::unordered_map<Kmer, std::vector<int>, KmerHash> all_kmap;
		// for each transcript in fasta file
		while ((l = kseq_read(seq)) > 0) {
			bool added = false;
			// if it is long enough
			if (seq->seq.l >= k) {
				KmerIterator kit(seq->seq.s), kit_end;
				// for each k-mer add to map
				for(;kit != kit_end; ++kit) {
					Kmer rep = kit->first.rep();
					kmcount[rep]++;
					all_kmap[rep].push_back(transid); // creates an entry if not found
					added = true;
				}
			}
			if (added) {
				transid++;
				if (transid % 100 == 1) {
					std::cerr << " " << transid << " size of k-mer map " << all_kmap.size() << std::endl;
				}
			}
		}

		num_trans = transid;
		std::cerr << "Found " << num_trans << " transcripts"
							<< std::endl
							<< "Size of k-mer map " << all_kmap.size() << std::endl;
			

		// for each transcript
		for (int i = 0; i < num_trans; i++ ) {
			// create its own eqs
			ecmap.insert({i,{i}});
			ecmapinv.insert({{i},i});
		}

		
		int eqs_id = num_trans;


		for (auto& kv : all_kmap) {
			auto search = ecmapinv.find(kv.second);
			// if we have seen this equivalence class
			if (search != ecmapinv.end()) {
				// update kmap
				kmap.insert({kv.first, search->second});
			} else {
				// else create a new equivalence class and update kmap
				ecmapinv.insert({kv.second,eqs_id});
				ecmap.insert({eqs_id, kv.second});
				kmap.insert({kv.first, eqs_id});
				eqs_id++;
			}
		}

		std::cerr << "Created " << ecmap.size() << " equivalence classes from " << num_trans << " transcripts" << std::endl;

		/* std::cout << "EqId\tTransIdList\n"; */
		/* for (auto &ekv : ecmap) { */
		/* 	std::cout << ekv.first; */
		/* 	for (auto el : ekv.second) { */
		/* 		std::cout << "\t" << el; */
		/* 	} */
		/* 	std::cout << "\n"; */
		/* } */
		/* std::cout.flush(); */


		int uniqueKmers;
		for (auto &kmv : kmap) {
			if (kmv.second < num_trans) {
				uniqueKmers++;
			}
		}
		std::cerr << "K-mer map has " << kmap.size() << " k-mers and " << uniqueKmers << " unique k-mers" << std::endl;
		std::map<int, int> histo;
		for (auto &kv : kmcount) {
			histo[kv.second]++;
		}
		int max = histo.rend()->first;
		for (int i = 1; i < max; i++) {
			std::cout << i << "\t" << histo[i] << "\n";
		}
		std::cout.flush();
			
		
	}
	
	int k; // k-mer size used
	int num_trans; // number of transcripts
	std::unordered_map<Kmer, int, KmerHash> kmap;
	std::map<int, std::vector<int>> ecmap;
	std::map<std::vector<int>, int> ecmapinv; 
	
};

#endif // KALLISTO_KMERINDEX_H
