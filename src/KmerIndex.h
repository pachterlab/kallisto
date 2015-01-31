#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H

#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
//#include <map>


#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"


#include "hash.hpp"


struct SortedVectorHasher {
        size_t operator()(const std::vector<int> &v) const {
                uint64_t r = 0;
                int i=0;
                for (auto x : v) {
                        uint64_t t;
                        MurmurHash3_x64_64(&x,sizeof(x), 0,&t);
                        t = (x>>i) | (x<<(64-i));
                        r = r ^ t;
                        i = (i+1)%64;
                }
                return r;
        }
};

struct KmerIndex
{
    KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0) {
		//LoadTranscripts(opt.transfasta);
	}

	~KmerIndex() {}


	// use:  match(s,l,v)
	// pre:  v is initialized
	// post: v contains all equiv classes for the k-mers in s
	void match(const char *s, int l, std::vector<int> & v) const {
		KmerIterator kit(s), kit_end;
		for (;kit != kit_end; ++kit) {
			Kmer rep = kit->first.rep();
			auto search = kmap.find(rep);
			if (search != kmap.end()) {
				// if k-mer founc
				v.push_back(search->second); // add equivalence class
			}
		}
	}

	// use:  res = intersect(ec,v)
	// pre:  ec is in ecmap, v is a vector of valid transcripts
	//       v is sorted in increasing order
	// post: res contains the intersection  of ecmap[ec] and v sorted increasing
	//       res is empty if ec is not in ecmap
	std::vector<int> intersect(int ec, const std::vector<int>& v) const {
		std::vector<int> res;
		auto search = ecmap.find(ec);
		if (search != ecmap.end()) {
			auto &u = search->second;
			res.reserve(v.size());

			auto a = u.begin();
			auto b = v.begin();

			while (a != u.end() && b != v.end()) {
				if (*a < *b) {
					++a;
				} else if (*b < *a) {
					++b;
				} else {
					// match
					res.push_back(*a);
					++a;
					++b;
				}
			}
		}
		return res;
	}


	void BuildTranscripts(const std::string& fasta) {
		// TODO: add code to check if binary file exists and load it directly
		// FIXME: check if FASTA file actually exists
		// If it doesn't, will just hang
		int l;
		std::cerr << "Loading fasta file " << fasta
							<< std::endl;
        std::cerr << "k: " << k << std::endl;
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
				if (transid % 1000 == 1) {
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


		std::cerr << "K-mer map has " << kmap.size() << " k-mers and " << std::endl;
		kseq_destroy(seq);
		gzclose(fp);
	}

    void write(const std::string& index_out)
    {
        std::ofstream out;
        out.open(index_out, std::ios::out | std::ios::binary);

        if (!out.is_open()) {
            // TODO: better handling
            std::cerr << "index output file could not be open!";
            exit(1);
        }

        // TODO: add version to index
        out.write((char*)&k, sizeof(k));
        out.write((char*)&num_trans, sizeof(num_trans));

        // print out the size of the map
        size_t kmap_size = kmap.size();
        /* std::cerr << "orig kmap size " << kmap_size << '\t' << kmap.size() << std::endl; */
        out.write((char*)&kmap_size, sizeof(kmap_size));


        for (auto& kv : kmap) {
            out.write((char*)&kv.first, sizeof(kv.first));
            out.write((char*)&kv.second, sizeof(kv.second));
        }

        out.flush();

        size_t tmp_size;
        tmp_size = ecmap.size();
        out.write((char*)&tmp_size, sizeof(tmp_size));

        for (auto& kv : ecmap) {
            out.write((char*)&kv.first, sizeof(kv.first));

            tmp_size = kv.second.size();
            out.write((char*)&tmp_size, sizeof(tmp_size));

            for (auto& val: kv.second) {
                out.write((char*)&val, sizeof(val));
            }
        }


        out.flush();

        out.close();
    }

	void load(const std::string& index_in)
    {
        std::ifstream in;

        in.open(index_in, std::ios::in | std::ios::binary);

        if (!in.is_open()) {
            // TODO: better handling
            std::cerr << "index input file could not be open!";
            exit(1);
        }

        in.read((char*)&k, sizeof(k));
        in.read((char*)&num_trans, sizeof(num_trans));

        size_t kmap_size;
        in.read((char*)&kmap_size, sizeof(kmap_size));

        std::cerr << "[index] k: " << k << std::endl;
        std::cerr << "[index] num_trans read: " << num_trans << std::endl;
        std::cerr << "[index] kmap size: " << kmap_size << std::endl;

        // read in the kmap
	kmap.clear();
	kmap.reserve(kmap_size);
        Kmer tmp_kmer;
        int tmp_val;
        for (size_t i = 0; i < kmap_size; ++i)
        {
            in.read((char*)&tmp_kmer, sizeof(tmp_kmer));
            in.read((char*)&tmp_val, sizeof(tmp_val));

            // if (i < 15)
            //     std::cout << "\t\tval:\t" << tmp_kmer << "\t" << tmp_val << std::endl;
            kmap.insert({tmp_kmer, tmp_val});
        }

        // read in the ecmap
        size_t ecmap_size;
        in.read((char*)&ecmap_size, sizeof(ecmap_size));

        std::cerr << "[index] ecmap size: " << ecmap_size << std::endl;

        int tmp_id;
        size_t vec_size;

        for (size_t i = 0; i < ecmap_size; ++i) {
            in.read((char*)&tmp_id, sizeof(tmp_id));
            in.read((char*)&vec_size, sizeof(vec_size));

            std::vector<int> tmp_vec(vec_size);
            for (size_t j = 0; j < vec_size; ++j )
            {
                in.read((char*)&tmp_val, sizeof(tmp_val));
                tmp_vec.push_back(tmp_val);
            }
            ecmap.insert({tmp_id, tmp_vec});
            ecmapinv.insert({tmp_vec, tmp_id});
        }

        in.close();
    }

	int k; // k-mer size used
	int num_trans; // number of transcripts
	std::unordered_map<Kmer, int, KmerHash> kmap;
	std::unordered_map<int, std::vector<int>> ecmap;
	std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;

	// TODO: include lengths of transctipts
    std::vector<unsigned int> trans_lens_;
};

#endif // KALLISTO_KMERINDEX_H
