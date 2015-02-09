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
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>

//#include <map>


#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"

#include "KmerHashTable.h"

#include "hash.hpp"

using EcMap = std::unordered_map<int, std::vector<int>>;

struct SortedVectorHasher {
  size_t operator()(const std::vector<int>& v) const {
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

struct KmerIndex {
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0), skip(opt.skip) {
    //LoadTranscripts(opt.transfasta);
  }

  KmerIndex() : num_trans(0) {
      // minimal constructor for loading python
  }

  ~KmerIndex() {}


  // use:  match(s,l,v)
  // pre:  v is initialized
  // post: v contains all equiv classes for the k-mers in s
  void match(const char *s, int l, std::vector<int>& v) const {
    KmerIterator kit(s), kit_end;
    for (int i = 0; kit != kit_end; ++kit,++i) {
      if (i==skip) {
        i=0;
      }
      if (i==0) {
        Kmer rep = kit->first.rep();
        auto search = kmap.find(rep);
        if (search != kmap.end()) {
          // if k-mer found
          v.push_back(search->second); // add equivalence class
        }
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
      auto& u = search->second;
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
      target_names_.push_back(seq->name.s);
      trans_lens_.push_back(seq->seq.l);
      
      // if it is long enough
      if (seq->seq.l >= k) {
        KmerIterator kit(seq->seq.s), kit_end;
        // for each k-mer add to map
        for(; kit != kit_end; ++kit) {
          Kmer rep = kit->first.rep();
          kmcount[rep]++;
          auto search = all_kmap.find(rep);
          if (search == all_kmap.end()) {
            // new k-mer
            all_kmap.insert({rep, {transid}});
          } else {
            // seen before
            std::vector<int>& v = search->second;
            if (*v.rbegin() < transid) {
              // but new transcript
              v.push_back(transid);
            }
          }
        }
      }

      transid++;
      if (transid % 1000 == 0) {
        std::cerr << " " << transid << " size of k-mer map " << all_kmap.size() << std::endl;
      }

    }

    num_trans = transid;
    std::cerr << "Found " << num_trans << " transcripts"
              << std::endl
              << "Size of k-mer map " << all_kmap.size() << std::endl;


    // for each transcript
    for (int i = 0; i < num_trans; i++ ) {
      // create its own eqs
      std::vector<int> single(1,i);
      ecmap.insert({i,single});
      ecmapinv.insert({single,i});
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

  void write(const std::string& index_out, bool writeKmerTable = true) {
    std::ofstream out;
    out.open(index_out, std::ios::out | std::ios::binary);

    if (!out.is_open()) {
      // TODO: better handling
      std::cerr << "Error: index output file could not be opened!";
      exit(1);
    }

    // 1. write index
    out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

    // 2. write k
    out.write((char *)&k, sizeof(k));

    // 3. write number of transcripts
    out.write((char *)&num_trans, sizeof(num_trans));

    // 4. write out transcript lengths
    for (int tlen : trans_lens_) {
      out.write((char *)&tlen, sizeof(tlen));
    }

    size_t kmap_size = kmap.size();

    if (writeKmerTable) {
      // 5. write number of k-mers in map
      out.write((char *)&kmap_size, sizeof(kmap_size));

      // 6. write kmer->ec values
      for (auto& kv : kmap) {
        out.write((char *)&kv.first, sizeof(kv.first));
        out.write((char *)&kv.second, sizeof(kv.second));
      }
    } else {
      // 5. write fake k-mer size
      kmap_size = 0;
      out.write((char *)&kmap_size, sizeof(kmap_size));

      // 6. write none of the kmer->ec values
    }
    // 7. write number of equivalence classes
    size_t tmp_size;
    tmp_size = ecmap.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 8. write out each equiv class
    for (auto& kv : ecmap) {
      out.write((char *)&kv.first, sizeof(kv.first));

      // 8.1 write out the size of equiv class
      tmp_size = kv.second.size();
      out.write((char *)&tmp_size, sizeof(tmp_size));
      // 8.2 write each member
      for (auto& val: kv.second) {
        out.write((char *)&val, sizeof(val));
      }
    }

    // 9. Write out target ids
    // XXX: num_trans should equal to target_names_.size(), so don't need
    // to write out again.
    assert(num_trans == target_names_.size());
    for (auto& tid : target_names_) {
      // 9.1 write out how many bytes
      // XXX: Note: this doesn't actually encore the max targ id size.
      // might cause problems in the future
      tmp_size = tid.size();
      out.write((char *)&tmp_size, sizeof(tmp_size));

      // 9.2 write out the actual string
      out.write(tid.c_str(), tid.size());
    }

    out.flush();
    out.close();
  }

  // note opt is not const
  void load(ProgramOptions& opt, bool loadKmerTable = true) {

    std::string& index_in = opt.index;
    std::ifstream in;

    in.open(index_in, std::ios::in | std::ios::binary);

    if (!in.is_open()) {
      // TODO: better handling
      std::cerr << "Error: index input file could not be opened!";
      exit(1);
    }

    // 1. read version
    size_t header_version = 0;
    in.read((char *)&header_version, sizeof(header_version));

    if (header_version != INDEX_VERSION) {
      std::cerr << "Error: Incompatiple indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
                << "Rerun with index to regenerate!";
      exit(1);
    }

    // 2. read k
    in.read((char *)&k, sizeof(k));
    if (Kmer::k == 0) {
      //std::cerr << "[index] no k has been set, setting k = " << k << std::endl;
      Kmer::set_k(k);
      opt.k = k;
    } else if (Kmer::k == k) {
      //std::cerr << "[index] Kmer::k has been set and matches" << k << std::endl;
      opt.k = k;
    } else {
      std::cerr << "Error: Kmer::k was already set to = " << Kmer::k << std::endl
                << "       conflicts with value of k  = " << k << std::endl;
      exit(1);
    }

    // 3. read number of transcripts
    in.read((char *)&num_trans, sizeof(num_trans));

    // 4. read number of transcripts
    trans_lens_.clear();
    trans_lens_.reserve(num_trans);

    for (int i = 0; i < num_trans; i++) {
      int tlen;
      in.read((char *)&tlen, sizeof(tlen));
      trans_lens_.push_back(tlen);
    }

    // 5. read number of k-mers
    size_t kmap_size;
    in.read((char *)&kmap_size, sizeof(kmap_size));

    std::cerr << "[index] k: " << k << std::endl;
    std::cerr << "[index] num_trans read: " << num_trans << std::endl;
    std::cerr << "[index] kmap size: " << kmap_size << std::endl;

    kmap.clear();
    if (loadKmerTable) {
      kmap.reserve(kmap_size);
    }

    // 6. read kmer->ec values
    Kmer tmp_kmer;
    int tmp_val;
    for (size_t i = 0; i < kmap_size; ++i) {
      in.read((char *)&tmp_kmer, sizeof(tmp_kmer));
      in.read((char *)&tmp_val, sizeof(tmp_val));

      if (loadKmerTable) {
        kmap.insert({tmp_kmer, tmp_val});
      }
    }

    // 7. read number of equivalence classes
    size_t ecmap_size;
    in.read((char *)&ecmap_size, sizeof(ecmap_size));

    std::cerr << "[index] ecmap size: " << ecmap_size << std::endl;

    int tmp_id;
    size_t vec_size;
    // 8. read each equiv class
    for (size_t i = 0; i < ecmap_size; ++i) {
      in.read((char *)&tmp_id, sizeof(tmp_id));

      // 8.1 read size of equiv class
      in.read((char *)&vec_size, sizeof(vec_size));

      // 8.2 read each member
      std::vector<int> tmp_vec;
      tmp_vec.reserve(vec_size);
      for (size_t j = 0; j < vec_size; ++j ) {
        in.read((char *)&tmp_val, sizeof(tmp_val));
        tmp_vec.push_back(tmp_val);
      }
      ecmap.insert({tmp_id, tmp_vec});
      ecmapinv.insert({tmp_vec, tmp_id});
    }

    // 9. read in target ids
    target_names_.clear();
    target_names_.reserve(num_trans);

    size_t tmp_size;
    char buffer[1024]; // if your target_name is longer than this, screw you.
    for (auto i = 0; i < num_trans; ++i) {
      // 9.1 read in the size
      in.read((char *)&tmp_size, sizeof(tmp_size));

      // 9.2 read in the character string
      in.read(buffer, tmp_size);

      std::string tmp_targ_id( buffer );
      target_names_.push_back(std::string( buffer ));

      // clear the buffer for next string
      memset(buffer,0,strlen(buffer));
    }

    in.close();
  }

  int k; // k-mer size used
  int num_trans; // number of transcripts
  int skip;
  //std::unordered_map<Kmer, int, KmerHash> kmap;
  KmerHashTable<int, KmerHash> kmap;

  EcMap ecmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  const size_t INDEX_VERSION = 5; // increase this every time you change the fileformat

  std::vector<int> trans_lens_;

  std::vector<std::string> target_names_;
};



#endif // KALLISTO_KMERINDEX_H
