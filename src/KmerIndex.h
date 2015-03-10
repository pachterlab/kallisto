#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include <seqan/sequence.h>
#include <seqan/index.h>


//#include <map>


#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"

#include "KmerHashTable.h"

#include "hash.hpp"

typedef seqan::Index<seqan::StringSet<seqan::CharString>, seqan::IndexSa<>> TIndex;
namespace SEQAN_NAMESPACE_MAIN
{
  
template<>
struct SAValue<TIndex>
{
  typedef Pair<unsigned, unsigned, Pack> Type;
};

}

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

  ~KmerIndex() {}


  // use:  match(s,l,v)
  // pre:  v is initialized
  // post: v contains all equiv classes for the k-mers in s
  void match(const char *s, int l, std::vector<std::pair<int, int>>& v) const {
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
          v.push_back({search->second, kit->second}); // add equivalence class, and position
        }
      }
    }
  }

  int mapPair(const char* s1, int l1, const char* s2, int l2, int ec) const;

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

  void BuildTranscripts(const ProgramOptions& opt);

  void write(const std::string& index_out, bool writeKmerTable = true);
  void load(ProgramOptions& opt, bool loadKmerTable = true);
  bool loadSuffixArray(const ProgramOptions& opt);
  void clearSuffixArray();
  // note opt is not const

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

  // Suffix array indices

  typedef seqan::Finder<TIndex> TFinder;

  TIndex index;
  seqan::StringSet<seqan::CharString> seqs;
  
};





#endif // KALLISTO_KMERINDEX_H
