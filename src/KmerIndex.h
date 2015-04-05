#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <stdint.h>

#include <seqan/sequence.h>
#include <seqan/index.h>


//#include <map>


#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"

#include "KmerHashTable.h"

#include "hash.hpp"

typedef seqan::Index<seqan::StringSet<seqan::CharString>, seqan::IndexSa<>> TIndex;
typedef seqan::Finder<TIndex> TFinder;

namespace seqan {

template<>
//struct SAValue<TIndex>
struct SAValue<seqan::StringSet<seqan::CharString>> {
  typedef Pair<uint32_t, uint32_t> Type;
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

struct KmerEntry {
  int32_t id;    // id of equivalence class
  int16_t fdist; // 0-based forward distance to EC-junction
  int16_t bdist; // 0-based backward distance to EC-junction

  KmerEntry() : id(-1), fdist(-1), bdist(-1) {}
  KmerEntry(int _id) : id(_id), fdist(-1), bdist(-1) {}
  KmerEntry(int _id, int _fdist, int _bdist) : id(_id), fdist(_fdist), bdist(_bdist) {}
};



struct KmerIndex {
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0), skip(opt.skip) {
    //LoadTranscripts(opt.transfasta);
  }

  ~KmerIndex() {}

  void match(const char *s, int l, std::vector<std::pair<int, int>>& v) const;
  int mapPair(const char *s1, int l1, const char *s2, int l2, int ec, TFinder& finder) const;
  std::vector<int> intersect(int ec, const std::vector<int>& v) const;




  void BuildTranscripts(const ProgramOptions& opt);
  bool fwStep(Kmer km, Kmer& end, int ec) const;
  void write(const std::string& index_out, bool writeKmerTable = true);
  // note opt is not const
  void load(ProgramOptions& opt, bool loadKmerTable = true);


  bool loadSuffixArray(const ProgramOptions& opt);
  void clearSuffixArray();




  int k; // k-mer size used
  int num_trans; // number of transcripts
  int skip;
  //std::unordered_map<Kmer, int, KmerHash> kmap;
  KmerHashTable<KmerEntry, KmerHash> kmap;

  EcMap ecmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  const size_t INDEX_VERSION = 6; // increase this every time you change the fileformat

  std::vector<int> trans_lens_;

  std::vector<std::string> target_names_;

  // Suffix array indices
  TIndex index;
  //TFinder finder;

  seqan::StringSet<seqan::CharString> seqs;

};





#endif // KALLISTO_KMERINDEX_H
