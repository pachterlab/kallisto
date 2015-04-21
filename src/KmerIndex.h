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

struct TRInfo {
  int trid;
  int start;
  int stop; //exclusive [start,stop)
};

using EcMap = std::vector<std::vector<int>>; //std::unordered_map<int, std::vector<int>>;

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
  int32_t contig; // id of contig
  uint32_t _pos; // 0-based forward distance to EC-junction

  KmerEntry() : contig(-1), _pos(0xFFFFFFF) {}
  KmerEntry(int id, int pos, bool isFw) : contig(id) {
    setPos(pos);
    setDir(isFw);
  }

  inline int getPos() {return (_pos & 0x0FFFFFFF);}
  inline int isFw()   {return (_pos & 0xF0000000) == 0; }
  inline void setPos(int p) {_pos = (_pos & 0xF0000000) | (p & 0x0FFFFFFF);}
  inline void setDir(bool _isFw) {_pos = (_pos & 0x0FFFFFFF) | ((_isFw) ? 0 : 0xF0000000);}
};

struct Contig {
  int id; // internal id
  int length; // number of k-mers
  std::string seq; // sequence

  
  int getDist(KmerEntry x, bool fw) const {
    if (x.isFw() == fw ) {
      return (length -1 - x.getPos());
    } else {
      return x.getPos();
    }
  }
};

struct DBGraph {
  std::vector<int> ecs; // contig id -> ec-id
  std::vector<Contig> contigs; // contig id -> contig
//  std::vector<pair<int, bool>> edges; // contig id -> edges
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
  void BuildDeBruijnGraph(const ProgramOptions& opt);
  void BuildEquivalenceClasses(const ProgramOptions& opt);
  void FixSplitContigs(const ProgramOptions& opt, const std::vector<std::vector<TRInfo>>& trinfos);
  bool fwStep(Kmer km, Kmer& end) const;
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
  DBGraph dbGraph;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  const size_t INDEX_VERSION = 7; // increase this every time you change the fileformat

  std::vector<int> trans_lens_;

  std::vector<std::string> target_names_;

  // Suffix array indices
  TIndex index;
  //TFinder finder;

  seqan::StringSet<seqan::CharString> seqs;

};








#endif // KALLISTO_KMERINDEX_H
