#ifndef GENE_MODEL_HPP
#define GENE_MODEL_HPP

#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <unordered_map>
#include "common.h"
#include "KmerIndex.h"



struct ExonModel {
  int chr;
  int start,stop;
  bool strand;
};

struct TranscriptModel {
  int id;
  int chr;
  int start,stop;
  int length;
  std::string name; // like ENST  
  std::vector<ExonModel> exons;
  bool strand;
  int gene_id;
};

struct TranscriptAlignment {
  int chr; // chr number i based on bam header
  int chrpos; // 0-based position of left-most bp
  bool strand;
  std::vector<uint32_t> cigar; // the cigar ops that determine the alignment

  TranscriptAlignment() : chr(-1), chrpos(-1), strand(true) {}

  inline bool operator==(const TranscriptAlignment& o) const {
    bool comparison = (chr == o.chr) && (chrpos == o.chrpos) && (strand == o.strand) && (cigar.size() == o.cigar.size());
    if (!comparison) return false;
    for (int i = 0; i < cigar.size(); i++) {
      if (cigar[i] != o.cigar[i]) {
        return false;
      }
    }
    return true;
  }
};

namespace std {
  template<>
  struct hash<TranscriptAlignment> {
      inline size_t operator()(const TranscriptAlignment& x) const {
          uint32_t value = std::hash<int>()((x.chr+1)*x.chrpos + ((int) x.strand));
          for (uint32_t i = 0; i < x.cigar.size(); i++) {
            value ^= std::hash<uint32_t>()((i+1)*(x.cigar[i]+1));
          }
          return value;
      }
  };
  
  template<>
  struct hash<std::pair<TranscriptAlignment, TranscriptAlignment>> {
      inline size_t operator()(const std::pair<TranscriptAlignment, TranscriptAlignment>& x) const {
          return std::hash<TranscriptAlignment>()(x.first) ^ (std::hash<TranscriptAlignment>()(x.second)<<32);
      }
  };
}

struct GeneModel {
  int id;
  std::string name; // like ENSG
  std::string commonName; // like ALB
  std::vector<int> transcripts;  
  int chr;
  int start,stop; // start 0-based, 1 bp past end for stop
  bool strand;
};

struct Chromosome {
  int len;
  std::string name;
};

struct Transcriptome {
  std::vector<TranscriptModel> transcripts;
  std::vector<GeneModel> genes;
  std::vector<Chromosome> chr;

  
  std::unordered_map<std::string, int> chrNameToId;
  std::unordered_map<std::string, int> trNameToId;
  std::unordered_map<std::string, int> geneNameToId;
  
  // maps transcript tr and 0-based position trpos
  //      to chr:chrpos in genome mapping 
  //bool translateTrPosition(const std::string &tr, const int _trpos, std::string &chr, int& chrpos, std::string &gene_id);
  void loadChromosomes(const std::string &chrom_fn);
  bool translateTrPosition(const int tr, const int pos, const int length,  bool strand, TranscriptAlignment &aln) const;
  void parseGTF(const std::string &gtf_fn, const KmerIndex& index, const ProgramOptions& options, bool guessChromosomes);
  int addGTFLine(const std::string& line, const KmerIndex& index, bool guessChromosome);
  void loadTranscriptome(const KmerIndex& index, std::istream &in, const ProgramOptions& options);

};







#endif // GENE_MODEL_HPP