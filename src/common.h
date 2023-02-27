#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.49.0"

// NOTE: MAKE SURE THIS FILE GETS INCLUDED FIRST IN ALL OTHER FILES AND BEFORE ANY EXTERNAL LIBRARIES

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "kseq.h"

#ifdef _WIN64
typedef unsigned int uint;
#endif

extern std::unordered_map<char, std::string>cfc_aa_map;
extern std::unordered_map<std::string, std::string>cfc_map;
extern std::string revcomp(const std::string s);

struct BUSOptionSubstr {
  BUSOptionSubstr() : fileno(-1), start(0), stop(0) {}
  BUSOptionSubstr(int f, int a, int b) : fileno(f), start(a), stop(b) {}
  int fileno;
  int start;
  int stop;
};

struct BUSOptions {
  int nfiles;

  std::vector<BUSOptionSubstr> umi;
  std::vector<BUSOptionSubstr> bc;
  std::vector<BUSOptionSubstr> seq;

  bool paired;
  bool aa;

  int getBCLength() const {
    int r =0 ;
    if (!bc.empty()) {
      for (auto& b : bc) {
        if (b.start < 0) {
          return 0;
        } else if (b.stop == 0) {
          return 0;
        } else {
          r += b.stop - b.start;
        }
      }
    }
    return r;
  }

  int getUMILength() const {
    int r =0 ;
    if (!umi.empty()) {
      for (auto& u : umi) {
        if (u.start < 0) {
          return 0;
        } else if (u.stop == 0) {
          return 0;
        } else {
          r += u.stop - u.start;
        }
      }
    }
    return r;
  }
};

struct ProgramOptions {
  bool verbose;
  bool aa;
  int threads;
  std::string index;
  int k;
  int g;
  int max_ec_size;
  int iterations;
  std::string output;
  int skip;
  size_t seed;
  double fld;
  double sd;
  int min_range;
  int bootstrap;
  int max_num_reads;
  std::vector<std::string> transfasta;
  bool batch_mode;
  bool pseudo_read_files_supplied;
  bool bus_mode;
  BUSOptions busOptions;
  bool pseudo_quant;
  bool bam;
  bool num;
  std::string batch_file_name;
  std::vector<std::vector<std::string> > batch_files;
  std::vector<std::string> batch_ids;
  std::vector<std::string> files;
  std::vector<std::string> umi_files;
  std::vector<std::string> d_list;
  bool plaintext;
  bool write_index;
  bool single_end;
  bool strand_specific;
  bool peek; // only used for H5Dump
  bool bias;
  bool pseudobam;
  bool genomebam;
  bool make_unique;
  bool fusion;
  enum class StrandType {None, FR, RF};
  StrandType strand;
  bool umi;
  bool batch_bus_write;
  bool batch_bus;
  std::string gfa; // used for inspect
  bool inspect_thorough;
  bool single_overhang;
  bool record_batch_bus_barcode;
  int input_interleaved_nfiles;
  std::string gtfFile;
  std::string chromFile;
  std::string bedFile;
  std::string technology;
  std::string tagsequence;
  std::string tccFile;
  std::string ecFile;
  std::string fldFile;
  std::string transcriptsFile;
  std::string genemap;

ProgramOptions() :
  verbose(false),
  threads(1),
  k(31),
  g(0),
  max_ec_size(0),
  iterations(500),
  skip(1),
  seed(42),
  fld(0.0),
  sd(0.0),
  min_range(1),
  bootstrap(0),
  max_num_reads(0),
  input_interleaved_nfiles(0),
  record_batch_bus_barcode(false),
  batch_mode(false),
  pseudo_read_files_supplied(false),
  bus_mode(false),
  pseudo_quant(false),
  bam(false),
  num(false),
  plaintext(false),
  write_index(false),
  single_end(false),
  strand_specific(false),
  peek(false),
  bias(false),
  pseudobam(false),
  genomebam(false),
  make_unique(false),
  fusion(false),
  strand(StrandType::None),
  umi(false),
  batch_bus_write(false),
  batch_bus(false),
  inspect_thorough(false),
  single_overhang(false),
  aa(false)
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(int64_t num);
std::string pretty_num(int num);


#ifdef KALLISTO_USE_ZLIB_NG
#include "zlib-ng/zlib.h"
#else
#include <zlib.h>
#endif

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif


#endif // KALLISTO_COMMON_H
