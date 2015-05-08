#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.42"

#include <string>
#include <vector>

struct ProgramOptions {
  bool verbose;
  int threads;
  std::string index;
  int k;
  int iterations;
  std::string output;
  int skip;
  size_t seed;
  double fld;
  int min_range;
  int bootstrap;
  std::vector<std::string> transfasta;
  std::vector<std::string> files;
  bool plaintext;
  bool write_index;
  bool single_end;
  bool peek; // only used for H5Dump

ProgramOptions() :
  verbose(false),
  threads(1),
  k(31),
  iterations(500),
  skip(1),
  seed(42),
  fld(0.0),
  min_range(1),
  bootstrap(0),
  plaintext(false),
  write_index(false),
  single_end(false),
  peek(false)
  {}
};

#endif // KALLISTO_COMMON_H
