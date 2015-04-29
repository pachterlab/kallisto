#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.3"

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
  std::string transfasta;
  std::vector<std::string> files;
  bool plaintext;

ProgramOptions() :
  verbose(false),
  seed(42),
  threads(1),
  k(21),
  iterations(500),
  skip(1),
  min_range(2*k+1),
  fld(0.0),
  bootstrap(0),
  plaintext(false)
  {}
};

#endif // KALLISTO_COMMON_H
