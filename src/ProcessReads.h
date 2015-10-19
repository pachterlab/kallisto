#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "MinCollector.h"

#include "common.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  fp1(0),fp2(0),seq1(0),seq2(0),l1(0),l2(0),
  paired(!opt.single_end), files(opt.files),
  current_file(0), state(false) {}
  bool empty();
  ~SequenceReader();

  bool fetchSequences(char *buf, const int limit,  std::vector<std::pair<const char*, int>>& seqs);
  bool fetchSequencesFull(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals);

private:
  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  int l1,l2;
  bool paired;
  const std::vector<std::string>& files;
  int current_file;
  bool state; // is the file open
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc)
    : tc(tc), index(index), opt(opt), SR(opt), numreads(0), nummapped(0) {}

  std::mutex reader_lock;
  std::mutex writer_lock;

  SequenceReader SR;
  MinCollector& tc;
  KmerIndex& index;
  const ProgramOptions& opt;
  int numreads;
  int nummapped;
  void processReads();
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, int n);
};

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp);
  ~ReadProcessor();
  char *buffer;
  size_t bufsize;
  bool paired;
  const MinCollector& tc;
  const KmerIndex& index;
  MasterProcessor& mp;
  int numreads;

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<std::vector<int>> newEcs;

  std::vector<int> counts;

  void operator()();
  void clear();
};




#endif // KALLISTO_PROCESSREADS_H
