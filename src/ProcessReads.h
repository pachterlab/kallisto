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
#include <atomic>
#include <condition_variable>

#include "MinCollector.h"

#include "common.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc);
int ProcessBatchReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, std::vector<std::vector<int>> &batchCounts);
int findFirstMappingKmer(const std::vector<std::pair<KmerEntry,int>> &v,KmerEntry &val);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(!opt.single_end), files(opt.files),
  f_umi(new std::ifstream{}),
  current_file(0), state(false) {}
  SequenceReader() :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(false), 
  f_umi(new std::ifstream{}),
  current_file(0), state(false) {}
  SequenceReader(SequenceReader&& o);
  
  bool empty();
  ~SequenceReader();

  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<std::string>& umis, 
                      bool full=false);

public:
  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  int l1,l2,nl1,nl2;
  bool paired;
  std::vector<std::string> files;
  std::vector<std::string> umi_files;
  std::unique_ptr<std::ifstream> f_umi;
  int current_file;
  bool state; // is the file open
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc)
    : tc(tc), index(index), opt(opt), SR(opt), numreads(0)
    ,nummapped(0), num_umi(0), tlencount(0), biasCount(0), maxBiasCount((opt.bias) ? 1000000 : 0) { 
      if (opt.batch_mode) {
        batchCounts.resize(opt.batch_ids.size(), {});
        
        for (auto &t : batchCounts) {
          t.resize(tc.counts.size(),0);
        }
        newBatchECcount.resize(opt.batch_ids.size());
        newBatchECumis.resize(opt.batch_ids.size());
        batchUmis.resize(opt.batch_ids.size());
      }
      if (opt.fusion) {
        ofusion.open(opt.output + "/fusion.txt");
        ofusion << "TYPE\tNAME1\tSEQ1\tKPOS1\tNAME2\tSEQ2\tKPOS2\tINFO\tPOS1\tPOS2\n";
      }

    }

  std::mutex reader_lock;
  std::mutex writer_lock;

  SequenceReader SR;
  MinCollector& tc;
  KmerIndex& index;
  const ProgramOptions& opt;
  int numreads;
  int nummapped;
  int num_umi;
  std::atomic<int> tlencount;
  std::atomic<int> biasCount;
  std::vector<std::vector<int>> batchCounts;
  const int maxBiasCount;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> newECcount;
  std::ofstream ofusion;
  void outputFusion(const std::stringstream &o);
  std::vector<std::unordered_map<std::vector<int>, int, SortedVectorHasher>> newBatchECcount;
  std::vector<std::vector<std::pair<int, std::string>>> batchUmis;
  std::vector<std::vector<std::pair<std::vector<int>, std::string>>> newBatchECumis;
  void processReads();

  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, int n, std::vector<int>& flens, std::vector<int> &bias, int id = -1);
};

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  size_t bufsize;
  bool paired;
  const MinCollector& tc;
  std::vector<std::pair<int, std::string>> ec_umi;
  std::vector<std::pair<std::vector<int>, std::string>> new_ec_umi;
  const KmerIndex& index;
  MasterProcessor& mp;
  SequenceReader batchSR;
  int numreads;
  int id;

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<std::string> umis;
  std::vector<std::vector<int>> newEcs;
  std::vector<int> flens;
  std::vector<int> bias5;

  std::vector<int> counts;

  void operator()();
  void processBuffer();
  void clear();
};




#endif // KALLISTO_PROCESSREADS_H
