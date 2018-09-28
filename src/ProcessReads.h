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
#include "PseudoBam.h"
#include "EMAlgorithm.h"
#include "GeneModel.h"


#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

class MasterProcessor;

int ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);
int ProcessBatchReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, std::vector<std::vector<int>> &batchCounts);
int findFirstMappingKmer(const std::vector<std::pair<KmerEntry,int>> &v,KmerEntry &val);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(!opt.single_end), files(opt.files),
  f_umi(new std::ifstream{}),
  current_file(0), state(false), readbatch_id(-1) {}
  SequenceReader() :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(false), 
  f_umi(new std::ifstream{}),
  current_file(0), state(false), readbatch_id(-1) {}
  SequenceReader(SequenceReader&& o);
  
  bool empty();
  void reset();
  ~SequenceReader();

  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<std::string>& umis, int &readbatch_id,
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
  int readbatch_id = -1;
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc, const Transcriptome& model)
    : tc(tc), index(index), model(model), bamfp(nullptr), bamfps(nullptr), bamh(nullptr), opt(opt), SR(opt), numreads(0)
    ,nummapped(0), num_umi(0), bufsize(1ULL<<23), tlencount(0), biasCount(0), maxBiasCount((opt.bias) ? 1000000 : 0), last_pseudobatch_id (-1) { 
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
      if (opt.pseudobam) {
        pseudobatchf_out.open(opt.output + "/pseudoaln.bin", std::ios::out | std::ios::binary);
      }

    }

  ~MasterProcessor() {
    if (bamfp) {
      hts_close(bamfp);
      bamfp = nullptr;
    }
    if (bamfps) {
      for (int i = 0; i < numSortFiles; i++) {
        if (bamfps[i]) {
          hts_close(bamfps[i]);
          bamfps[i] = nullptr;
        }
      } 
      bamfps = nullptr;
    }
    if (bamh) {
      bam_hdr_destroy(bamh);
      bamh = nullptr;
    }    
  }

  std::mutex reader_lock;
  std::mutex writer_lock;


  SequenceReader SR;
  MinCollector& tc;
  KmerIndex& index;
  const Transcriptome& model;
  htsFile *bamfp;
  const int numSortFiles = 32;
  htsFile **bamfps;

  bam_hdr_t *bamh;
  const ProgramOptions& opt;
  int numreads;
  int nummapped;
  int num_umi;
  size_t bufsize;
  std::atomic<int> tlencount;
  std::atomic<int> biasCount;
  std::vector<std::vector<int>> batchCounts;
  const int maxBiasCount;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> newECcount;
  std::ofstream ofusion;
  std::ofstream pseudobatchf_out;
  std::ifstream pseudobatchf_in;
  std::vector<PseudoAlignmentBatch> pseudobatch_stragglers;
  int last_pseudobatch_id;
  void outputFusion(const std::stringstream &o);
  std::vector<std::unordered_map<std::vector<int>, int, SortedVectorHasher>> newBatchECcount;
  std::vector<std::vector<std::pair<int, std::string>>> batchUmis;
  std::vector<std::vector<std::pair<std::vector<int>, std::string>>> newBatchECumis;
  void processReads();
  void processAln(const EMAlgorithm& em, bool useEM);
  void writePseudoBam(const std::vector<bam1_t> &bv);
  void writeSortedPseudobam(const std::vector<std::vector<bam1_t>> &bvv);
  std::vector<uint64_t> breakpoints;
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, int n, std::vector<int>& flens, std::vector<int> &bias, const PseudoAlignmentBatch& pseudobatch, int id = -1);



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
  PseudoAlignmentBatch pseudobatch;
  

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


class AlnProcessor {
public:
  AlnProcessor(const KmerIndex& index, const ProgramOptions& opt, MasterProcessor& mp, const EMAlgorithm& em, const Transcriptome& model, bool useEM, int id = -1);
  AlnProcessor(AlnProcessor && o);
  ~AlnProcessor();
  char *buffer;
  char *bambuffer;
  size_t bufsize;
  size_t bambufsize;
  bool paired;
  std::vector<std::pair<int, std::string>> ec_umi;
  const KmerIndex& index;
  const EMAlgorithm& em;
  MasterProcessor& mp;
  SequenceReader batchSR;
  int numreads;
  int id;
  PseudoAlignmentBatch pseudobatch;
  const Transcriptome& model;
  bool useEM;
  


  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<std::string> umis;

  void operator()();
  void processBufferTrans();
  void processBufferGenome();
  void clear();
};


int fillBamRecord(bam1_t &b, uint8_t* buf, const char *seq, const char *name, const char *qual, int slen, int nlen, bool unmapped, int auxlen);
void fixCigarStringTrans(bam1_t &b, int rlen, int softclip, int overhang);
void fixCigarStringGenome(bam1_t &b, const TranscriptAlignment& tra);
void reverseComplementSeqInData(bam1_t &b);

#endif // KALLISTO_PROCESSREADS_H
