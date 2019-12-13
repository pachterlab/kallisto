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
#include "KmerIndex.h"
#include "common.h"
#include "PseudoBam.h"
#include "EMAlgorithm.h"
#include "GeneModel.h"
#include "BUSData.h"
#include "BUSTools.h"
#include <htslib/sam.h>


#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

class MasterProcessor;

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);
int64_t ProcessBatchReads(MasterProcessor& MP, const ProgramOptions& opt);
int64_t ProcessBUSReads(MasterProcessor& MP, const ProgramOptions& opt);
int findFirstMappingKmer(const std::vector<std::pair<KmerEntry,int>> &v,KmerEntry &val);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  readbatch_id(-1) {};
  SequenceReader() : state(false), readbatch_id(-1) {};
  virtual ~SequenceReader() {}
//  SequenceReader(SequenceReader&& o);
  
  virtual bool empty() = 0;
  virtual void reset();
  virtual void reserveNfiles(int n) = 0;
  virtual bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false) = 0;


public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:

  FastqSequenceReader(const ProgramOptions& opt) : SequenceReader(opt),
  current_file(0), paired(!opt.single_end), files(opt.files),
  f_umi(new std::ifstream{}) {
    SequenceReader::state = false;

    if (opt.bus_mode) {
      nfiles = opt.busOptions.nfiles;      
    } else {
      nfiles = paired ? 2 : 1;
    }
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : SequenceReader(), 
  paired(false), 
  f_umi(new std::ifstream{}),
  current_file(0) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();

  bool empty();  
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false);

public:
  int nfiles = 1;
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  bool paired;
  std::vector<std::string> files;
  std::vector<std::string> umi_files;
  std::unique_ptr<std::ifstream> f_umi;
  int current_file;
  std::vector<kseq_t*> seq;
};

class BamSequenceReader : public SequenceReader {
public:

  BamSequenceReader(const ProgramOptions& opt) :
  SequenceReader(opt) {
    SequenceReader::state = true;

    fp = bgzf_open(opt.files[0].c_str(), "rb");
    head = bam_hdr_read(fp);
    rec = bam_init1();
    
    err = bam_read1(fp, rec);
    eseq = bam_get_seq(rec);
    l_seq = rec->core.l_qseq;

    bc = bam_aux2Z(bam_aux_get(rec, "CR"));
    l_bc = 0;
    for (char *pbc = bc; *pbc != '\0'; ++pbc) {
      ++l_bc;
    }

    umi = bam_aux2Z(bam_aux_get(rec, "UR"));
    l_umi = 0;
    for (char *pumi = umi; *pumi != '\0'; ++pumi) {
      ++l_umi;
    }
  }
  BamSequenceReader() : SequenceReader() {};
  BamSequenceReader(BamSequenceReader &&o);
  ~BamSequenceReader();
  
  bool empty();
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false);

public:
  BGZF *fp;
  bam_hdr_t *head;
  bam1_t *rec;
  uint8_t *eseq;
  int32_t l_seq;
  char *bc;
  int l_bc;
  char *umi;
  int l_umi;
  int err;
  
  
private:
  static const std::string seq_enc;
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc, const Transcriptome& model)
    : tc(tc), index(index), model(model), bamfp(nullptr), bamfps(nullptr), bamh(nullptr), opt(opt), numreads(0)
    ,nummapped(0), num_umi(0), bufsize(1ULL<<23), tlencount(0), biasCount(0), maxBiasCount((opt.bias) ? 1000000 : 0), last_pseudobatch_id (-1) { 
      if (opt.bam) {
        SR = new BamSequenceReader(opt);
      } else {
        SR = new FastqSequenceReader(opt);
      }

      if (opt.batch_mode) {
        memset(&bus_bc_len[0],0,33);
        memset(&bus_umi_len[0],0,33);
        batchCounts.assign(opt.batch_ids.size(), {});
        tlencounts.assign(opt.batch_ids.size(), 0);
        batchFlens.assign(opt.batch_ids.size(), std::vector<int>(1000,0));
        /*for (auto &t : batchCounts) {
          t.resize(tc.counts.size(),0);
        }*/
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
      if (opt.bus_mode) {
        busf_out.open(opt.output + "/output.bus", std::ios::out | std::ios::binary);
        
        writeBUSHeader(busf_out, opt.busOptions.getBCLength(), opt.busOptions.getUMILength());
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
    delete SR;
  }

  std::mutex reader_lock;
  std::mutex writer_lock;


  SequenceReader *SR;
  MinCollector& tc;
  KmerIndex& index;
  const Transcriptome& model;
  htsFile *bamfp;
  const int numSortFiles = 32;
  htsFile **bamfps;

  bam_hdr_t *bamh;
  const ProgramOptions& opt;
  int64_t numreads;
  int64_t nummapped;
  int64_t num_umi;
  size_t bufsize;

  int bus_bc_len[33];
  int bus_umi_len[33];

  std::atomic<int> tlencount;
  std::vector<int> tlencounts;
  std::atomic<int> biasCount;
  std::vector<std::vector<int>> batchFlens;
  std::vector<std::vector<std::pair<int32_t, int32_t>>> batchCounts;
  std::vector<std::vector<int32_t>> tmp_bc;
  const int maxBiasCount;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> newECcount;
  //  std::vector<std::pair<BUSData, std::vector<int32_t>>> newB;  
  EcMap bus_ecmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> bus_ecmapinv;


  std::ofstream ofusion;
  std::ofstream pseudobatchf_out;
  std::ofstream busf_out;
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
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, int n, std::vector<int>& flens, std::vector<int> &bias, const PseudoAlignmentBatch& pseudobatch, std::vector<BUSData> &bv, std::vector<std::pair<BUSData, std::vector<int32_t>>> newB, int *bc_len, int *umi_len,   int id = -1, int local_id = -1);  
};

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1, int local_id = -1);
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
  FastqSequenceReader batchSR;
  int64_t numreads;
  int id;
  int local_id;
  PseudoAlignmentBatch pseudobatch;
  

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  std::vector<std::string> umis;
  std::vector<std::vector<int>> newEcs;
  std::vector<int> flens;
  std::vector<int> bias5;

  std::vector<int> counts;

  void operator()();
  void processBuffer();
  void clear();
};


class BUSProcessor {
public:
  BUSProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1, int local_id = -1);
  BUSProcessor(BUSProcessor && o);
  ~BUSProcessor();
  char *buffer;
  
  size_t bufsize;
  bool paired;
  bool bam;
  bool num;
  const MinCollector& tc;
  const KmerIndex& index;
  MasterProcessor& mp;
  int64_t numreads;
  int id;
  int local_id;
  PseudoAlignmentBatch pseudobatch;

  int bc_len[33];
  int umi_len[33];
  

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;

  std::vector<std::vector<int>> newEcs;
  std::vector<int> flens;
  std::vector<int> bias5;
  std::vector<int> counts;
  std::vector<BUSData> bv;
  std::vector<std::pair<BUSData, std::vector<int32_t>>> newB;

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
  FastqSequenceReader batchSR;
  int64_t numreads;
  int id;
  PseudoAlignmentBatch pseudobatch;
  const Transcriptome& model;
  bool useEM;
  


  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
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
