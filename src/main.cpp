#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>

#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif


#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"
#include "EMAlgorithm.h"
#include "weights.h"


using namespace std;


void ParseOptionsIndex(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;

  const char *opt_string = "i:k:f:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    // short args
    {"index", required_argument, 0, 'i'},
    {"kmer-size", required_argument, 0, 'k'},
    {"trans-fasta", required_argument, 0, 'f'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    case 'f': {
      opt.transfasta = optarg;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
}


void ParseOptionsEM(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;

  const char *opt_string = "t:i:s:o:n:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    // short args
    {"threads", required_argument, 0, 't'},
    {"index", required_argument, 0, 'i'},
    {"skip", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.skip;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.iterations;
      break;
    }
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
}

void ParseOptionsEMOnly(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;

  const char *opt_string = "t:s:o:n:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    // short args
    {"threads", required_argument, 0, 't'},
    {"seed", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.iterations;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
}


bool CheckOptionsIndex(ProgramOptions& opt) {

  bool ret = true;

  if (opt.k <= 0 || opt.k >= Kmer::MAX_K) {
    cerr << "Error: invalid k-mer size " << opt.k << ", maximum is " << (Kmer::MAX_K -1) << endl;
    ret = false;
  }

  // we want to generate the index, check k, index and transfasta
  struct stat stFileInfo;
  auto intStat = stat(opt.transfasta.c_str(), &stFileInfo);
  if (intStat != 0) {
    cerr << "Error: transcript fasta file not found " << opt.transfasta << endl;
    ret = false;
  }

  if (opt.index.empty()) {
    cerr << "Error: need to specify index name" << endl;
    ret = false;
  }

  return ret;
}

bool CheckOptionsEM(ProgramOptions& opt, bool emonly = false) {

  bool ret = true;


  // check for index
  if (!emonly) {
    if (opt.index.empty()) {
      cerr << "Error: index file missing" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.index.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: index file not found " << opt.index << endl;
        ret = false;
      }
    }
  }

  // check for read files
  if (!emonly) {
    if (opt.files.size() == 0) {
      cerr << "Error: Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      for (auto& fn : opt.files) {
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << "Error: file not found " << fn << endl;
          ret = false;
        }
      }
    }

    if (!(opt.files.size() == 1 || opt.files.size() == 2)) {
      cerr << "Error: Input files should be 1 or 2 files only" << endl;
      ret = false;
    }

    if (opt.skip <= 0) {
      cerr << "Error: skip has to be a positive integer" << endl;
      ret = false;
    }
  }


  if (opt.iterations <= 0) {
    cerr << "Error: Invalid number of iterations " << opt.iterations << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      } else if (emonly) {
        // check for directory/counts.txt
        struct stat stCountInfo;
        auto intcountstat = stat((opt.output + "/counts.txt" ).c_str(), &stCountInfo);
        if (intcountstat != 0) {
          cerr << "Error: could not find file " << opt.output << "/counts.txt" << endl;
          ret = false;
        }

        // check for directory/index.saved
        struct stat stIndexInfo;
        auto intindexstat = stat((opt.output + "/index.saved").c_str(), &stIndexInfo);
        if (intindexstat != 0) {
          cerr << "Error: could not find index " << opt.output << "/index.saved" << endl;
          ret = false;
        }
        opt.index = (opt.output + "/index.saved");
      }
    } else {
      if (emonly) {
        cerr << "Error: output directory needs to exist, run em first" << endl;
        ret = false;
      } else {
        // create directory
        if (mkdir(opt.output.c_str(), 0777) == -1) {
          cerr << "Error: could not create directory " << opt.output << endl;
          ret = false;
        }
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
  }


  return ret;
}


void PrintCite() {
  cout << "The paper describing this software has not been published." << endl;
  //  cerr << "When using this program in your research, please cite" << endl << endl;
}

void PrintVersion() {
  cout << "Kallisto, version: " << 	KALLISTO_VERSION << endl;
}

void usage() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto CMD [options] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    index         Builds the index "<< endl
       << "    em            Runs the EM algorithm " << endl
       << "    em-only       Runs the EM algorithm without mapping" << endl
       << "    cite          Prints citation information " << endl
       << "    version       Prints version information"<< endl << endl;
}


void usageIndex() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto index [options]" << endl << endl
       << "-k, --kmer-size=INT         Size of k-mers, default (21), max value is " << (Kmer::MAX_K-1) << endl
       << "-i, --index=STRING             Filename for index to be constructed " << endl
       << "-f, --trans-fasta=STRING       FASTA file containing reference transcriptome " << endl
       << "    --verbose               Print lots of messages during run" << endl;
}

void usageEM() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto em [options] FASTQ-files" << endl << endl
       << "-t, --threads=INT           Number of threads to use (default value 1)" << endl
       << "-i, --index=INT             Filename for index " << endl
       << "-s, --seed=INT              Seed value for randomness (default value 0, use time based randomness)" << endl
       << "-n, --iterations=INT        Number of iterations of EM algorithm (default value 500)" << endl
       << "-o, --output-dir=STRING        Directory to store output to" << endl
       << "    --verbose               Print lots of messages during run" << endl;
}

void usageEMOnly() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto em-only [options]" << endl << endl
       << "-t, --threads=INT           Number of threads to use (default value 1)" << endl
       << "-s, --seed=INT              Seed value for randomness (default value 0, use time based randomness)" << endl
       << "-n, --iterations=INT        Number of iterations of EM algorithm (default value 500)" << endl
       << "-o, --output-dir=STRING        Directory to store output to" << endl
       << "    --verbose               Print lots of messages during run" << endl;
}


int main(int argc, char *argv[]) {

  if (argc < 2) {
    usage();
    exit(1);
  } else {
    ProgramOptions opt;
    string cmd(argv[1]);
    if (cmd == "version") {
      PrintVersion();
    } else if (cmd == "cite") {
      PrintCite();
    } else if (cmd == "index") {
      if (argc==2) {
        usageIndex();
        return 0;
      }
      ParseOptionsIndex(argc-1,argv+1,opt);
      if (!CheckOptionsIndex(opt)) {
        usageIndex();
        exit(1);
      } else {
        // create an index
        Kmer::set_k(opt.k);
        KmerIndex index(opt);
        std::cerr << "Building index from: " << opt.transfasta << std::endl;
        index.BuildTranscripts(opt.transfasta);
        index.write(opt.index);
      }
    } else if (cmd == "em") {
      if (argc==2) {
        usageEM();
        return 0;
      }
      ParseOptionsEM(argc-1,argv+1,opt);
      if (!CheckOptionsEM(opt)) {
        usageEM();
        exit(1);
      } else {
        // run the em algorithm
        KmerIndex index(opt);
        index.load(opt);
        auto collection = ProcessReads<KmerIndex, MinCollector<KmerIndex>>(index, opt);
        // save modified index for future use
        index.write((opt.output+"/index.saved"), false);
        // compute mean frag length somewhere?
        auto eff_lens = calc_eff_lens(index.trans_lens_, 30.0);
        auto weights = calc_weights (collection.counts, index.ecmap, eff_lens);
        EMAlgorithm<KmerIndex> em(opt, index, collection.counts, eff_lens, weights);
        em.run();
        em.compute_rho();
        em.write(opt.output);
      }
    } else if (cmd == "em-only") {
      if (argc==2) {
        usageEMOnly();
        return 0;
      }
      ParseOptionsEMOnly(argc-1,argv+1,opt);
      if (!CheckOptionsEM(opt, true)) {
        usageEMOnly();
        exit(1);
      } else {
        // run the em algorithm
        KmerIndex index(opt);
        index.load(opt, false); // skip the k-mer map
        MinCollector<KmerIndex> collection(index, opt);
        collection.loadCounts(opt);
        // compute mean frag length somewhere?
        auto eff_lens = calc_eff_lens(index.trans_lens_, 30.0);
        auto weights = calc_weights (collection.counts, index.ecmap, eff_lens);
        EMAlgorithm<KmerIndex> em(opt, index, collection.counts, eff_lens, weights);
        em.run();
        em.compute_rho();
        em.write(opt.output);
      }
    } else {
      cerr << "Did not understand command " << cmd << endl;
      usage();
      exit(1);
    }

  }
  return 0;
}
