#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>

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
#include "Inspect.h"
#include "Bootstrap.h"
#include "H5Writer.h"


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

void ParseOptionsInspect(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;

  const char *opt_string = "i:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    // short args
    {"index", required_argument, 0, 'i'},
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
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
}



void ParseOptionsEM(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;

  const char *opt_string = "t:i:s:l:o:n:m:d:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"index", required_argument, 0, 'i'},
    {"skip", required_argument, 0, 's'},
    {"fragment-length", required_argument, 0, 'l'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"min-range", required_argument, 0, 'm'},
    {"bootstrap-samples", required_argument, 0, 'b'},
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
    case 'l': {
      stringstream(optarg) >> opt.fld;
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
    case 'm': {
      stringstream(optarg) >> opt.min_range;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
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

  if (plaintext_flag) {
    opt.plaintext = true;
  }
}

void ParseOptionsEMOnly(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;

  const char *opt_string = "t:s:l:o:n:m:d:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"fragment-length", required_argument, 0, 'l'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"min-range", required_argument, 0, 'm'},
    {"bootstrap-samples", required_argument, 0, 'b'},
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
    case 'l': {
      stringstream(optarg) >> opt.fld;
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
    case 'm': {
      stringstream(optarg) >> opt.min_range;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }

  if (plaintext_flag) {
    opt.plaintext = true;
  }
}

bool ParseOptionsH5Dump(int argc, char **argv) {
  int peek_flag = 0;
  const char *opt_string = "";
  static struct option long_options[] = {
    // long args
    {"peek", no_argument, &peek_flag, 1},
    {0,0,0,0}
  };

  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    default: break;
    }
  }

  return static_cast<bool>(peek_flag);
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

  if (opt.fld == 0.0) {
    // In the future, if we have single-end data we should require this
    // argument
    cerr << "[quant] Mean fragment length not provided. Will estimate from data" << endl;
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.iterations <= 0) {
    cerr << "Error: Invalid number of iterations " << opt.iterations << endl;
    ret = false;
  }

  if (opt.min_range <= 0) {
    cerr << "Error: Invalid value for minimum range " << opt.min_range << endl;
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
        cerr << "Error: output directory needs to exist, run quant first" << endl;
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

  if (opt.bootstrap < 0) {
    cerr << "Error: number of bootstrap samples must be a non-negative integer." << endl;
    ret = false;
  }

  return ret;
}


bool CheckOptionsInspect(ProgramOptions& opt) {

  bool ret = true;
  // check for index
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
       << "Usage: kallisto CMD [options] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    index         Builds the index "<< endl
       << "    quant         Runs the quantification algorithm " << endl
       << "    h5dump        Get info from quantification and optionally output to plaintext " << endl
       << "    cite          Prints citation information " << endl
       << "    version       Prints version information"<< endl << endl;
}


void usageIndex() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto index [options]" << endl << endl
       << "-k, --kmer-size=INT         Size of k-mers, default (21), max value is " << (Kmer::MAX_K-1) << endl
       << "-i, --index=STRING             Filename for index to be constructed " << endl
       << "-f, --trans-fasta=STRING       FASTA file containing reference transcriptome " << endl;
}

void usageInspect() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto inspect [options]" << endl << endl
       << "-i, --index=STRING             Filename for index to inspect " << endl;
}

void usageEM() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto quant [options] FASTQ-files" << endl << endl
       << "-i, --index=INT               Filename for index " << endl
       << "-l, --fragment-length=DOUBLE  Estimated fragment length (default values are estimated from data)" << endl
       << "-m, --min-range               Minimum range to assign a read to a transcript (default value 2*k+1)" << endl
       << "-n, --iterations=INT          Number of iterations of EM algorithm (default value 500)" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples to perform (default value 0)" << endl
       << "--seed=INT                    Seed for bootstrap samples (default value 42)" << endl
       << "-o, --output-dir=STRING       Directory to store output to" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl;
}

void usageEMOnly() {
  cout << "Kallisto " << endl
       << "Does transcriptome stuff" << endl << endl
       << "Usage: Kallisto quant-only [options]" << endl << endl
       << "-s, --skip=INT                Number of k-mers to skip (default value 1)" << endl
       << "-l, --fragment-length=DOUBLE  Estimated fragment length (default values are estimated from data)" << endl
       << "-m, --min-range               Minimum range to assign a read to a transcript (default value 2*k+1)" << endl
       << "-n, --iterations=INT          Number of iterations of EM algorithm (default value 500)" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples to perform (default value 0)" << endl
       << "--seed=INT                    Seed for bootstrap samples (default value 42)" << endl
       << "-o, --output-dir=STRING       Directory to store output to" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl;
}

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    usage();
    exit(1);
  } else {
    auto start_time(get_local_time());
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
        index.BuildTranscripts(opt);
        index.write(opt.index);
      }
    } else if (cmd == "inspect") {
      if (argc==2) {
        usageInspect();
        return 0;
      }
      ParseOptionsInspect(argc-1, argv+1, opt);
      if (!CheckOptionsInspect(opt)) {
        usageInspect();
        exit(1);
      } else {
        KmerIndex index(opt);
        index.load(opt);
        InspectIndex(index);
      }

    } else if (cmd == "quant") {
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

        auto firstFile = opt.files[0];
        MinCollector collection(index, opt);
        if (firstFile.size() >= 4 && firstFile.compare(firstFile.size()-4,4,".bam") == 0) {
          ProcessBams<KmerIndex, MinCollector>(index, opt, collection);
        } else {
          ProcessReads<KmerIndex, MinCollector>(index, opt, collection);
        }
        // save modified index for future use
        index.write((opt.output+"/index.saved"), false);

        // if mean FL not provided, estimate
        auto mean_fl = (opt.fld > 0.0) ? opt.fld : get_mean_frag_len(collection);
        std::cerr << "Estimated mean fragment length: " << mean_fl << std::endl;
        auto eff_lens = calc_eff_lens(index.trans_lens_, mean_fl);
        auto weights = calc_weights (collection.counts, index.ecmap, eff_lens);
        EMAlgorithm em(index.ecmap, collection.counts, index.target_names_,
                       eff_lens, weights);
        em.run();

        std::string call = argv_to_string(argc, argv);

        H5Writer writer;
        if (!opt.plaintext) {
          writer.init(opt.output + "/abundance.h5", opt.bootstrap, 6,
              index.INDEX_VERSION, call, start_time);
          writer.write_main(em, index.target_names_, index.trans_lens_);
        } else {
          plaintext_aux(
              opt.output + "/run_info.json",
              std::string(std::to_string(eff_lens.size())),
              std::string(std::to_string(opt.bootstrap)),
              KALLISTO_VERSION,
              std::string(std::to_string(index.INDEX_VERSION)),
              start_time,
              call);

          plaintext_writer(opt.output + "/abundance.txt", em.target_names_,
              em.alpha_, em.eff_lens_, index.trans_lens_);
        }

        if (opt.bootstrap > 0) {
          std::cerr << "Bootstrapping!" << std::endl;
          auto B = opt.bootstrap;
          std::mt19937_64 rand;
          rand.seed( opt.seed );

          std::vector<size_t> seeds;
          for (auto s = 0; s < B; ++s) {
            seeds.push_back( rand() );
          }

          for (auto b = 0; b < B; ++b) {
            Bootstrap bs(collection.counts, index.ecmap,
                         index.target_names_, eff_lens, seeds[b]);
            std::cerr << "Running EM bootstrap: " << b << std::endl;
            auto res = bs.run_em(em);

            if (!opt.plaintext) {
              writer.write_bootstrap(res, b);
            } else {
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".txt",
                  em.target_names_, res.alpha_, em.eff_lens_, index.trans_lens_);
            }
          }

        }
      }
    } else if (cmd == "quant-only") {
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
        MinCollector collection(index, opt);
        collection.loadCounts(opt);
        // if mean FL not provided, estimate
        auto mean_fl = (opt.fld > 0.0) ? opt.fld : get_mean_frag_len(collection);
        std::cerr << "Estimated mean fragment length: " << mean_fl << std::endl;
        auto eff_lens = calc_eff_lens(index.trans_lens_, mean_fl);
        auto weights = calc_weights (collection.counts, index.ecmap, eff_lens);

        EMAlgorithm em(index.ecmap, collection.counts, index.target_names_,
                       eff_lens, weights);

        em.run();

        std::string call = argv_to_string(argc, argv);
        H5Writer writer;

        if (!opt.plaintext) {
          writer.init(opt.output + "/abundance.h5", opt.bootstrap, 6,
              index.INDEX_VERSION, call, start_time);
          writer.write_main(em, index.target_names_, index.trans_lens_);
        } else {
          plaintext_aux(
              opt.output + "/run_info.json",
              std::string(std::to_string(eff_lens.size())),
              std::string(std::to_string(opt.bootstrap)),
              KALLISTO_VERSION,
              std::string(std::to_string(index.INDEX_VERSION)),
              start_time,
              call);

          plaintext_writer(opt.output + "/abundance.txt", em.target_names_,
              em.alpha_, em.eff_lens_, index.trans_lens_);
        }

        if (opt.bootstrap > 0) {
          std::cerr << "Bootstrapping!" << std::endl;
          auto B = opt.bootstrap;
          std::mt19937_64 rand;
          rand.seed( opt.seed );

          std::vector<size_t> seeds;
          for (auto s = 0; s < B; ++s) {
            seeds.push_back( rand() );
          }

          for (auto b = 0; b < B; ++b) {
            Bootstrap bs(collection.counts, index.ecmap,
                         index.target_names_, eff_lens, seeds[b]);
            std::cerr << "Running EM bootstrap: " << b << std::endl;
            auto res = bs.run_em(em);

            if (!opt.plaintext) {
              writer.write_bootstrap(res, b);
            } else {
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".txt",
                  em.target_names_, res.alpha_, em.eff_lens_, index.trans_lens_);
            }
          }
        }
      }
    } else if (cmd == "h5dump") {

      if (argc != 4) {
        cerr << "Usage:\n\tkallisto h5dump /path/to/abundance.h5 OUTPUT_DIR" << endl;
        cerr << "\t\tOR\n"
          << "\tkallisto h5dump --peek /path/to/abundance.h5" << endl;
        exit(1);
      }

      auto peek = ParseOptionsH5Dump(argc-1,argv+1);

      std::string h5file;
      std::string out_dir;

      if (!peek) {
        h5file = argv[2];
        out_dir = argv[3];

        struct stat stFileInfo;
        auto intStat = stat(out_dir.c_str(), &stFileInfo);
        if (intStat == 0) {
          // file/dir exits
          if (!S_ISDIR(stFileInfo.st_mode)) {
            cerr << "Error: Tried to open " << out_dir << " but another file already exists there" << endl;
            exit(1);
          }
        } else if (mkdir(out_dir.c_str(), 0777) == -1) {
          cerr << "Error: could not create directory " << out_dir << endl;
          exit(1);
        }
      } else {
        h5file = argv[3];
        out_dir = "";
      }

      struct stat stFileInfo;
      auto intStat = stat(h5file.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: couldn't find file " << h5file << endl;
        exit(1);
      }

      H5Converter h5conv(h5file, out_dir);
      if (!peek) {
        h5conv.write_aux();
        h5conv.convert();
      }
    }  else {
      cerr << "Did not understand command " << cmd << endl;
      usage();
      exit(1);
    }

  }
  return 0;
}
