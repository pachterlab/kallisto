#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);


#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"


using namespace std;


void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
	int version_flag = 0;

  const char* opt_string = "t:i:k:s:o:n:f:";
  static struct option long_options[] =
  {
		// long args
    {"verbose", no_argument, &verbose_flag, 1},
		{"version", no_argument, &version_flag, 1},
		// short args
    {"threads", required_argument, 0, 't'},
		{"index", required_argument, 0, 'i'},
    {"kmer-size", required_argument, 0, 'k'},
    {"seed", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
		{"iterations", required_argument, 0, 'n'},
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
		case 'i':
		{
			opt.index = optarg;
			break;
		}
		case 'n':
		{
			stringstream(optarg) >> opt.iterations;
			break;
		}
    case 'k':
		{
      stringstream(optarg) >> opt.k;
      break;
		}
    case 'o':
      opt.output = optarg;
      break;
    case 's':
			stringstream(optarg) >> opt.seed;
      break;
    case 't':
			stringstream(optarg) >> opt.threads;
      break;
		case 'f':
			opt.transfasta = optarg;
			break;
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
	if (version_flag) {
		opt.version = true;
	}
}


bool CheckOptions(ProgramOptions& opt) {

	bool ret = true;

	if (opt.k <= 0 || opt.k > Kmer::MAX_K) {
		cerr << "Error: invalid k-mer size " << opt.k << endl;
		ret = false;
	}
	
	if (!opt.transfasta.empty()) {
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

	} else {
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

		// check for read files
		if (opt.files.size() == 0) {
			cerr << "Error: Missing read files" << endl;
			ret = false;
		} else {
			struct stat stFileInfo;
			for (auto &fn : opt.files) {
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
				if (!S_ISDIR(stFileInfo.st_mode)) {
					cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
					ret = false;
				}
			}
		}

		if (opt.threads <= 0) {
			cerr << "Error: invalid number of threads " << opt.threads << endl;
			ret = false;
		} else {
			unsigned int n = std::thread::hardware_concurrency();
			if (n != 0 && n < opt.threads) {
				cerr << "Warning: you asked for " << opt.threads << ", but only " << n << " cores on the machine" << endl;
			}
		}
	}
	
	return ret;
}

void version()
{
  cout << "Kallisto, version: " << 	KALLISTO_VERSION << endl;
}

void usage()
{
	cout << "Kallisto " << endl
			 << "Does transcriptome stuff" << endl << endl
			 << "Usage: Kallisto [options]   FASTQ-file" << endl << endl
			 << "-k, --kmer-size=INT         Size of k-mers" << endl
			 << "-t, --threads=INT           Number of threads to use (default value 1)" << endl
			 << "-s, --seed=INT              Seed value for randomness (default value 0, use time based randomness)" << endl
			 << "-i, --index=INT             File with fragment length distribution " << endl
			 << "-n, --iterations=INT        Number of iterations of EM algorithm (default value 500)" << endl
			 << "-f, --trans-fasta=INT       FASTA file containing reference transcriptome " << endl
			 << "-o, --output-dir=INT        Directory to store output to" << endl
			 << "    --verbose               Print lots of messages during run" << endl 
			 << "    --version               Display version info" << endl << endl;
}


int main(int argc, char *argv[])
{
	ProgramOptions opt;
  ParseOptions(argc,argv,opt);
	if (opt.version) {
		version();
		exit(1);
	}
  if (!CheckOptions(opt)) {
		cout << endl;
    usage();
    exit(1);
  }
	Kmer::set_k(opt.k);
	std::cerr << "setting k = " << opt.k << std::endl;
	KmerIndex index(opt);
	ProcessReads<KmerIndex, MinCollector>(index, opt);
	
	return 0;
}
