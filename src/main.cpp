#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include <getopt.h>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);


#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
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
	return !opt.files.empty();
}


void usage()
{
	// TODO: write me !
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
  if (!CheckOptions(opt)) {
    usage();
    exit(1);
  }
	Kmer::set_k(opt.k);
	std::cerr << "setting k = " << opt.k << std::endl;
	KmerIndex index(opt);
	ProcessReads<KmerIndex, MinCollector>(index, opt);
	
	return 0;
}
