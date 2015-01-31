#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

#include <zlib.h>
#include "kseq.h"
#include <string>

#include <iostream>
#include <fstream>

#include "common.h"


template<typename Index, typename TranscriptCollector>
TranscriptCollector ProcessReads(Index& index, const ProgramOptions& opt) {
	
	// need to receive an index map
	std::ios_base::sync_with_stdio(false);


	bool paired = (opt.files.size() == 2);
	
	gzFile fp1 = 0, fp2 = 0;
	kseq_t *seq1 = 0, *seq2;
	std::vector<int> v;
	v.reserve(1000);

	int l1,l2; // length of read
	size_t nreads = 0;

	TranscriptCollector tc(index, opt);

	// for each file
	
	fp1 = gzopen(opt.files[0].c_str(), "r");
	seq1 = kseq_init(fp1);
	if (paired) {
		fp2 = gzopen(opt.files[1].c_str(),"r");
		seq2 = kseq_init(fp2);
	}


	// for each read
	while (true) {
		l1 = kseq_read(seq1);
		if (paired) {
			l2 = kseq_read(seq2);
		}
		if (l1 <= 0) {
			break;
		}
		if (paired && l2 <= 0) {
			break;
		}

		nreads++;
		v.clear();
		// process read
		index.match(seq1->seq.s, seq1->seq.l, v);
		if (paired) {
			index.match(seq2->seq.s, seq2->seq.l, v);
		}
		
		// collect the transcript information
		tc.collect(v);
		if (nreads % 100 == 0 ) {
			std::cerr << "Processed " << nreads << std::endl;
		}
	}
	gzclose(fp1);
	if (paired) {
		gzclose(fp2);
	}

	kseq_destroy(seq1);
	if (paired) {
		kseq_destroy(seq2);
	}

	// write output to outdir
	std::string outfile = opt.output + "/counts.txt"; // figure out filenaming scheme
	std::ofstream of;
	of.open(outfile.c_str(), std::ios::out);
	tc.write(of);
	of.close();
	
	return tc;
}


#endif // KALLISTO_PROCESSREADS_H
