#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H


#include <string>
#include <iostream>
#include <fstream>


template<typename Index, typename TranscriptCollector>
void ProcessReads(const Index& index, const ProgramOptions& opt) {
	
	// need to receive an index map
	std::ios_base::sync_with_stdio(false);

	gzFile fp = 0;
	kseq_t *seq = 0;
	std::vector<int> v;
	v.reserve(1000);

	int l; // length of read
	size_t nreads = 0;

	TranscriptCollector tc(opt);

	// for each file
	for (auto& file : opt.files) {
		fp = gzopen(file.c_str(), "r");
		seq = kseq_init(fp);
		// for each read
		while ((l = kseq_read(seq)) > 0) {
			nreads++;
			v.clear();
			// process read
			index.match(seq->seq.s, seq->seq.l, v);
			// iterate over all k-mers in read to find list of transcripts
			
			// collect the transcript information
			tc.collect(v);
		}
		gzclose(fp);
	}

	kseq_destroy(seq);

	// write output to outdir
	std::string outfile = opt.output + "/file"; // figure out filenaming scheme
	std::ofstream of;
	of.open(outfile.c_str(), std::ios::out);
	tc.write(of);
	of.close();
	
	return;
}


#endif // KALLISTO_PROCESSREADS_H
