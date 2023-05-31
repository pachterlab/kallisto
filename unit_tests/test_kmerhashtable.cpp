#include "catch.hpp"

#include <random>
#include <string>
#include <map>
#include <unordered_map>
#include <stdio.h>
#include <iostream>

#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "KmerHashTable.h"

#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif


using namespace std;

TEST_CASE("Build table", "[build_table]")
{

    ProgramOptions global_opts;
    std::string short_reads {"../test/input/short_reads.fastq"};

		std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(1, 100);
		
    // TODO: figure out a better way to deal with initialization
    Kmer::set_k(global_opts.k);

		unordered_map<Kmer, int, KmerHash> kmap1;
		map<Kmer, int> kmap2;
		KmerHashTable<int, KmerHash> kmap3;
		vector<Kmer> v;


		std::ios_base::sync_with_stdio(false);
		gzFile fp = gzopen(short_reads.c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		int val = 0;
		while ((l = kseq_read(seq)) >= 0) {
			KmerIterator kit(seq->seq.s), kit_end;
			for (; kit != kit_end; ++kit) {
				Kmer rep = kit->first.rep();
				//int val = dis(gen);
				v.push_back(rep);
				kmap1.insert({rep,val});
				kmap2.insert({rep,val});
				kmap3.insert({rep,val});
				val++;
			}
			
		}

		gzclose(fp);
		kseq_destroy(seq);

		REQUIRE(kmap1.size() == kmap2.size());
		REQUIRE(kmap3.size() == kmap1.size());
		cerr << kmap1.size() << endl;
		for (auto & km : v) {
			auto s1 = kmap1.find(km);
			REQUIRE(s1 != kmap1.end());
			REQUIRE(s1->first == km);
			//REQUIRE(s1->second == i);


			auto s2 = kmap2.find(km);
			REQUIRE(s2 != kmap2.end());
			REQUIRE(s2->first == km);
			REQUIRE(s2->second == s1->second);

			auto s3 = kmap3.find(km);
			REQUIRE(s3 != kmap3.end());
			REQUIRE(s3->second == s1->second);
		}
		
}
