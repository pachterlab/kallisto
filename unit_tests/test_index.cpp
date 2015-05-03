#include "catch.hpp"

#include "common.h"
#include "KmerIndex.h"
#include "KmerIterator.hpp"

#include <string>

#include <stdio.h>


// TEST_CASE("Build index", "[build_index]")
// {
//
// 	ProgramOptions global_opts, local_opts;
//     std::string short_trans {"../test/input/10_trans_gt_500_bp.fasta"};
//
//     global_opts.k = 21;
//
//     // TODO: figure out a better way to deal with initialization
//     Kmer::set_k(global_opts.k);
//
//     KmerIndex kidx(global_opts);
//     kidx.BuildTranscripts(short_trans);
//     kidx.write("tmp.idx");
//
//
//
// 		local_opts.index = "tmp.idx";
// 		local_opts.k = 0;
//     KmerIndex jidx(local_opts);
//     jidx.load(local_opts);
//
//
//     remove("tmp.idx");
//
// 		REQUIRE( kidx.k == local_opts.k );
//     REQUIRE( kidx.k == jidx.k );
//     REQUIRE( kidx.num_trans == jidx.num_trans );
//     REQUIRE( kidx.kmap.size() == jidx.kmap.size() );
//     REQUIRE( kidx.ecmap.size() == jidx.ecmap.size() );
//
//     // TODO: write tests to compare actual maps
// }
