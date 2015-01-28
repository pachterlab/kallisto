#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "common.h"
#include "KmerIndex.h"
#include "KmerIterator.hpp"

#include <string>
#include <stdio.h>

// #include "common.h"
// #include "KmerIndex.h"


TEST_CASE("Build index", "[build_index]")
{
    ProgramOptions global_opts;
    std::string short_trans {"../test/input/10_trans_gt_500_bp.fasta"};

    TestStruct ts(global_opts);
    ts.read_trans(short_trans);

    global_opts.k = 21;
    KmerIndex kidx(global_opts);
    kidx.BuildTranscripts(short_trans);
}

TEST_CASE("kmer iterator", "[kmer_iter]")
{
    std::cout << "******************************" << std::endl;
    std::cout << "testing kmer iterator" << std::endl;
    std::cout << "******************************" << std::endl;
    KmerIterator kit("aacgacgacgacgacgacgacgacgacgacgacgacgacgcgt"), kit_end;
    for(; kit != kit_end; ++kit) {
        std::cout << "kmer\t";
    }
    std::cout << std::endl;
}
