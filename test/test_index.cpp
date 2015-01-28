#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "common.h"
#include "KmerIndex.h"
#include "KmerIterator.hpp"

#include <string>
#include <stdio.h>


TEST_CASE("Build index", "[build_index]")
{
    ProgramOptions global_opts;
    std::string short_trans {"../test/input/10_trans_gt_500_bp.fasta"};

    global_opts.k = 21;
    KmerIndex kidx(global_opts);
    kidx.BuildTranscripts(short_trans);
}
