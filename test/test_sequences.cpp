#include "catch.hpp"

#include "sequences.h"

TEST_CASE("Test kmers", "[kmer]")
{
    std::string seq1 {"ACGGT"};
    KmerView kmers1 (&seq1, 3);
    std::string res;
    REQUIRE( kmers1.next(res) );
    REQUIRE( res == "ACG" );
    REQUIRE( kmers1.next(res) );
    REQUIRE( res == "CGG" );
    REQUIRE( kmers1.next(res) );
    REQUIRE( res == "GGT" );
    REQUIRE( kmers1.done() );

    std::string seq2 {"ACG"};
    KmerView kmers2 (&seq2, 4);
    REQUIRE( kmers2.done() );
    REQUIRE( not kmers2.next(res) );

    std::string seq3 {"ACG"};
    KmerView kmers3 (&seq3, 3);
    REQUIRE( !kmers3.done() );
    REQUIRE( kmers3.next(res) );
    REQUIRE( res == "ACG" );
    REQUIRE( kmers3.done() );
}

TEST_CASE("Test reverse complement", "[reverse_complement]")
{
    std::string seq1 {"ACGGT"};
    REQUIRE( reverse_complement(seq1) == "ACCGT" );
    // TODO: add more tests
}
