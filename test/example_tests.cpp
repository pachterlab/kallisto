#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "example.h"
#include "sequences.h"

TEST_CASE("Square example", "[square]")
{
    REQUIRE( int(sq(3)) == 9 );
}

TEST_CASE("Test kmers", "[kmer]")
{
    std::string seq1 {"ACGGT"};
    KmerView kmers1 (&seq1, 3);
    REQUIRE( !kmers1.done() );
    REQUIRE( kmers1.next() == "ACG" );
    REQUIRE( !kmers1.done() );
    REQUIRE( kmers1.next() == "CGG" );
    REQUIRE( !kmers1.done() );
    REQUIRE( kmers1.next() == "GGT" );
    REQUIRE( kmers1.done() );

    std::string seq2 {"ACG"};
    KmerView kmers2 (&seq2, 4);
    REQUIRE( kmers2.done() );

    std::string seq3 {"ACG"};
    KmerView kmers3 (&seq3, 3);
    REQUIRE( !kmers3.done() );
    REQUIRE( kmers3.next() == "ACG" );
    REQUIRE( kmers3.done() );
}
