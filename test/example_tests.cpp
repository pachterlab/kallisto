#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "example.h"

TEST_CASE("Square example", "[square]")
{
    REQUIRE( int(sq(3)) == 9 );
}
