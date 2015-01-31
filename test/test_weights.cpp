#include "catch.hpp"

#include <vector>

#include "weights.h"

TEST_CASE("effective lengths", "[eff_lens]")
{
    std::vector<int> lens(3, 19);
    auto res = calc_eff_lens(lens, 3.0);
    for (auto& val: res)
    {
        REQUIRE(val == static_cast<double>(19) - 3.0 + 1);
    }
}
