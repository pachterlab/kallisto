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


// TEST_CASE("calc weights", "[weights]")
// {
//     EcMap tmp_map;
//     tmp_map.insert( {0, {0}} );
//     tmp_map.insert( {1, {1}} );
//     tmp_map.insert( {2, {2}} );
//     tmp_map.insert( {3, {0, 2}} );
//     tmp_map.insert( {4, {0, 1}} );
//
//     std::vector<int> counts {3, 1, 0, 10, 7};
//     std::vector<double> eff_lens {307.4, 500.4, 302.0};
//
//     auto w = calc_weights(counts, tmp_map, eff_lens);
//     REQUIRE( w.size() == 5 );
//
//     std::cout << std::endl;
//
//     // TODO: actually compute the assertions
//     for (auto& kv : w) {
//         std::cout << "ecid:\t" << kv.first;
//         for (auto& a_weight : kv.second) {
//             std::cout << "\t" << a_weight;
//         }
//         std::cout << std::endl;
//     }
//
//     // TODO: more rigorous test
// }
