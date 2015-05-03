#include "catch.hpp"

#include <iostream>
#include <vector>

#include "Multinomial.hpp"

TEST_CASE("multinomial", "[multinomial]")
{
    std::default_random_engine generator(42);
    std::vector<int> x {5, 5, 10, 0};
    std::discrete_distribution<int> dd(x.begin(), x.end());

    for (auto p : dd.probabilities()) {
        std::cout << p << std::endl;
    }

    Multinomial mult(x);
    for (auto i = 0; i < 1000; ++i) {
        auto samp = mult.sample();
        int cur_n {0};
        for (auto c : samp) {
            cur_n += c;
        }

        REQUIRE(cur_n == mult.n());
        REQUIRE(samp[3] == 0);
    }
}
