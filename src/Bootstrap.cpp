#include "Bootstrap.h"
#include "weights.h"
#include "EMAlgorithm.h"

EMAlgorithm Bootstrap::run_em(const EMAlgorithm& em_start) {
    auto counts = mult_.sample();
    EMAlgorithm em(counts, index_, tc_, mean_fl);

    //em.set_start(em_start);
    em.run(10000, 50, false, false);
    /* em.compute_rho(); */

    return em;
}
