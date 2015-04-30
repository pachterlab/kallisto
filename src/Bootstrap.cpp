#include "Bootstrap.h"
#include "weights.h"
#include "EMAlgorithm.h"

EMAlgorithm Bootstrap::run_em(const EMAlgorithm& em_start) {
    auto counts = mult_.sample();
    auto weights = calc_weights(counts, ecmap_, eff_lens_);
    EMAlgorithm em(ecmap_, counts, target_names_, eff_lens_, weights);

    em.set_start(em_start);
    em.run(10000, 20);
    em.compute_rho();

    return em;
}
