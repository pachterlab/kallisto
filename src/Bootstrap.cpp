#include "Bootstrap.h"
#include "weights.h"
#include "EMAlgorithm.h"

EMAlgorithm Bootstrap::run_em() {
    auto counts = mult_.sample();
    auto weights = calc_weights(counts, ecmap_, eff_lens_);
    EMAlgorithm em(ecmap_, counts, target_names_, eff_lens_, weights);

    em.run();
    em.compute_rho();

    return em;
}
