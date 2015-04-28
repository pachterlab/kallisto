#include "PlaintextWriter.h"

std::vector<double> counts_to_tpm(const std::vector<double>& est_counts,
    const std::vector<double>& eff_lens) {
  assert( est_counts.size() == eff_lens.size() );
  const double MILLION {1e6};

  std::vector<double> tpm(est_counts.size(), 0.0);

  // TODO: get rid of this thing after testing
  assert(tpm.size() == est_counts.size());

  double total_mass = 0.0;

  for (size_t i = 0; i < est_counts.size(); ++i) {
    if (eff_lens[i] < 1.0) {
      std::cerr << "Why is this eff_len < 1.0? id: " << i << std::endl;
    }
    tpm[i] = (est_counts[i] / eff_lens[i]);
    total_mass += tpm[i];
  }

  for (size_t i = 0; i < est_counts.size(); ++i) {
    tpm[i] = (tpm[i] * MILLION) / total_mass;
  }

  return tpm;
}

void plaintext_writer(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const std::vector<int>& lens
    ){

  // TODO: convert to TPM
  // writeout in format:
  std::ofstream of;
  of.open( out_name );

  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << out_name << std::endl;

    exit(1);
  }

  auto tpm = counts_to_tpm(alpha, eff_lens);

  of << "target_id" << "\t"
    << "kallisto_id" << "\t"
    << "length" << "\t"
    << "eff_length" << "\t"
    << "est_counts" << "\t"
    << "tpm" << std::endl;

  for (auto i = 0; i < alpha.size(); ++i) {
    of << targ_ids[i] << '\t'
      << i << '\t'
      << lens[i] << '\t'
      << eff_lens[i] << '\t'
      << alpha[i] << '\t'
      << tpm[i] << std::endl;
  }

  of.close();
}
