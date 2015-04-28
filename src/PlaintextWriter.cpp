#include "PlaintextWriter.h"

void plaintext_writer(
    const std::string& out_name,
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

  of << "target_id" << "\t"
    << "kallisto_id" << "\t"
    << "length" << "\t"
    << "eff_length" << "\t"
    << "est_counts" << "\t"
    << "tpm" << std::endl;

  of.close();
}
