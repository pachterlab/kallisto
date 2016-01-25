#include "PlaintextWriter.h"

std::vector<double> counts_to_tpm(const std::vector<double>& est_counts,
    const std::vector<double>& eff_lens) {
  assert( est_counts.size() == eff_lens.size() );
  const double MILLION {1e6};

  std::vector<double> tpm(est_counts.size(), 0.0);

  double total_mass = 0.0;

  for (size_t i = 0; i < est_counts.size(); ++i) {
    if (eff_lens[i] < 1.0) {
      std::cerr << "Why is this eff_len < 1.0? id: " << i << std::endl;
    }
    tpm[i] = (est_counts[i] / eff_lens[i]);
    total_mass += tpm[i];
  }

  for (size_t i = 0; i < est_counts.size(); ++i) {
    tpm[i] = (tpm[i] / total_mass) * MILLION;
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

  std::ofstream of;
  of.open( out_name );

  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << out_name << std::endl;

    exit(1);
  }

  auto tpm = counts_to_tpm(alpha, eff_lens);

  of << "target_id" << "\t"
    /* << "kallisto_id" << "\t" */
    << "length" << "\t"
    << "eff_length" << "\t"
    << "est_counts" << "\t"
    << "tpm" << std::endl;

  for (auto i = 0; i < alpha.size(); ++i) {
    of << targ_ids[i] << '\t'
      /* << i << '\t' */
      << lens[i] << '\t'
      << eff_lens[i] << '\t'
      << alpha[i] << '\t'
      << tpm[i] << std::endl;
  }

  of.close();
}

std::string to_json(const std::string& id, const std::string& val, bool quote,
    bool comma, int level) {
  std::string out;

  for (auto i = 0; i < level; ++i) {
    out += "\t";
  }

  out += '"';
  out += id;
  out += "\": ";
  if (quote) {
    out += '"';
  }
  out += val;
  if (quote) {
    out += '"';
  }
  if (comma) {
    out += ',';
  }

  return out;
}

void plaintext_aux(
    const std::string& out_name,
    const std::string& n_targs,
    const std::string& n_bootstrap,
    const std::string& n_processed,
    const std::string& version,
    const std::string& index_v,
    const std::string& start_time,
    const std::string& call) {
  std::ofstream of;
  of.open( out_name );

  of << "{" << std::endl <<
    to_json("n_targets", n_targs, false) << std::endl <<
    to_json("n_bootstraps", n_bootstrap, false) << std::endl <<
    to_json("n_processed", n_processed, false) << std::endl <<
    to_json("kallisto_version", version, true) << std::endl <<
    to_json("index_version", index_v, false) << std::endl <<
    to_json("start_time", start_time, true) << std::endl <<
    to_json("call", call, true, false) << std::endl <<
    "}" << std::endl;

  of.close();
}


void writeBatchMatrix(
  const std::string &prefix,
  const KmerIndex &index,
  const std::vector<std::string> &ids,
  std::vector<std::vector<int>> &counts) {

    std::string ecfilename = prefix + ".ec";
    std::string countsfilename = prefix + ".tsv";

    std::ofstream ecof, countsof;
    ecof.open(ecfilename.c_str(), std::ios::out);
    // output equivalence classes in the form "EC TXLIST";
    for (int i = 0; i < index.ecmap.size(); i++) {
      ecof << i << "\t";
      // output the rest of the class
      const auto &v = index.ecmap[i];
      bool first = true;
      for (auto x : v) {
        if (!first) {
          ecof << ",";
        } else {
          first = false;
        }
        ecof << x;
      }
      ecof << "\n";
    }
    ecof.close();

    countsof.open(countsfilename.c_str(), std::ios::out);
    for (int j = 0; j < ids.size(); j++) {
      countsof << "\t" << ids[j];
    }
    countsof << "\n";
    if (!counts.empty()) {
      // write out the NxM matrix, N is # of ecs, M is number of samples
      int M = counts.size();
      int N = 0;
      for (int j = 0; j < M; j++) {
        if (N < counts[j].size()) {
          N = counts[j].size();
        }
      }

      for (int i = 0; i < N; i++) {
        countsof << i;
        for (int j = 0; j < M; j++) {
          if (counts[j].size() <= i) {
            countsof << "\t0";
          } else {
            countsof << "\t" << counts[j][i];
          }
        }
        countsof << "\n";
      }
    }
    countsof.close();

}
