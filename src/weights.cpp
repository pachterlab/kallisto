#include "weights.h"

#include <cmath>

const double MIN_ALPHA = 1e-8;

std::vector<double> get_frag_len_means(const std::vector<int>& lengths,
    const std::vector<double>& mean_frag_len_trunc) {

  std::vector<double> frag_len_means;
  frag_len_means.reserve( lengths.size() );

  // we'll just assume we don't see any fragments longer or equal to
  // MAX_FRAG_LEN
  double marginal_mean = mean_frag_len_trunc[MAX_FRAG_LEN - 1];

  for (size_t i = 0; i < lengths.size(); ++i) {

    if (lengths[i] >= MAX_FRAG_LEN) {
      frag_len_means.push_back(marginal_mean);
    } else {
      frag_len_means.push_back(mean_frag_len_trunc[ lengths[i] ]);
    }

  }

  return frag_len_means;
}

// XXX: DEPRECATED
std::vector<double> calc_eff_lens(const std::vector<int>& lengths, double mean)
{
  // for now do the total naive thing and subtract mean frag length
  std::vector<double> eff_lens;
  eff_lens.reserve(lengths.size());

  auto n_too_short = 0;

  for (auto& cur_len: lengths) {
      double cur_len_dbl = static_cast<double>(cur_len);
      double cur_eff_len = cur_len_dbl - mean + 1.0;
      if (cur_eff_len < 1.0) {
          // cur_eff_len = 1.0;
          cur_eff_len = cur_len_dbl;
          ++n_too_short;
      }
      eff_lens.push_back( cur_eff_len );
  }

  if (n_too_short > 0) {
    //std::cerr << "[quant] total targets with effective length less than " << mean <<
    //      ":\t" << n_too_short  << std::endl;
  }

  return eff_lens;
}

std::vector<double> calc_eff_lens(const std::vector<int>& lengths,
    const std::vector<double>& means) {
  std::vector<double> eff_lens;
  eff_lens.reserve( lengths.size() );

  assert( lengths.size() == means.size() );

  // Sir Too $hort comin' straight from Oakland
  auto n_too_short = 0;

  for (size_t i = 0; i < lengths.size(); ++i) {
    double cur_len_dbl = static_cast<double>(lengths[i]);
    double cur_eff_len = cur_len_dbl - means[i] + 1;
    if (cur_eff_len < 1.0) {
      cur_eff_len = cur_len_dbl;
      ++n_too_short;
    }
    eff_lens.push_back( cur_eff_len );
  }

  return eff_lens;
}

inline int update_hexamer(int hex, char c, bool revcomp) {
  if (!revcomp) {
    hex = ((hex & 0x3FF) << 2);
    switch (c & 0xDF) {
    case 'C': hex += 1; break;
    case 'G': hex += 2; break;
    case 'T': hex += 3; break;
    }
  } else {
    hex = (hex >> 2);
    switch (c & 0xDF) {
    case 'A': hex += 3072; break;
    case 'C': hex += 2048; break;
    case 'G': hex += 1024; break;
    }
  }
  return hex;
}

std::vector<double> update_eff_lens(
    const std::vector<double>& means,
    const MinCollector& tc,
    const KmerIndex &index,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    std::vector<double>& dbias5,
    const ProgramOptions& opt
    ) {

  double biasDataNorm = 0.0;
  double biasAlphaNorm = 0.0;
  const int num6mers = 4096;
  for (int i = 0; i < num6mers; i++) {
    biasDataNorm += tc.bias5[i];
  }

  dbias5.clear();
  dbias5.resize(num6mers, 0.0); // clear the bias

  index.loadTranscriptSequences();
  

  for (int i = 0; i < index.num_trans; i++) {
    if (index.target_lens_[i] < means[i]) {
      // this should never happen.. but I'll sleep better at night with this
      // condition -HP
      continue;
    }

    if (alpha[i] < MIN_ALPHA) {
      continue;
    }

    double contrib = 0.5*alpha[i]/eff_lens[i];
    if (opt.strand_specific) {
      contrib = alpha[i]/eff_lens[i];
    }
    int seqlen = index.target_seqs_[i].size();
    const char* cs = index.target_seqs_[i].c_str();

    if (!opt.strand_specific || (opt.strand == ProgramOptions::StrandType::FR)) {
      int hex = hexamerToInt(cs,false);
      int fwlimit = (int) std::max(seqlen - means[i] - 6, 0.0);
      for (int j = 0; j < fwlimit; j++) {
        dbias5[hex] += contrib;
        hex = update_hexamer(hex,*(cs+j+6),false);
      } 
    }

    if (!opt.strand_specific || (opt.strand == ProgramOptions::StrandType::RF)) {
      int bwlimit = (int) std::max(means[i] - 6, 0.0);
      int hex = hexamerToInt(cs+bwlimit,true);
      for (int j = bwlimit; j < seqlen - 6; j++) {
        dbias5[hex] += contrib;
        if (j < seqlen - 6) {
          hex = update_hexamer(hex,*(cs+j+6),true);
        }
      }
    }
  }

  for (int i = 0; i < num6mers; i++) {
    biasAlphaNorm += dbias5[i];
  }

  std::vector<double> biaslens(index.num_trans);

  for (int i = 0; i < index.num_trans; i++) {
    double efflen = 0.0;
    if (index.target_lens_[i] >= means[i] && alpha[i] >= MIN_ALPHA) {

      int seqlen = index.target_seqs_[i].size();
      const char* cs = index.target_seqs_[i].c_str();

      // forward direction
      if (!opt.strand_specific || (opt.strand == ProgramOptions::StrandType::FR)) {
        int hex = hexamerToInt(cs,false);
        int fwlimit = (int) std::max(seqlen - means[i] - 6, 0.0);
        for (int j = 0; j < fwlimit; j++) {
          //int hex = hexamerToInt(cs+j,false);
          //efflen += 0.5*(tc.bias5[hex]/biasDataNorm) / (dbias5[hex]/biasAlphaNorm );
          efflen += tc.bias5[hex] / dbias5[hex];
          hex = update_hexamer(hex,*(cs+j+6),false);
        }
      }
      if (!opt.strand_specific || (opt.strand == ProgramOptions::StrandType::RF)) {
        int bwlimit = (int) std::max(means[i] - 6 , 0.0);
        int hex = hexamerToInt(cs+bwlimit,true);
        for (int j = bwlimit; j < seqlen - 6; j++) {
          efflen += tc.bias5[hex] / dbias5[hex];
          if (j < seqlen-6) {
            hex = update_hexamer(hex,*(cs+j+6),true);
          }
        }
      }
      
      
      if (!opt.strand_specific) {
        efflen *= 0.5*biasAlphaNorm/biasDataNorm;
      } else {
        efflen *= biasAlphaNorm/biasDataNorm;
      }
    }


    if (efflen > means[i]) {
      //efflen *= ((seqlen-mean) / ((double) (seqlen-mean-6));
      biaslens[i] = efflen;
    } else {
      biaslens[i] = eff_lens[i]; // just for unexpressed sequences
    }
    //std::cout << index.target_names_[i] << "\t" << eff_lens[i] << "\t" << biaslens[i] << "\t" << efflen << "\n";
  }



  return biaslens;
}

WeightMap calc_weights(
  const std::vector<int>& counts,
  const EcMap& ecmap,
  const std::vector<double>& eff_lens)
{

  // TODO: throw some assertions in here to make sure the length of counts
  // and ec map are correct... as well as eff_lens size is reasonable

  // weights are stored _exactly_ in the same orientation as the ec map
  WeightMap weights(ecmap.size());

  for (size_t ec = 0; ec < ecmap.size(); ec++) {
    auto& v = ecmap[ec];
    //std::cout << ec;
    std::vector<double> trans_weights;
    trans_weights.reserve(v.size());

    for (auto& trans_id : v) {
      trans_weights.push_back( static_cast<double>(counts[ec]) /
                               eff_lens[trans_id] );
    }

    weights[ec] = trans_weights;
  }


  return weights;
}

std::vector<double> trunc_gaussian_fld(int start, int stop, double mean,
    double sd) {
  size_t n = stop - start;
  std::vector<double> mean_fl(n, 0.0);

  double total_mass = 0.0;
  double total_density = 0.0;

  for (size_t i = 0; i < n; ++i) {
    double x = static_cast<double>(start + i);
    x = (x - mean) / sd;
    // XXX: this isn't normalized, but it doesn't matter since it gets
    // normalized below
    double cur_density = std::exp( - 0.5 * x * x ) / sd;

    total_mass += cur_density * i;
    total_density += cur_density;
    if (total_mass > 0) {
      mean_fl[i] = total_mass / total_density;
    }
  }

  return mean_fl;
}

std::vector<int> trunc_gaussian_counts(int start, int stop, double mean,
    double sd, int total_count) {
  size_t n = stop - start;
  std::vector<int> obs_fl(n, 0);

  double total_mass = 0.0;


  for (size_t i = 0; i < n; ++i) {
    double x = static_cast<double>(start + i);
    x = (x - mean) / sd;
    double cur_density = std::exp( - 0.5 * x * x ) / sd;
    total_mass += cur_density;
  }

  for (size_t i = 0; i < n; ++i) {
    double x = static_cast<double>(start + i);
    x = (x - mean) / sd;
    double cur_density = std::exp( - 0.5 * x * x ) / sd;
    obs_fl[i] = (int) std::round(cur_density * total_count / total_mass);
  }

  return obs_fl;
}
