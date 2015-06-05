#include "weights.h"
#include <cmath>

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

std::vector<double> update_eff_lens(double mean, const MinCollector& tc,  const KmerIndex &index, const std::vector<double> alpha, const std::vector<double> eff_lens) {

  double biasDataNorm = 0.0;
  double biasAlphaNorm = 0.0;
  const int num6mers = 4096;
  for (int i = 0; i < num6mers; i++) {
    biasDataNorm += tc.bias5[i];
  }

  std::vector<double> dbias5(num6mers);

  index.loadTranscriptSequences();

  for (int i = 0; i < index.num_trans; i++) {
    if (index.trans_lens_[i] < mean) {
      continue;
    }

    if (alpha[i] < 1e-8) {
      continue;
    }
    
    double contrib = 0.5*alpha[i]/eff_lens[i];
    int seqlen = index.target_seqs_[i].size();
    const char* cs = index.target_seqs_[i].c_str();

    int hex = hexamerToInt(cs,false);
    int fwlimit = (int) (seqlen - mean - 6);
    for (int j = 0; j < fwlimit; j++) {
      dbias5[hex] += contrib;
      hex = update_hexamer(hex,*(cs+j+6),false);
    }

    int bwlimit = (int) (mean - 6);
    hex = hexamerToInt(cs+bwlimit,true);
    for (int j = bwlimit; j < seqlen - 6; j++) {
      dbias5[hex] += contrib;
      if (j < seqlen - 6) {
        hex = update_hexamer(hex,*(cs+j+6),true);
      }
    }
  }

  for (int i = 0; i < num6mers; i++) {
    biasAlphaNorm += dbias5[i];
  }

  std::vector<double> biaslens(index.num_trans);

  for (int i = 0; i < index.num_trans; i++) {
    double efflen = 0.0;
    if (index.trans_lens_[i] >= mean && alpha[i] >= 1e-8) {

      int seqlen = index.target_seqs_[i].size();
      const char* cs = index.target_seqs_[i].c_str();

      // forward direction
      int hex = hexamerToInt(cs,false);
      int fwlimit = (int) seqlen - mean-6;
      for (int j = 0; j < fwlimit; j++) {
        //int hex = hexamerToInt(cs+j,false);
        efflen += 0.5*(tc.bias5[hex]/biasDataNorm) / (dbias5[hex]/biasAlphaNorm );
        hex = update_hexamer(hex,*(cs+j+6),false);
      }
      int bwlimit = (int) std::max(mean-6,0.0);
      hex = hexamerToInt(cs+bwlimit,true);
      for (int j = bwlimit; j < seqlen - 6; j++) {
        efflen += 0.5*(tc.bias5[hex]/biasDataNorm) / (dbias5[hex]/biasAlphaNorm );
        if (j < seqlen-6) {
          hex = update_hexamer(hex,*(cs+j+6),true);
        }
      }
    }

    
    if (efflen > mean) {
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

//  for (auto& kv : ecmap) {
  for (int ec = 0; ec < ecmap.size(); ec++) {
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
