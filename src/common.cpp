#include <algorithm>
#include <unordered_map>
#include "common.h"
#include <iostream>

// Look-up maps for comma-free code
// amino acid -> comma-free nuc_seq
std::unordered_map<char, std::string>cfc_aa_map = {
    {'F', "TAG"},
    {'L', "GGG"},
    {'I', "GTG"},
    {'M', "CAG"},
    {'V', "AGC"},
    {'S', "TTC"},
    {'P', "AGT"},
    {'T', "GAC"},
    {'A', "TCG"},
    {'Y', "TAC"},
    {'H', "GCA"},
    {'Q', "TGA"},
    {'N', "GTA"}, 
    {'K', "CAT"},
    {'D', "CTA"},
    {'E', "CGA"},
    {'C', "CGT"},
    {'W', "ACT"},
    {'R', "ATG"},
    {'G', "GCT"},
    {'X', "NNN"},  // Amino acid not known
    {'B', "CTA"},  // Represents either N or D - will translate as D here (N is only off by one base)
    {'J', "GGG"},  // Represents either L or I - will translate as L here (I is only off by one base)
    {'Z', "CGA"}   // Represents either E or Q - will translate as E here (Q is only off by one base)
};

// nucleotide seq -> reverse complement
std::string revcomp(const std::string s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
      switch(c) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case 'a': return 'T';
      case 'c': return 'G';
      case 'g': return 'C';
      case 't': return 'A';
      default: return 'N';
      }
      return 'N';
    });
  return r;
}

std::string pretty_num(int num) {
  return pretty_num(static_cast<size_t>(num));
}

std::string pretty_num(int64_t num) {
  if (num < 0) {
    return "-" + pretty_num(static_cast<size_t>(num));
  } else {
    return pretty_num(static_cast<size_t>(num));
  }
}

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");

  if (s.size() <= 3) {
    return s;
  }

  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }

  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }

  ret += s.substr(start_pos, 3);

  return ret;
}
