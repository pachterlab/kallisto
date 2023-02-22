#include <algorithm>
#include <unordered_map>
#include "common.h"

// Look-up maps for comma-free code
// amino acid -> comma-free triplet
std::unordered_map<char, std::string>cfc_aa_map = {
    {'F', "ACC"},
    {'L', "ACA"},
    {'I', "ATA"},
    {'M', "ATC"},
    {'V', "ATT"},
    {'S', "CTA"},
    {'P', "CTC"},
    {'T', "CTT"},
    {'A', "AGA"},
    {'Y', "AGC"},
    {'H', "AGT"},
    {'Q', "AGG"},
    {'N', "CGA"}, 
    {'K', "CGC"},
    {'D', "CGT"},
    {'E', "CGG"},
    {'C', "TGA"},
    {'W', "TGC"},
    {'R', "TGT"},
    {'G', "TGG"},
    {'X', "NNN"},  // Amino acid not known
    {'B', "CGT"},  // Represents either N or D - will translate as D here (N is only off by one base)
    {'J', "ACA"},  // Represents either L or I - will translate as L here (I is only off by one base)
    {'Z', "CGG"}   // Represents either E or Q - will translate as E here (Q is only off by one base)
};

// nucleotide triplet -> comma-free triplet
std::unordered_map<std::string, std::string>cfc_map = {
  {"TTT", "ACC"},
  {"TTC", "ACC"},
  {"TTA", "ACA"},
  {"TTG", "ACA"},
  {"CTT", "ACA"},
  {"CTC", "ACA"},
  {"CTA", "ACA"},
  {"CTG", "ACA"},
  {"ATT", "ATA"},
  {"ATC", "ATA"},
  {"ATA", "ATA"},
  {"ATG", "ATC"},
  {"GTT", "ATT"},
  {"GTC", "ATT"},
  {"GTA", "ATT"},
  {"GTG", "ATT"},
  {"TCT", "CTA"},
  {"TCC", "CTA"},
  {"TCA", "CTA"},
  {"TCG", "CTA"},
  {"AGT", "CTA"},
  {"AGC", "CTA"},
  {"CCT", "CTC"},
  {"CCC", "CTC"},
  {"CCA", "CTC"},
  {"CCG", "CTC"},
  {"ACT", "CTT"},
  {"ACC", "CTT"},
  {"ACA", "CTT"},
  {"ACG", "CTT"},
  {"GCT", "AGA"},
  {"GCC", "AGA"},
  {"GCA", "AGA"},
  {"GCG", "AGA"},
  {"TAT", "AGC"},
  {"TAC", "AGC"},
  {"CAT", "AGT"},
  {"CAC", "AGT"},
  {"CAA", "AGG"},
  {"CAG", "AGG"},
  {"AAT", "CGA"},
  {"AAC", "CGA"},
  {"AAA", "CGC"},
  {"AAG", "CGC"},
  {"GAT", "CGT"},
  {"GAC", "CGT"},
  {"GAA", "CGG"},
  {"GAG", "CGG"},
  {"TGT", "TGA"},
  {"TGC", "TGA"},
  {"TGG", "TGC"},
  {"CGT", "TGT"},
  {"CGC", "TGT"},
  {"CGA", "TGT"},
  {"CGG", "TGT"},
  {"AGA", "TGT"},
  {"AGG", "TGT"},
  {"GGT", "TGG"},
  {"GGC", "TGG"},
  {"GGA", "TGG"},
  {"GGG", "TGG"},
  {"TGA", "NNN"}, // stop codon
  {"TAA", "NNN"}, // stop codon
  {"TAG", "NNN"}  // stop codon
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
