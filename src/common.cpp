#include "common.h"

// Added by Laura
#include <map>
// Create map cfc_map as a look-up table for the comma-free code
// nucleotide triplet -> comma-free triplet
std::map<std::string, std::string>cfc_map = {
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
  {"GGG", "TGG"}
};

// function to transform nucleotide seq to its complement
char complement(char n) {   
    switch(n)
    {   
    case 'A':
        return 'T';
        break;
    case 'T':
        return 'A';
        break;
    case 'G':
        return 'C';
        break;
    case 'C':
        return 'G';
        break;
    case 'N':
        return 'N';
        break;
    case 'a':
        return 't';
        break;
    case 't':
        return 'a';
        break;
    case 'g':
        return 'c';
        break;
    case 'c':
        return 'g';
        break;
    case 'n':
        return 'n';
        break;
    default:
        return 'N';
    }
};
// End Laura

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
