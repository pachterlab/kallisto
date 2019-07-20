#include <zlib.h>
#include <string>

#ifndef __BAM_UTIL_HPP__
#define __BAM_UTIL_HPP__

#define BAMBUFSIZE 960

typedef struct BAMHeader {
  char magic[4];
  int32_t l_text;
  std::string text;
  int32_t n_ref;
} BAMHeader;

typedef struct BAMRefInfo {
  int32_t l_name;
  std::string name;
  int32_t l_ref;
} BAMRefInfo;

typedef struct BAMAuxData {
  char tag[2];
  char val_type;
} BAMAuxData;

typedef struct BAMAlignment {
  int32_t block_size;
  int32_t refID;
  int32_t pos;
  uint8_t l_read_name;
  uint8_t mapq;
  uint16_t bin;
  uint16_t n_cigar_op;
  uint16_t flag;
  int32_t l_seq;
  int32_t next_pos;
  int32_t tlen;
  char pad[4];
} BAMAlignment;

typedef struct BAMSequence {
/*std::string seq;
  std::string qual;
  std::string name;
  std::string umi;
  std::string bc;*/
  char bc[32];
  char umi[32];
  char seq[BAMBUFSIZE];
  uint8_t nh;
  int l_bc;
  int l_umi;
  int l_seq;
} BAMSequence;

void parseBAMHeader(gzFile &f, BAMHeader &h);
void parseBAMRefInfo(gzFile &f, std::vector<BAMRefInfo> &v, int32_t n_ref);
int getBAMSequence(gzFile &f, BAMSequence &s, char *buf);

#endif

