#ifndef KALLISTO_PSEUDOBAM_H
#define KALLISTO_PSEUDOBAM_H


#include <vector>
#include <iostream>
#include <utility>
#include <htslib/sam.h>

#include "KmerIndex.h"

void outputPseudoBam(const KmerIndex &index, const std::vector<int> &u,
                    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
                    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
                    bool paired, bam_hdr_t *h, samFile *fp);
void revseq(char *b1, char *b2, const char *s, const char *q, int n);
void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength);
bam_hdr_t* createPseudoBamHeader(const KmerIndex& index);


struct PseudoAlignmentInfo {
  int32_t id; // id local to batch
  bool paired;  // 1
  bool r1empty; // 2
  bool r2empty; // 4  
  int32_t ec_id;
  std::vector<int32_t> u;
  PseudoAlignmentInfo() : id(-1), r1empty (true), r2empty(true), paired(true), ec_id(-1) {}
};


struct PseudoAlignmentBatch {
  int32_t batch_id;
  std::vector<PseudoAlignmentInfo> aln;
  PseudoAlignmentBatch() : batch_id(-1) {}
};



void writePseudoAlignmentBatch(std::ofstream& of, const PseudoAlignmentBatch& batch);
void readPseudoAlignmentBatch(std::ifstream& in, PseudoAlignmentBatch& batch);


#endif // KALLISTO_PSEUDOBAM_H