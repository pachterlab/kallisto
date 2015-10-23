#include <vector>
#include <iostream>
#include <utility>

#include "KmerIndex.h"

void outputPseudoBam(const KmerIndex &index, const std::vector<int> &u,
                    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
                    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
                    bool paired);
void revseq(char *b1, char *b2, const char *s, const char *q, int n);
void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength);
