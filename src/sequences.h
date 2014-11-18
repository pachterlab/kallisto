#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <string>

inline char complement(const char& s);

std::string reverse_complement(const std::string& seq);

class KmerView {
    public:
        KmerView(std::string* seq, size_t k) : k_(k)
        {
            i_ = 0;
            seq_ = seq;
        }
        ~KmerView() {}
        bool next(std::string& s);
        bool done();
    private:
            std::string* seq_;
            size_t k_;
            size_t i_;
};

#endif
