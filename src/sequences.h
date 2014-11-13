#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <string>

class KmerView {
    public:
        KmerView(std::string* seq, size_t k) : k_(k)
        {
            i_ = 0;
            seq_ = seq;
        }
        ~KmerView() {}
        std::string next();
        bool done();
    private:
            std::string* seq_;
            size_t k_;
            size_t i_;
};

#endif
