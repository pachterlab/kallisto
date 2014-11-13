#include <sequences.h>

inline char complement(const char& s)
{
    switch( s )
    {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }

}

std::string reverse_complement(const std::string& seq)
{
    std::string rev_seq;
    rev_seq.resize(seq.size());

    size_t i = 0;

    for (auto c = seq.crbegin(); c != seq.crend(); ++c)
    {
        rev_seq[i] = complement(*c);
        ++i;
    }


    return rev_seq;
}

std::string KmerView::next()
{
    if (done())
        throw std::out_of_range("No more kmers left.");

    std::string ret( seq_->substr(i_, k_) );
    ++i_;

    return ret;
}

bool KmerView::done()
{
    if (i_ + k_ > seq_->size())
        return true;
    return false;
}
