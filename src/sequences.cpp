#include <sequences.h>

std::string KmerView::next()
{
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
