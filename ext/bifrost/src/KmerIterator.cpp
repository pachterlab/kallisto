#include <iterator>
#include <utility>
#include "Kmer.hpp"
#include "KmerIterator.hpp"

KmerIterator& KmerIterator::operator++() {

    if (!invalid) {

        while (str[pos_e] != '\0') {

            const char c = str[pos_e] & 0xDF; // mask lowercase bit

            if (isDNA(c)) {

                if (pos_s + Kmer::k - 1 == pos_e){

                    if (pos_s == p.second + 1) p.first.selfForwardBase(c);
                    else p.first = Kmer(str + pos_s);

                    p.second = pos_s;

                    ++pos_s;
                    ++pos_e;

                    return *this;
                }
            }
            else pos_s = pos_e + 1;

            ++pos_e;
        }

        invalid = true;
    }

    return *this;
}

KmerIterator& KmerIterator::operator+=(const int len){

    if (!invalid) {

        if (len == 1) operator++();
        else if (len > 1) {

            const int next_pos_e = pos_e + len - 1;

            while ((pos_e < next_pos_e) && (str[pos_e] != '\0')) ++pos_e;

            if (str[pos_e] != '\0') {

                pos_s = pos_e - Kmer::k + 1;
                pos_e = pos_s;

                operator++();
            }
            else invalid = true;
        }
    }

    return *this;
}
