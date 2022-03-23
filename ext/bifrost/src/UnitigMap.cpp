#include "UnitigMap.hpp"

UnitigMapBase::UnitigMapBase(const size_t length) : len(length), dist(0), size(0), strand(true), isEmpty(true) {}

UnitigMapBase::UnitigMapBase(const size_t start, const size_t length, const size_t unitig_sz, const bool strand) :
                            dist(start), len(length), size(unitig_sz), strand(strand), isEmpty(false) {}

bool UnitigMapBase::operator==(const UnitigMapBase& o) const {

    return (len == o.len) && (dist == o.dist) && (size == o.size) && (strand == o.strand) && (isEmpty == o.isEmpty);
}

bool UnitigMapBase::operator!=(const UnitigMapBase& o) const {

    return (len != o.len) || (dist != o.dist) || (size != o.size) || (strand != o.strand) || (isEmpty != o.isEmpty);
}
