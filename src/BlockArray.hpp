#ifndef BLOCK_ARRAY_H
#define BLOCK_ARRAY_H

#include <vector>
#include <algorithm>

template<typename T = int>
struct block {
    size_t lb, ub;
    T val;

    block(size_t idx) : lb(idx) {}
    block(size_t lb, size_t ub, T val) : lb(lb), ub(ub), val(val) {}

    bool operator<(const block<T>& rhs) const {
        return (this->lb < rhs.lb);
    }
};

template<typename T = int>
class BlockArray {
    private:
        std::vector<block<T> > blocks;

    public:
        BlockArray() {}

        void insert(size_t lb, size_t ub, T val) {

            // TODO: Make sure block does not overlap with existing block
            blocks.emplace_back(lb, ub, val);
            std::sort(blocks.begin(), blocks.end());
        }

        const T& operator[](size_t idx) const {

            // TODO: Find upper bound without having to create temp struct
            block<T> tmp(idx);
            auto ub = std::upper_bound(blocks.begin(), blocks.end(), tmp);
            if (blocks.size() > 0) {
                block<T> bb = blocks[blocks.size() - 1];
            }

            if (ub == blocks.begin()) {
                // No elements with lb <= idx
                // TODO: Handle this gracefully
                std::cout << "not found????" << std::endl;
            }
            return (--ub)->val;
        }

        std::pair<size_t, T> block_index(size_t idx) const {
            block<T> tmp(idx);
            auto ub = std::upper_bound(blocks.begin(), blocks.end(), tmp);
            --ub;
            return make_pair(ub - blocks.begin(), (ub)->val);
        }

        // Returns the ECs of the block up to and including the EC at idx
        std::vector<T> get_leading_vals(size_t idx) const {
            std::vector<T> vals;
            vals.reserve(blocks.size());
            for (const auto& b : blocks) {
                if (b.lb <= idx) {
                    vals.push_back(b.val);
                }
            }
            return vals;
        }

        std::pair<size_t, size_t> get_block_at(size_t idx) const {
            block<T> tmp(idx);
            auto ub = std::upper_bound(blocks.begin(), blocks.end(), tmp);

            if (ub == blocks.begin()) {
                // No elements with lb <= idx
                return std::make_pair(-1, -1);
            }

            --ub;

            return std::make_pair(ub->lb, ub->ub);
        }

        void reserve(size_t sz) {
            blocks.reserve(sz);
        }

        void clear() {
            blocks.clear();
        }

        size_t size() const {
            return blocks.size();
        }

        size_t length() const {
            return blocks[blocks.size()-1].ub;
        }

        void get_vals(std::vector<T>& vals) const {
            vals.clear();
            vals.reserve(blocks.size());
            for (const auto& b : blocks) {
                vals.push_back(b.val);
            }
        }

        typename std::vector<block<T> >::iterator begin() {
            return blocks.begin();
        }

        typename std::vector<block<T> >::iterator end() {
            return blocks.end();
        }

        void print() const {
            for (const auto& b : blocks) {
                std::cout << "[" << b.lb << ", " << b.ub << "): " << b.val << ", ";
            }
            std::cout << std::endl;
        }

        // Extracts the slice [lb, ub) from *this and shifts it such that it
        // occupies the range [0, ub-lb).
        BlockArray<T> get_slice(size_t lb, size_t ub) const {
            BlockArray<T> slice;

            block<T> tmp(lb);
            auto upper = std::upper_bound(blocks.begin(), blocks.end(), tmp);

            if (upper == blocks.begin()) {
                // No elements with lower bound <= lb
                // TODO: Handle this
            }

            --upper;

            do {
                slice.insert(
                        std::max(lb, upper->lb)-lb,
                        std::min(ub, upper->ub)-lb,
                        upper->val
                );
                ++upper;
            } while(upper != blocks.end() && upper->lb < ub);

            return slice;
        }


        void serialize(std::ostream& out) const {

            size_t tmp_size = blocks.size();
            out.write((char *)&tmp_size, sizeof(tmp_size));
            for (const auto b : blocks) {
                out.write((char *)&b.lb, sizeof(b.lb));
                out.write((char *)&b.ub, sizeof(b.ub));
                char* buffer = new char[b.val.getSizeInBytes()];
                tmp_size = b.val.write(buffer);
                out.write((char *)&tmp_size, sizeof(tmp_size));
                out.write(buffer, tmp_size);
                delete[] buffer;
                buffer = nullptr;
            }
        }

        void deserialize(std::istream& in) {

            blocks.clear();
            size_t tmp_size, roaring_size, lb, ub;
            T val;
            in.read((char *)&tmp_size, sizeof(tmp_size));
            blocks.reserve(tmp_size);

            for (size_t i = 0; i < tmp_size; ++i) {
                in.read((char *)&lb, sizeof(lb));
                in.read((char *)&ub, sizeof(ub));
                in.read((char *)&roaring_size, sizeof(roaring_size));
                char* buffer = new char[roaring_size];
                in.read(buffer, roaring_size);
                val = val.read(buffer);
                delete[] buffer;
                buffer = nullptr;
                blocks.emplace_back(lb, ub, val);
            }
        }
};

#endif // BLOCK_ARRAY_H
