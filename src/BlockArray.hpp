#ifndef BLOCK_ARRAY_H
#define BLOCK_ARRAY_H

#include <vector>
#include <algorithm>

template<typename T = int>
struct block {
    uint32_t lb, ub;
    T val;

    block() : lb(0), ub(0) {}
    block(uint32_t idx) : lb(idx) {}
    block(uint32_t lb, uint32_t ub, T val) : lb(lb), ub(ub), val(val) {}
    ~block() {}

    bool operator<(const block<T>& rhs) const {
        return (this->lb < rhs.lb);
    }
};

template<typename T = int>
class BlockArray {
    private:
        union {
            block<T> mono;
            std::vector<block<T> > poly;
        };
        // empty, mono, poly
        enum {EMPTY, MONO, POLY} flag;

    public:
        BlockArray() : flag(EMPTY), mono(0) {
        }
        ~BlockArray() {
            switch(flag) {
                case EMPTY:
                case MONO:
                    mono.~block();
                    break;
                case POLY:
                    poly.~vector();
            }
        }

        BlockArray(const BlockArray& o) {

            flag = o.flag;
            switch(flag) {
                case EMPTY:
                case MONO:
                    new (&mono) auto(o.mono);
                    break;
                case POLY:
                    new (&poly) auto(o.poly);
            }
        }

        BlockArray(BlockArray&& o) {

            flag = o.flag;
            o.flag = EMPTY;

            switch(flag) {
                case EMPTY:
                case MONO:
                    new (&mono) auto(o.mono);
                    o.mono = block<T>();
                    break;
                case POLY:
                    new (&poly) auto(o.poly);
                    o.poly.clear();
            }
        }

        BlockArray& operator=(const BlockArray& o) {

            clear();
            flag = o.flag;
            switch(flag) {
                case EMPTY:
                case MONO:
                    new (&mono) auto(o.mono);
                    break;
                case POLY:
                    new (&poly) auto(o.poly);
            }
            return *this;
        }

        BlockArray& operator=(BlockArray&& o){

            if (this != &o) {

                flag = o.flag;
                o.flag = EMPTY;

                switch(flag) {
                    case EMPTY:
                    case MONO:
                        new (&mono) auto(o.mono);
                        o.mono = block<T>();
                        break;
                    case POLY:
                        new (&poly) auto(o.poly);
                        o.poly.clear();
                }
            }

            return *this;
        }

        void insert(uint32_t lb, uint32_t ub, T val) {

            // TODO: Make sure block does not overlap with existing block
            if (flag == EMPTY) {
                mono.lb = lb;
                mono.ub = ub;
                mono.val = val;
                flag = MONO;

                return;
            } else if (flag == MONO) {
                block<T> b = mono;
                mono.~block();
                new (&poly) std::vector<block<T> >;
                poly.push_back(b);
                flag = POLY;
            }

            poly.emplace_back(lb, ub, val);
            std::sort(poly.begin(), poly.end());
        }

        const T& operator[](size_t idx) const {

            if (flag == MONO) {
                return mono.val;
            }

            // TODO: Find upper bound without having to create temp struct
            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);
            // if (poly.blocks.size() > 0) {
            //     block<T> bb = poly.blocks[poly.blocks.size() - 1];
            // }

            if (ub == poly.begin()) {
                // No elements with lb <= idx
                // TODO: Handle this gracefully
                std::cout << "not found????" << std::endl;
            }
            return (--ub)->val;
        }

        std::pair<uint32_t, T> block_index(uint32_t idx) const {
            if (flag == MONO) {
                return std::make_pair(0, mono.val);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);
            --ub;
            return std::make_pair(ub - poly.begin(), (ub)->val);
        }

        // Returns the ECs of the block up to and including the EC at idx
        std::vector<T> get_leading_vals(size_t idx) const {
            std::vector<T> vals;
            if (flag == MONO) {
                vals.push_back(mono.val);
            } else {
                vals.reserve(poly.size());
                for (const auto& b : poly) {
                    if (b.lb <= idx) {
                        vals.push_back(b.val);
                    }
                }
            }

            return vals;
        }

        std::pair<uint32_t, uint32_t> get_block_at(size_t idx) const {
            if (flag == MONO) {
                return std::make_pair(mono.lb, mono.ub);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);

            if (ub == poly.begin()) {
                // No elements with lb <= idx
                return std::make_pair(-1, -1);
            }

            --ub;

            return std::make_pair(ub->lb, ub->ub);
        }

        void reserve(uint32_t sz) {
            if (sz < 2) return;
            flag = POLY;
            poly.reserve(sz);
        }

        void clear() {
            switch(flag) {
                case EMPTY:
                case MONO:
                    mono.~block();
                    break;
                case POLY:
                    poly.~vector();
            }
            flag = EMPTY;
        }

        size_t size() const {
            if (flag == EMPTY) return 0;
            if (flag == MONO) return 1;
            return poly.size();
        }

        size_t length() const {
            if (flag == MONO) return mono.ub;
            return poly[poly.size()-1].ub;
        }

        void get_vals(std::vector<T>& vals) const {
            vals.clear();
            if (flag == MONO) {
                vals.push_back(mono.val);
            } else {
                vals.reserve(poly.size());
                for (const auto& b : poly) {
                    vals.push_back(b.val);
                }
            }
        }

        typename std::vector<block<T> >::iterator begin() {
            return poly.begin();
        }

        typename std::vector<block<T> >::iterator end() {
            return poly.end();
        }

        void print() const {
            for (const auto& b : poly) {
                std::cout << "[" << b.lb << ", " << b.ub << "): " << b.val << ", ";
            }
            std::cout << std::endl;
        }

        // Extracts the slice [lb, ub) from *this and shifts it such that it
        // occupies the range [0, ub-lb).
        BlockArray<T> get_slice(uint32_t lb, uint32_t ub) const {
            BlockArray<T> slice;

            if (flag == MONO) {
                slice.insert(lb, ub, mono.val);
            } else if (flag == POLY) {

                block<T> tmp(lb);
                auto upper = std::upper_bound(poly.begin(), poly.end(), tmp);

                if (upper == poly.begin()) {
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
                } while(upper != poly.end() && upper->lb < ub);
            }

            return slice;
        }


        void serialize(std::ostream& out) const {

            out.write((char *)&flag, sizeof(flag));
            if (flag == EMPTY) return;
            else if (flag == MONO) {

                out.write((char *)&mono.lb, sizeof(mono.lb));
                out.write((char *)&mono.ub, sizeof(mono.ub));
                char* buffer = new char[mono.val.getSizeInBytes()];
                size_t tmp_size = mono.val.write(buffer);
                out.write((char *)&tmp_size, sizeof(tmp_size));
                out.write(buffer, tmp_size);
                delete[] buffer;
            } else {

                size_t tmp_size = poly.size();
                out.write((char *)&tmp_size, sizeof(tmp_size));
                for (const auto b : poly) {
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
        }

        void deserialize(std::istream& in) {

            clear();

            size_t roaring_size;
            in.read((char *)&flag, sizeof(flag));
            if (flag == EMPTY) return;
            else if (flag == MONO) {

                in.read((char *)&mono.lb, sizeof(mono.lb));
                in.read((char *)&mono.ub, sizeof(mono.ub));
                in.read((char *)&roaring_size, sizeof(roaring_size));
                char* buffer = new char[roaring_size];
                in.read(buffer, roaring_size);
                mono.val.read(buffer);
                delete[] buffer;
                buffer = nullptr;
            } else {

                size_t tmp_size, lb, ub;
                T val;
                in.read((char *)&tmp_size, sizeof(tmp_size));
                poly.reserve(tmp_size);

                for (size_t i = 0; i < tmp_size; ++i) {
                    in.read((char *)&lb, sizeof(lb));
                    in.read((char *)&ub, sizeof(ub));
                    in.read((char *)&roaring_size, sizeof(roaring_size));
                    char* buffer = new char[roaring_size];
                    in.read(buffer, roaring_size);
                    val = val.read(buffer);
                    delete[] buffer;
                    buffer = nullptr;
                    poly.emplace_back(lb, ub, val);
                }
            }
        }
};

#endif // BLOCK_ARRAY_H
