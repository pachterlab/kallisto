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

    bool operator<(const block<T>& rhs) const {
        return (this->lb < rhs.lb);
    }
};

template<typename T = int>
class BlockArray {
    private:
        union {
            struct {
                block<T> b;
            } mono;

            struct {
                std::vector<block<T> > blocks;
            } poly;
        };
        // empty, mono, poly
        uint8_t flags;

    public:
        BlockArray() : flags(0) {}
        ~BlockArray() {
            clear();
        }

        BlockArray(const BlockArray& o) {
            flags = o.flags;
            if (flags == 1) {
                mono.b = o.mono.b;
            } else if (flags == 2) {
                poly.blocks = o.poly.blocks;
            }
        }

        BlockArray(BlockArray&& o) {
            flags = o.flags;
            o.flags = 0;
            if (flags == 1) {
                mono.b = std::move(o.mono.b);
            } else if (flags == 2) {
                poly.blocks = std::move(o.poly.blocks);
            }
        }

        BlockArray& operator=(const BlockArray& o) {

            clear();
            flags = o.flags;
            if (flags == 1) {
                mono.b = o.mono.b;
            } else if (flags == 2) {
                poly.blocks = o.poly.blocks;
            }
            return *this;
        }

        BlockArray& operator=(BlockArray&& o){

            if (this != &o) {

                clear();
                flags = o.flags;
                o.flags = 0;

                if (flags == 1) {
                    mono.b = std::move(o.mono.b);
                } else if (flags == 2){
                    poly.blocks = std::move(o.poly.blocks);
                }
            }

            return *this;
        }

        void insert(uint32_t lb, uint32_t ub, T val) {

            // TODO: Make sure block does not overlap with existing block
            if (flags == 0) {
                mono.b = block<T>(lb, ub, val);
                ++flags;
                return;
            } else if (flags == 1) {
                block<T> b = mono.b;
                poly.push_back(std::move(b));
                ++flags;
            }

            poly.blocks.emplace_back(lb, ub, val);
            std::sort(poly.blocks.begin(), poly.blocks.end());
        }

        const T& operator[](size_t idx) const {

            if (flags == 1) {
                return mono.b.val;
            }

            // TODO: Find upper bound without having to create temp struct
            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.blocks.begin(), poly.blocks.end(), tmp);
            // if (poly.blocks.size() > 0) {
            //     block<T> bb = poly.blocks[poly.blocks.size() - 1];
            // }

            if (ub == poly.blocks.begin()) {
                // No elements with lb <= idx
                // TODO: Handle this gracefully
                std::cout << "not found????" << std::endl;
            }
            return (--ub)->val;
        }

        std::pair<uint32_t, T> block_index(uint32_t idx) const {
            if (flags == 1) {
                return std::make_pair(0, mono.b.val);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.blocks.begin(), poly.blocks.end(), tmp);
            --ub;
            return std::make_pair(ub - poly.blocks.begin(), (ub)->val);
        }

        // Returns the ECs of the block up to and including the EC at idx
        std::vector<T> get_leading_vals(size_t idx) const {
            std::vector<T> vals;
            if (flags == 1) {
                vals.push_back(mono.b.val);
            } else {
                vals.reserve(poly.blocks.size());
                for (const auto& b : poly.blocks) {
                    if (b.lb <= idx) {
                        vals.push_back(b.val);
                    }
                }
            }

            return vals;
        }

        std::pair<uint32_t, uint32_t> get_block_at(size_t idx) const {
            if (flags == 1) {
                return std::make_pair(mono.b.lb, mono.b.ub);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.blocks.begin(), poly.blocks.end(), tmp);

            if (ub == poly.blocks.begin()) {
                // No elements with lb <= idx
                return std::make_pair(-1, -1);
            }

            --ub;

            return std::make_pair(ub->lb, ub->ub);
        }

        void reserve(uint32_t sz) {
            if (sz < 2) return;
            flags = 2;
            poly.blocks.reserve(sz);
        }

        void clear() {
            if (flags > 1) {
                poly.blocks.clear();
            }
            flags = 0;
        }

        size_t size() const {
            if (flags < 2) return flags;
            return poly.blocks.size();
        }

        size_t length() const {
            if (flags == 1) return mono.b.ub;
            return poly.blocks[poly.blocks.size()-1].ub;
        }

        void get_vals(std::vector<T>& vals) const {
            vals.clear();
            if (flags == 1) {
                vals.push_back(mono.b.val);
            } else {
                vals.reserve(poly.blocks.size());
                for (const auto& b : poly.blocks) {
                    vals.push_back(b.val);
                }
            }
        }

        typename std::vector<block<T> >::iterator begin() {
            return poly.blocks.begin();
        }

        typename std::vector<block<T> >::iterator end() {
            return poly.blocks.end();
        }

        void print() const {
            for (const auto& b : poly.blocks) {
                std::cout << "[" << b.lb << ", " << b.ub << "): " << b.val << ", ";
            }
            std::cout << std::endl;
        }

        // Extracts the slice [lb, ub) from *this and shifts it such that it
        // occupies the range [0, ub-lb).
        BlockArray<T> get_slice(uint32_t lb, uint32_t ub) const {
            BlockArray<T> slice;

            if (flags == 1) {
                slice.insert(lb, ub, mono.b.val);
            } else if (flags > 1) {

                block<T> tmp(lb);
                auto upper = std::upper_bound(poly.blocks.begin(), poly.blocks.end(), tmp);

                if (upper == poly.blocks.begin()) {
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
                } while(upper != poly.blocks.end() && upper->lb < ub);
            }

            return slice;
        }


        void serialize(std::ostream& out) const {

            out.write((char *)&flags, sizeof(flags));
            if (flags == 0) return;
            else if (flags == 1) {

                out.write((char *)&mono.b.lb, sizeof(mono.b.lb));
                out.write((char *)&mono.b.ub, sizeof(mono.b.ub));
                char* buffer = new char[mono.b.val.getSizeInBytes()];
                size_t tmp_size = mono.b.val.write(buffer);
                out.write((char *)&tmp_size, sizeof(tmp_size));
                out.write(buffer, tmp_size);
                delete[] buffer;
            } else {

                size_t tmp_size = poly.blocks.size();
                out.write((char *)&tmp_size, sizeof(tmp_size));
                for (const auto b : poly.blocks) {
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
            in.read((char *)&flags, sizeof(flags));
            if (flags == 0) return;
            else if (flags == 1) {

                in.read((char *)&mono.b.lb, sizeof(mono.b.lb));
                in.read((char *)&mono.b.ub, sizeof(mono.b.ub));
                in.read((char *)&roaring_size, sizeof(roaring_size));
                char* buffer = new char[roaring_size];
                in.read(buffer, roaring_size);
                mono.b.val.read(buffer);
                delete[] buffer;
                buffer = nullptr;
            } else {

                size_t tmp_size, lb, ub;
                T val;
                in.read((char *)&tmp_size, sizeof(tmp_size));
                poly.blocks.reserve(tmp_size);

                for (size_t i = 0; i < tmp_size; ++i) {
                    in.read((char *)&lb, sizeof(lb));
                    in.read((char *)&ub, sizeof(ub));
                    in.read((char *)&roaring_size, sizeof(roaring_size));
                    char* buffer = new char[roaring_size];
                    in.read(buffer, roaring_size);
                    val = val.read(buffer);
                    delete[] buffer;
                    buffer = nullptr;
                    poly.blocks.emplace_back(lb, ub, val);
                }
            }
        }
};

#endif // BLOCK_ARRAY_H
