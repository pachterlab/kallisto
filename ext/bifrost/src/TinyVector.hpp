#ifndef BIFROST_TINYVECTOR_HPP
#define BIFROST_TINYVECTOR_HPP

#include <vector>
#include <algorithm>
#include <stdint.h>

template<class T, int N = (1 + (48 / sizeof(T)))>
class tiny_vector {

    public:

        tiny_vector() : short_(true) { arr.size = 0; }

        tiny_vector(std::initializer_list<T> l) {

            short_ = true;
            arr.size = 0;

            reserve(l.size());

            for (auto &x : l) push_back(x);
        }

        ~tiny_vector() { clear(); }

        union {

            struct {
                T data[N];
                uint8_t size;
            } arr;

            struct {
                T *data;
                size_t cap;
                size_t size; // 24 bytes + ~24 malloc overhead
            } vec;
        };

        tiny_vector(const tiny_vector& o) {

            short_ = true;
            arr.size = 0;

            reserve(o.capacity());

            for (auto &x : o) push_back(x);

        }

        tiny_vector(const std::vector<T>& o) {

            short_ = true;
            arr.size = 0;

            reserve(o.capacity());

            for (auto &x : o) push_back(x);
        }

        tiny_vector(tiny_vector&& o) {

            if (o.isShort()) {

                const size_t sz = o.size();

                for (size_t i = 0; i < sz; i++) arr.data[i] = std::move(o.arr.data[i]);

                o.arr.size = 0;

                arr.size = sz;
                short_ = true;
            }
            else {

                vec.data = o.vec.data;
                vec.size = o.vec.size;
                vec.cap  = o.vec.cap;

                o.vec.data = nullptr;
                o.arr.size = 0;
                o.short_ = true;

                short_ = false;
            }
        }

        tiny_vector& operator=(const tiny_vector& o) {

            clear();
            reserve(o.capacity());

            for (auto &x : o) push_back(x);

            return *this;
        }

        tiny_vector& operator=(tiny_vector&& o){

            if (this != &o) {

                clear();

                if (o.isShort()) {

                    const size_t sz = o.size();

                    for (size_t i = 0; i < sz; i++) arr.data[i] = std::move(o.arr.data[i]);

                    o.arr.size = 0;
                    arr.size = sz;
                    short_ = true;
                }
                else {

                    vec.data = o.vec.data;
                    vec.size = o.vec.size;
                    vec.cap  = o.vec.cap;
                    o.vec.data = nullptr;
                    o.arr.size = 0;
                    o.short_ = true;
                    short_ = false;
                }
            }

            return *this;
        }

        bool short_;

        typedef T* iterator;
        typedef T const* const_iterator;

        iterator begin() { return getPointer(); }
        const_iterator begin() const { return getPointer(); }
        iterator end() { return getPointer() + size(); }
        const_iterator end() const { return (getPointer() + size()); }

        inline bool isShort() const { return short_; }

        inline T* getPointer() {

            if (isShort()) return &arr.data[0];
            return vec.data;
        }

        inline const T* getPointer() const {

            if (isShort()) return &arr.data[0];
            return vec.data;
        }

        inline size_t size() const {

            if (isShort()) return arr.size;
            return vec.size;
        }

        inline bool empty() const { return (size()==0); }

        bool operator==(const tiny_vector& o) const {

            if (size() == o.size()) {

                for (auto it=begin(), oit = o.begin(); it != end(); ++it, ++oit) {

                    if (*it != *oit) return false;
                }

                return true;
            }

            return false;
        }

        bool operator!=(const tiny_vector& o) const {

            return !operator==(o);
        }

        inline size_t capacity() const {

            if (isShort()) return N;
            return vec.cap;
        }

        T& operator[](size_t i) {

            return *(begin() + i);
        }

        const T& operator[](size_t i) const {

            return *(begin() + i);
        }

        void push_back(const T& value) {

            if (size() >= capacity()) _reserve_and_copy(std::max((size_t) 2, 3 * size() / 2));

            *(end()) = value;

            isShort() ? ++arr.size : ++vec.size;
        }

        void insert(const T& value, const size_t position) {

            if (size() >= capacity()) _reserve_and_copy(std::max((size_t) 2, 3 * size() / 2));

            T* data = getPointer();

            memmove(&data[position + 1], &data[position], (size() - position) * sizeof(T));

            data[position] = value;

            isShort() ? ++arr.size : ++vec.size;
        }

        void remove(const size_t position) {

            T* data = getPointer();

            if (position != size() - 1)
                memmove(&data[position], &data[position + 1], (size() - position - 1) * sizeof(T));

            isShort() ? --arr.size : --vec.size;
        }

        void clear() {

            if (!isShort()) {

                delete[] vec.data;
                vec.data = nullptr;
                vec.size = 0;
                vec.cap = 0;
            }
            else {

                for (auto it = begin(); it != end(); ++it) it->~T(); // manually call destructor
                arr.size = 0;
            }

            short_ = true;
        }

        inline void reserve(size_t sz) {

            if (sz > capacity()) _reserve_and_copy(sz);
        }

        void _reserve_and_copy(size_t sz) {

            if (sz <= capacity()) return;

            size_t old_size = size();

            T* newdata = new T[sz];
            // move old elements over
            size_t i = 0;
            for (auto it = begin(); it != end(); ++it, ++i) newdata[i] = std::move(*it);

            if (!isShort()) delete[] vec.data;

            short_ = false;
            vec.size = old_size;
            vec.cap = sz;
            vec.data = newdata;
        }
};

// A packed_tiny_vector is a special type of tiny_vector for size_t integer type.
// The particularity of this class is that it is not self-contained: an external flag must be used
// (and stored externally) to know how the structure is allocated internally.
// * Advantage is that if only one integer is stored, the structure takes only 9 bytes instead of
// 32 for the tiny_vector.
// * Numerous disadvantages:
// - the structure must be deallocated manually (no destructor)
// - Size and capacity can only be accessed with one cache-miss
// - No assignment operator
// - etc.
// TO USE WITH CAUTION!
class packed_tiny_vector {

    public:

        static const uint8_t FLAG_EMPTY = 0;
        static const uint8_t FLAG_STATIC_ALLOC = 1;
        static const uint8_t FLAG_DYNAMIC_ALLOC = 2;

        packed_tiny_vector(const size_t val = 0) : bits(val) {}

        inline void destruct(const uint8_t flag){

            if (flag == FLAG_DYNAMIC_ALLOC) delete[] ptr;
        }

        packed_tiny_vector(const packed_tiny_vector& o, const uint8_t flag) {

            if (flag == FLAG_DYNAMIC_ALLOC){

                const size_t capacity_sz = o.capacity(flag) + 2;

                ptr = new size_t[capacity_sz];

                memcpy(ptr, o.ptr, capacity_sz * sizeof(size_t));
            }
            else bits = o.bits;
        }

        /*packed_tiny_vector(packed_tiny_vector&& o) {

            bits = o.bits;
            o.clear();
        }*/

        packed_tiny_vector& copy(uint8_t& flag, const packed_tiny_vector& o, const uint8_t flag_o) {

            destruct(flag);

            if (flag_o == FLAG_DYNAMIC_ALLOC){

                const size_t capacity_sz = o.capacity(flag_o) + 2;

                ptr = new size_t[capacity_sz];

                memcpy(ptr, o.ptr, capacity_sz * sizeof(size_t));
            }
            else bits = o.bits;

            flag = flag_o;

            return *this;
        }

        packed_tiny_vector& move(uint8_t& flag, packed_tiny_vector&& o, uint8_t&& flag_o) {

            if (this != &o) {

                destruct(flag);

                bits = o.bits;
                o.clear();

                flag = flag_o;
                flag_o = FLAG_EMPTY;
            }

            return *this;
        }

        inline size_t size(const uint8_t flag) const {

            if (flag == FLAG_DYNAMIC_ALLOC) return ptr[0];
            if (flag == FLAG_STATIC_ALLOC) return 1;
            return 0;
        }

        inline size_t capacity(const uint8_t flag) const {

            if (flag == FLAG_DYNAMIC_ALLOC) return ptr[1];
            return 1;
        }

        inline bool empty(const uint8_t flag) const { return (flag == FLAG_EMPTY); }

        inline void clear() { bits = 0; }

        size_t& operator()(const size_t pos, const uint8_t flag) {

            if (flag == FLAG_DYNAMIC_ALLOC) return ptr[pos+2];
            return bits;
        }

        const size_t& operator()(const size_t pos, const uint8_t flag) const {

            if (flag == FLAG_DYNAMIC_ALLOC) return ptr[pos+2];
            return bits;
        }

        uint8_t push_back(const size_t val, const uint8_t flag) {

            const size_t sz = size(flag);

            uint8_t new_flag = flag;

            if (flag == FLAG_EMPTY) new_flag = FLAG_STATIC_ALLOC;
            else if (flag == FLAG_STATIC_ALLOC) new_flag = FLAG_DYNAMIC_ALLOC;

            if (sz >= capacity(flag)) _reserve_and_copy(std::max((size_t) 2, 3 * sz / 2), flag);

            operator()(sz, new_flag) = val;

            if (new_flag == FLAG_DYNAMIC_ALLOC) ++ptr[0];

            return new_flag;
        }

        uint8_t insert(const size_t val, const size_t pos, const uint8_t flag) {

            const size_t sz = size(flag);

            uint8_t new_flag = flag;

            if (flag == FLAG_EMPTY) new_flag = FLAG_STATIC_ALLOC;
            else if (flag == FLAG_STATIC_ALLOC) new_flag = FLAG_DYNAMIC_ALLOC;

            if (sz >= capacity(flag)) _reserve_and_copy(std::max((size_t) 2, 3 * sz / 2), flag);

            if (new_flag == FLAG_DYNAMIC_ALLOC){

                memmove(&ptr[pos + 3], &ptr[pos + 2], (sz - pos) * sizeof(size_t));

                ptr[pos + 2] = val;

                ++ptr[0];
            }
            else if (pos == 0) bits = val;

            return new_flag;
        }

        uint8_t remove(const size_t pos, const uint8_t flag) {

            if (flag == FLAG_DYNAMIC_ALLOC) {

                const size_t new_sz = size(flag) - 1;

                if (pos < new_sz) memmove(&ptr[pos + 2], &ptr[pos + 3], (new_sz - pos) * sizeof(size_t));

                if (new_sz == 1){

                    const size_t val = ptr[2];

                    delete[] ptr;

                    bits = val;

                    return FLAG_STATIC_ALLOC;
                }

                --ptr[0];

                return FLAG_DYNAMIC_ALLOC;
            }

            return FLAG_EMPTY;
        }

    private:

        union {

            size_t bits;
            size_t* ptr;
        };

        void _reserve_and_copy(const size_t sz, const uint8_t flag) {

            if (sz <= capacity(flag)) return;

            size_t old_sz = size(flag);
            size_t* newdata = new size_t[sz+2];

            for (size_t i = 0; i != old_sz; ++i) newdata[i+2] = operator()(i, flag);

            if (flag == FLAG_DYNAMIC_ALLOC) delete[] ptr;

            ptr = newdata;
            ptr[0] = old_sz;
            ptr[1] = sz;
        }
};

#endif
