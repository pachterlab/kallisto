#ifndef BIFROST_BIT_CONTAINER_HPP
#define BIFROST_BIT_CONTAINER_HPP

#include "roaring.hh"
#include "TinyBitmap.hpp"


class BitContainer {

    //Ensure that BitContainer::setPtrBmp is always allocated with an 8 bytes alignment
    struct alignas(8) Bitmap { Roaring r; };

    public:

        /** @class BitContainer_const_iterator
        * @brief See BitContainer::const_iterator
        */
        class BitContainer_const_iterator : public std::iterator<std::forward_iterator_tag, uint32_t> {

            friend class BitContainer;

            public:

                BitContainer_const_iterator();
                BitContainer_const_iterator(const BitContainer_const_iterator& o);

                ~BitContainer_const_iterator();

                BitContainer_const_iterator& operator=(const BitContainer_const_iterator& o);

                inline uint32_t operator*() const {

                    return static_cast<uint32_t>(ck_id);
                }


                BitContainer_const_iterator operator++(int);
                BitContainer_const_iterator& operator++();

                bool operator==(const BitContainer_const_iterator& o) const;
                bool operator!=(const BitContainer_const_iterator& o) const;

            private:

                const BitContainer* cs;

                size_t flag;

                size_t it_setBits;
                size_t cs_sz;

                uint64_t ck_id;

                const Roaring empty_roar;

                TinyBitmap t_bmp;

                Roaring::const_iterator it_roar;
                TinyBitmap::const_iterator it_t_bmp;

                BitContainer_const_iterator(const BitContainer* cs_, const bool beg);

                inline bool isInvalid() const {

                    return ((ck_id == 0xffffffffffffffff) || (it_setBits == cs_sz));
                }
        };

        typedef BitContainer_const_iterator const_iterator;

        BitContainer();
        BitContainer(const BitContainer& o); // Copy constructor
        BitContainer(BitContainer&& o); // Move  constructor

        ~BitContainer();

        BitContainer& operator=(const BitContainer& o);
        BitContainer& operator=(BitContainer&& o);

        bool operator==(const BitContainer& o) const;
        inline bool operator!=(const BitContainer& o) const { return !operator==(o); }

        BitContainer operator|(const BitContainer& rhs) const;
        BitContainer& operator|=(const BitContainer& rhs);

        BitContainer operator&(const BitContainer& rhs) const;
        BitContainer& operator&=(const BitContainer& rhs);

        BitContainer operator-(const BitContainer& rhs) const;
        BitContainer& operator-=(const BitContainer& rhs);

        void clear();

        void add(const size_t pair_id);
        void remove(const size_t pair_id);
        bool contains(const size_t pair_id) const;

        inline bool isEmpty() const { return (size() == 0); }

        size_t maximum() const;
        size_t minimum() const;

        size_t size() const;
        inline size_t cardinality() const { return size(); }

        bool write(ostream& stream_out) const;
        bool read(istream& stream_in);

        size_t getSizeInBytes() const;

        const_iterator begin() const;
        const_iterator end() const;

        void runOptimize();
        Roaring toRoaring() const;

        static BitContainer fastunion(const size_t sz, const BitContainer** p_id) {

            BitContainer p;

            if (sz > 0){

                p = *(p_id[0]);

                for (size_t i = 1; i < sz; ++i) p |= *(p_id[i]);
            }

            return p;
        }

    private:

        void addSortedVector(const vector<uint32_t>& v);
        void removeSortedVector(const vector<uint32_t>& v);

        inline void releaseMemory(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrBitmap) delete getPtrBitmap();
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.clear();
            }

            setBits = localBitVector;
        }

        inline void shrinkSize(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrBitmap) getPtrBitmap()->r.shrinkToFit();
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.shrinkSize();

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            }
        }

        inline bool isBitmap() const { return ((setBits & flagMask) == ptrBitmap); }
        inline bool isTinyBitmap() const { return ((setBits & flagMask) == localTinyBitmap); }

        inline Bitmap* getPtrBitmap() const {

            return reinterpret_cast<Bitmap*>(setBits & pointerMask);
        }

        inline const Bitmap* getConstPtrBitmap() const {

            return reinterpret_cast<const Bitmap*>(setBits & pointerMask);
        }

        inline uint16_t* getPtrTinyBitmap() const {

            return reinterpret_cast<uint16_t*>(setBits & pointerMask);
        }

        static const size_t maxBitVectorIDs; // 64 bits - 3 bits for the color set type = 61
        static const size_t shiftMaskBits; // 3 bits

        // asBits and asPointer represent:
        // Flag 0 - A TinyBitmap which can contain up to 65488 uint
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 integers
        // Flag 2 - A single integer
        // Flag 3 - A pointer to a CRoaring compressed bitmap which can contain up to 2^32 uint

        static const uintptr_t localTinyBitmap; // Flag 0
        static const uintptr_t localBitVector; // Flag 1
        static const uintptr_t localSingleInt; // Flag 2
        static const uintptr_t ptrBitmap; // Flag 3

        static const uintptr_t flagMask; // 0x7 (= 2^shiftMaskBits - 1)
        static const uintptr_t pointerMask; // 0xfffffffffffffff8 (= 2^64 - 1 - flagMask)

        uintptr_t setBits;
};

#endif
