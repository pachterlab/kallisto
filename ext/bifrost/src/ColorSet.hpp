#ifndef BIFROST_COLORSET_HPP
#define BIFROST_COLORSET_HPP

#include "roaring.hh"

#include "CompactedDBG.hpp"
#include "TinyBitmap.hpp"

/** @file src/ColorSet.hpp
* Interface for UnitigColors, the unitig container of k-mer color sets used in ColoredCDBG.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

/** @class UnitigColors
* @brief Represent the k-mer color sets of a unitig.
*/
class UnitigColors {

    //Ensure that UnitigColors::setPtrBmp is always allocated with an 8 bytes alignment
    struct alignas(8) Bitmap { Roaring r; };

    template<typename U> friend class ColoredCDBG;
    template<typename U> friend class DataAccessor;
    template<typename U> friend class DataStorage;

    public:

        typedef pair<UnitigColors, size_t> SharedUnitigColors;

        /** @class UnitigColors_const_iterator
        * @brief See UnitigColors::const_iterator
        */
        class UnitigColors_const_iterator : public std::iterator<std::forward_iterator_tag, pair<size_t, size_t>> {

            friend class UnitigColors;

            public:

                /**
                * Constructor of an empty iterator. The resulting iterator cannott be used as it is
                * because it is not associated with any UnitigColors.
                */
                UnitigColors_const_iterator();

                /**
                * Copy constructor. After the call to this function, the same iterator exists twice in memory.
                * @param o is a constant reference to the const_iterator to copy.
                * @return a reference to the const_iterator copy.
                */
                UnitigColors_const_iterator(const UnitigColors_const_iterator& o);

                /**
                * Destructor.
                */
                ~UnitigColors_const_iterator();

                /**
                * Copy assignment operator. After the call to this function, the same iterator exists twice in memory.
                * @param o is a constant reference to the const_iterator to copy.
                * @return a reference to the const_iterator copy.
                */
                UnitigColors_const_iterator& operator=(const UnitigColors_const_iterator& o);

                /**
                * Indirection operator.
                * @return a pair p of integers representing the position of a k-mer in the unitig (p.first)
                * and the ID of the color associated with the k-mer at the given position (p.second).
                */
                inline pair<size_t, size_t> operator*() const {

                    return make_pair(ck_id % um_sz, ck_id / um_sz);
                }

                /**
                * Get the k-mer position of the k-mer visited by the iterator. It is equal to (*it).first.
                * @return the k-mer position of the k-mer visited by the iterator. It is equal to (*it).first.
                */
                inline size_t getKmerPosition() const { return ck_id % um_sz; }

                /**
                * Get the color of the k-mer visited by the iterator. It is equal to (*it).second.
                * @return the color of the k-mer visited by the iterator. It is equal to (*it).second.
                */
                inline size_t getColorID() const { return ck_id / um_sz; }

                /**
                * Postfix increment operator: it iterates over the next k-mer of the unitig having the
                * current color or the first k-mer having the next color (if all k-mers having the current
                * color have already been visited by this iterator).
                * @return a copy of the iterator before the call to this operator.
                */
                UnitigColors_const_iterator operator++(int);

                /**
                * Prefix increment operator: it iterates over the next k-mer of the unitig having the
                * current color or the first k-mer having the next color (if all k-mers having the current
                * color have already been visited by this iterator).
                * @return a reference to the current iterator.
                */
                UnitigColors_const_iterator& operator++();

                /**
                * Color increment operator: it iterates over the first k-mer position of the next color.
                * @return a reference to the current iterator.
                */
                UnitigColors_const_iterator& nextColor();

                /**
                * Equality operator.
                * @return a boolean indicating if two iterators are the same (true) or not (false).
                */
                bool operator==(const UnitigColors_const_iterator& o) const;

                /**
                * Inequality operator.
                * @return a boolean indicating if two iterators are different (true) or not (false).
                */
                bool operator!=(const UnitigColors_const_iterator& o) const;

            private:

                const UnitigColors* cs;

                size_t flag;

                size_t it_setBits;
                size_t cs_sz;
                size_t um_sz;

                size_t start_pos;
                size_t end_pos;

                uint64_t ck_id;

                const Roaring empty_roar;

                TinyBitmap t_bmp;

                UnitigColors_const_iterator* it_uc;
                Roaring::const_iterator it_roar;
                TinyBitmap::const_iterator it_t_bmp;

                UnitigColors_const_iterator(const UnitigColors* cs_, const size_t start_pos_, const size_t end_pos_,
                                            const size_t len_unitig_km, const bool beg);

                inline uint64_t get_ID() const { return ck_id; }

                inline bool isInvalid() const {

                    return ((ck_id == 0xffffffffffffffff) || (it_setBits == cs_sz));
                }
        };

        /**
        * @typedef const_iterator
        * @brief Iterator for the colors of a unitig. The iterator iterates over the colors of a
        * unitig in ascending order. For each color, it iterates over the positions of the mapped k-mers
        * (UnitigMap parameter of UnitigColors::begin()) having this color, in ascending order.
        */
        typedef UnitigColors_const_iterator const_iterator;

        /**
        * Constructor (set up an empty container of k-mer color sets).
        */
        UnitigColors();

        /**
        * Copy constructor. After the call to this constructor, the same UnitigColors object
        * exists twice in memory.
        * @param o is the color set to copy.
        */
        UnitigColors(const UnitigColors& o); // Copy constructor

        /**
        * Move constructor. After the call to this constructor, the UnitigColors to move is empty:
        * its content has been transfered (moved) to a new UnitigColors.
        * @param o is the color set to move.
        */
        UnitigColors(UnitigColors&& o); // Move  constructor

        /**
        * Destructor.
        */
        ~UnitigColors();

        /**
        * Copy assignment operator. After the call to this operator, the same UnitigColors object
        * exists twice in memory.
        * @param o is the UnitigColors to copy.
        * @return a reference to the current UnitigColors which is a copy of o.
        */
        UnitigColors& operator=(const UnitigColors& o);

        /**
        * Move assignment operator. After the call to this operator, the UnitigColors to move is empty:
        * its content has been transfered (moved) to another UnitigColors.
        * @param o is the UnitigColors to move.
        * @return a reference to the current UnitigColors having (owning) the content of o.
        */
        UnitigColors& operator=(UnitigColors&& o);

        //bool operator==(const UnitigColors& o) const;
        // inline bool operator!=(const UnitigColors& o) const { return !operator==(o);

        /**
        * Empty a UnitigColors of its content.
        */
        void clear();

        /**
        * Check if two UnitigColors are equal.
        * @param um is a UnitigMapBase object representing the reference unitig to which the current
        * UnitigColors calling this function (this) is associated.
        * @param o is the UnitigColors to compare to the current UnitigColors calling this function (this).
        * @param um_o is a UnitigMapBase object representing the reference unitig to which the input
        * UnitigColors o (parameter) is associated.
        * @return a boolean indicating whether the two UnitigColors are equal.
        */
        bool isEqual(const UnitigMapBase& um, const UnitigColors& o, const UnitigMapBase& um_o) const;

        /**
        * Check if a UnitigColors is empty (no colors).
        * @return a boolean indicating if the UnitigColors is empty.
        */
        inline bool isEmpty() const { return (size() == 0); }

        /**
        * Add a color in the current UnitigColors to all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig for which the
        * color must be added. The color will be added only for the given mapping, i.e, unitig[um.dist..um.dist+um.len+k-1]
        * @param color_id is the ID of the color to add.
        */
        void add(const UnitigMapBase& um, const size_t color_id);

        /**
        * Remove a color in the current UnitigColors for all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig for which the
        * color must be removed. The color will be removed only for the given mapping, i.e, unitig[um.dist..um.dist+um.len+k-1]
        * @param color_id is the ID of the color to remove.
        */
        void remove(const UnitigMapBase& um, const size_t color_id);

        /**
        * Check if a color is present on all k-mers of a unitig mapping.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig. All k-mers of this mapping will be
        * checked for the presence of the color. If true is returned, all k-mers of the mapping have the color.
        * If false is returned, at least one k-mer of the mapping does not have the color.
        * @param color_id is the ID of the color to check the presence.
        * @return a boolean indicating if the color is present on all k-mers of the unitig mapping.
        */
        bool contains(const UnitigMapBase& um, const size_t color_id) const;

        /**
        * Get the number of pairs (k-mer position, color) of a reference unitig.
        * @return Number of pairs (k-mer position, color) of a reference unitig.
        */
        size_t size(const UnitigMapBase& um) const;

        /**
        * Get the number of k-mers of a reference unitig having a given color.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig.
        * @param color_id is the color index.
        * @return Number of k-mers of a reference unitig having a given color.
        */
        size_t size(const UnitigMapBase& um, const size_t color_id) const;

        /**
        * Get the largest color index of all k-mers of a reference unitig.
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig.
        * @return The largest color index of all k-mers of a reference unitig.
        */
        size_t colorMax(const UnitigMapBase& um) const;

        /**
        * Write a UnitigColors to a stream.
        * @param stream_out is an out stream to which the UnitigColors must be written. It must be
        * opened prior to the call of this function and it won't be closed by this function.
        * @return a boolean indicating if the write was successful.
        */
        bool write(ostream& stream_out, const bool copy_UnitigColors = true) const;

        /**
        * Read a UnitigColors from a stream.
        * @param stream_in is an in stream from which the UnitigColors must be read. It must be
        * opened prior to the call of this function and it won't be closed by this function.
        * @return a boolean indicating if the write was successful.
        */
        bool read(istream& stream_in);

        /**
        * Size of the UnitigColors in bytes.
        * @return Size of the UnitigColors in bytes.
        */
        size_t getSizeInBytes() const;

        /**
        * Create a constant iterator on all pairs (k-mer position, color) of the UnitigColors.
        * The iterator goes from the smallest to the largest color index. For each such color index,
        * the iterator goes from the smallest to the largest k-mer position.
        * @return a constant iterator to the smallest pair (k-mer position, color) of the UnitigColors.
        */
        const_iterator begin(const UnitigMapBase& um) const;

        /**
        * Create a constant iterator to the "past-the-last" pair (k-mer position, color) of the
        * UnitigColors.
        * @return a constant iterator to the "past-the-last" pair (k-mer position, color) of the
        * UnitigColors.
        */
        const_iterator end() const;

        /**
        * If possible, decrease the memory usage of the UnitigColors by optimizing the memory for "full
        * colors" (a color is "full" when it is present on all k-mers of the reference unitig).
        * @param um is a UnitigMapBase object representing a mapping to a reference unitig.
        * @return a boolean indicating if it was possible to optimize the memory usage of the UnitigColors.
        */
        bool optimizeFullColors(const UnitigMapBase& um);

        uint64_t hash(const size_t seed = 0) const;

    private:

        UnitigColors(SharedUnitigColors& o);
        UnitigColors(const UnitigColors& o, const SharedUnitigColors* old_ref_uc, const SharedUnitigColors* new_ref_uc);

        UnitigColors& operator=(SharedUnitigColors& o);

        void add(const size_t color_id);
        bool contains(const size_t color_km_id) const;

        inline void releaseMemory(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrUnitigColors) delete[] getPtrUnitigColors();
            else if (flag == ptrBitmap) delete getPtrBitmap();
            else if (flag == ptrSharedUnitigColors){

                SharedUnitigColors* s_uc = getPtrSharedUnitigColors();

                if (--(s_uc->second) == 0) s_uc->first.clear();
            }
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.clear();
            }

            setBits = localBitVector;
        }

        inline void shrinkSize(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrUnitigColors){

                UnitigColors* uc = getPtrUnitigColors();

                uc[0].shrinkSize();
                uc[1].shrinkSize();
            }
            else if (flag == ptrBitmap) getPtrBitmap()->r.shrinkToFit();
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.shrinkSize();

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            }
        }

        UnitigColors makeFullColors(const UnitigMapBase& um) const;
        UnitigColors getFullColors(const UnitigMapBase& um) const;
        UnitigColors getNonFullColors(const UnitigMapBase& um, const UnitigColors& full_uc) const;

        inline UnitigColors* getFullColorsPtr() {

            return (isUnitigColors() ? getPtrUnitigColors() : nullptr);
        }

        inline const UnitigColors* getFullColorsPtr() const {

            return (isUnitigColors() ? getPtrUnitigColors() : nullptr);
        }

        inline bool isBitmap() const { return ((setBits & flagMask) == ptrBitmap); }
        inline bool isTinyBitmap() const { return ((setBits & flagMask) == localTinyBitmap); }
        inline bool isUnitigColors() const { return ((setBits & flagMask) == ptrUnitigColors); }
        inline bool isSharedUnitigColors() const { return ((setBits & flagMask) == ptrSharedUnitigColors); }

        size_t size() const;

        UnitigColors reverse(const UnitigMapBase& um) const;

        const_iterator begin(const size_t start_pos, const size_t end_pos, const size_t len_km_sz) const;

        inline Bitmap* getPtrBitmap() const {

            return reinterpret_cast<Bitmap*>(setBits & pointerMask);
        }

        inline const Bitmap* getConstPtrBitmap() const {

            return reinterpret_cast<const Bitmap*>(setBits & pointerMask);
        }

        inline uint16_t* getPtrTinyBitmap() const {

            return reinterpret_cast<uint16_t*>(setBits & pointerMask);
        }

        inline UnitigColors* getPtrUnitigColors() const {

            return reinterpret_cast<UnitigColors*>(setBits & pointerMask);
        }

        inline const UnitigColors* getConstPtrUnitigColors() const {

            return reinterpret_cast<const UnitigColors*>(setBits & pointerMask);
        }

        inline SharedUnitigColors* getPtrSharedUnitigColors() const {

            return reinterpret_cast<SharedUnitigColors*>(setBits & pointerMask);
        }

        inline const SharedUnitigColors* getConstPtrSharedUnitigColors() const {

            return reinterpret_cast<const SharedUnitigColors*>(setBits & pointerMask);
        }

        static const size_t maxBitVectorIDs; // 64 bits - 3 bits for the color set type = 61
        static const size_t shiftMaskBits; // 3 bits

        // asBits and asPointer represent:
        // Flag 0 - A TinyBitmap which can contain up to 65488 uint
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 integers
        // Flag 2 - A single integer
        // Flag 3 - A pointer to a CRoaring compressed bitmap which can contain up to 2^32 uint
        // Flag 4 - A pointer to an array of 2 UnitigColors:
        //          1 - Contains "full" colors -> color is present on ALL k-mers of the unitig
        //          2 - Contains colors for k-mers if NOT full colors
        // Flag 5 - A pointer to a pair (UnitigColors, size_t) shared by multiple UnitigColors

        static const uintptr_t localTinyBitmap; // Flag 0
        static const uintptr_t localBitVector; // Flag 1
        static const uintptr_t localSingleInt; // Flag 2
        static const uintptr_t ptrBitmap; // Flag 3
        static const uintptr_t ptrUnitigColors; // Flag 4
        static const uintptr_t ptrSharedUnitigColors; // Flag 5

        static const uintptr_t flagMask; // 0x7 (= 2^shiftMaskBits - 1)
        static const uintptr_t pointerMask; // 0xfffffffffffffff8 (= 2^64 - 1 - flagMask)

        uintptr_t setBits;
};

struct UnitigColorsHash {

    size_t operator()(const UnitigColors& uc) const {

        return uc.hash();
    }
};

#endif
