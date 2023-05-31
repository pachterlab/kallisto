#ifndef BIFROST_DATA_ACCESSOR_HPP
#define BIFROST_DATA_ACCESSOR_HPP

#include "CompactedDBG.hpp"
#include "ColorSet.hpp"

/** @file src/DataAccessor.hpp
* Interface for the class DataAccessor. The purpose of a DataAccessor object is to provide access
* to the colors and the data associated with a unitig of a ColoredCDBG. Code snippets using this
* interface are provided in snippets/test.cpp.
*/

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

/** @class DataAccessor
* @brief Interface to access the colors and the data associated with a unitig of a ColoredCDBG.
* The class as one template parameter: the type of data associated with the unitigs of the graph.
*/
template<typename Unitig_data_t = void>
class DataAccessor : public CDBG_Data_t<DataAccessor<Unitig_data_t>, DataStorage<Unitig_data_t>> {

    typedef Unitig_data_t U;

    template<typename U> friend class DataStorage;

    public:

        /** Constructor (set up an empty DataAccessor).
        */
        DataAccessor(const uint8_t id = 0);

        /** Clear the colors and data associated with a colored unitig.
        * This function clears the k-mer color sets (of type UnitigColors) associated with the reference
        * unitig of the input parameter um (of type UnitigColorMap). It also clears the associated data
        * if there are some.
        */
        void clear(const UnitigColorMap<U>& um);

        /** Get the unitig data.
        * @param um is a constant reference to a const_UnitigColorMap object for which the reference unitig
        * used in the mapping is associated to the current DataAccessor object.
        * @return a constant pointer to the unitig data.
        */
        const U* getData(const const_UnitigColorMap<U>& um) const;

        /** Get the unitig data.
        * @param um is a constant reference to a UnitigColorMap object for which the reference unitig
        * used in the mapping is associated to the current DataAccessor object.
        * @return a pointer to the unitig data.
        */
        U* getData(const UnitigColorMap<U>& um) const;

        /** Get the colors of the reference unitig.
        * @param um is a constant reference to a const_UnitigColorMap object for which the colors of the
        * reference unitig used in the mapping must be obtained.
        * @return a constant pointer to a UnitigColors object representing the colors of the
        * reference unitig.
        */
        const UnitigColors* getUnitigColors(const const_UnitigColorMap<U>& um) const;

        /** Get the colors of the reference unitig.
        * @param um is a constant reference to a const_UnitigColorMap object for which the colors of the
        * reference unitig used in the mapping must be obtained.
        * @return a pointer to a UnitigColors object representing the colors of the reference unitig.
        */
        UnitigColors* getUnitigColors(const UnitigColorMap<U>& um) const;

        /** Create a new UnitigColors object for a unitig B corresponding to a unitig mapping to a
        * reference unitig A, i.e, B = A[um.dist..um.dist + um.len + k - 1] or
        * B = rev(A[um.dist..um.dist + um.len + k - 1]) if um.strand == false.
        * @param um is a constant reference to a const_UnitigColorMap object which is a unitig mapping from
        * which we want to obtain a new UnitigColors object.
        * @return a UnitigColors object for which the k-mers and colors match the input unitig mapping.
        */
        UnitigColors getSubUnitigColors(const const_UnitigColorMap<U>& um) const;

        /** Obtain the name of the colors present on AT LEAST one k-mer of a unitig mapping.
        * @param um is a constant reference to a const_UnitigColorMap object which is a unitig mapping.
        * @return a vector of string, each string is the name of a color.
        */
        vector<string> getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const;

        /**
        * Join data and colors of two colored unitigs which are going to be concatenated. Specifically, if A is the reference unitig
        * of the UnitigColorMap um_dest and B is the reference unitig of the UnitigColorMap um_src, then after the call to
        * this function, unitigs A and B will be removed and a unitig C = AB will be added to the graph. If um_dest.strand = false,
        * then the reverse-complement of A is going to be used in the concatenation. Reciprocally, if um_src.strand = false,
        * then the reverse-complement of B is going to be used in the concatenation.
        * @param um_dest is a UnitigColorMap object representing a unitig (the reference sequence of the mapping) to which
        * another unitig is going to be appended.
        * @param um_src is a UnitigColorMap object representing a unitig (and its data) that will be appended at the end of
        * the unitig represented by parameter um_dest.
        */
        void concat(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src);

        /**
        * Merge the data and colors of a sub-unitig B to the data and colors of a sub-unitig A.
        * The object calling this function represents the data associated with the reference unitig of um_dest.
        * The two unitigs A and B are NOT guaranteed to be from the same graph. The data of each unitig can be accessed
        * through the UnitigMap::getData.
        * @param um_dest is a UnitigColorMap object representing a sub-unitig (the mapped sequence of the mapping) A. The object
        * calling this function represents the data associated with the reference unitig of um_dest.
        * @param um_src is a UnitigColorMap object representing a sub-unitig (the mapped sequence of the mapping) for which the
        * data must be merged with the data of sub-unitig B (given by parameter um_dest).
        */
        void merge(const UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src);

        /** Extract data and colors corresponding to a sub-unitig of a unitig A. The extracted sub-unitig, called B in the following, is defined
        * as a mapping to A given by the input UnitigColorMap object um_src. Hence, B = A[um_src.dist, um_src.dist + um_src.len + k - 1]
        * or B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1]) if um_src.strand == false (B is reverse-complemented). After this
        * function returns, unitig A is deleted from the graph and B is inserted in the graph (along with their data and colors) IF the
        * input parameter last_extraction == true.
        * @param um_src is a UnitigColorMap object representing the mapping to a colored unitig A from which a new colored unitig B will be
        * extracted, i.e, B = A[um_src.dist, um_src.dist + um_src.len + k - 1] or B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1])
        * if um_src.strand == false.
        * @param last_extraction is a boolean indicating if this is the last call to this function on the reference unitig A used for the
        * mapping given by um_src. If last_extraction is true, the reference unitig A of um_src will be removed from the graph right after
        * this function returns. Also, all unitigs B extracted from the reference unitig A, along with their data and colors, will be
        * inserted in the graph.
        */
        void extract(const UnitigColorMap<U>& um_src, const bool last_extraction);

        /** Serialize the data to a GFA-formatted string. This function is used when the graph is written to disk in GFA format.
        * If the returned string is not empty, the string is appended as an optional field to the Segment line matching the unitig to which
        * this data is associated. Note that it is your responsability to add GFA-compatible tags matching your data in the string.
        * @param um_src is a const_UnitigColorMap object representing the (reference) unitig to which the data to serialize is
        * associated.
        * @return a GFA-formatted string which is the serialization of the data.
        */
        string serialize(const const_UnitigColorMap<U>& um_src) const;

    private:

        inline uint8_t get() const { return da_id; }

        inline void set(const uint8_t id) { da_id = id; }

        uint8_t da_id;
};

#endif
