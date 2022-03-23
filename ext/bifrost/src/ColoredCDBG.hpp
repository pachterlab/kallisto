#ifndef BIFROST_COLOREDCDBG_HPP
#define BIFROST_COLOREDCDBG_HPP

#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "CompactedDBG.hpp"
#include "DataManager.hpp"

#include "getRSS.h"

/** @file src/ColoredCDBG.hpp
* Interface for the Colored and Compacted de Bruijn graph API.
* Code snippets using this interface are provided in test/snippets.hpp.
*/

/** @struct CCDBG_Build_opt
* @brief This structure inherits from CDBG_Build_opt and introduces only a few new members which are
* color-related. See the documentation of CDBG_Build_opt for the inherited members.
* An example of using such a structure is shown in src/Bifrost.cpp.
* @var CCDBG_Build_opt::filename_colors_in
* String containing the name of a Bifrost color file to read in ColoredCDBG<U>::read(). Default is
* an empty string (no file specified).
* @var CCDBG_Build_opt::outputColors
* Boolean indicating if the graph should be colored or not. This member is not used by any function of
* ColoredCDBG<U> or CompactedDBG<U, T>. It is used by the Bifrost CLI. Default is true.
*/
struct CCDBG_Build_opt : CDBG_Build_opt {

    string filename_colors_in;

    bool outputColors;

    CCDBG_Build_opt() : outputColors(true) {}
};

template<typename U = void> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U = void> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

/** @class CCDBG_Data_t
* @brief If data are to be associated with the unitigs of the colored and compacted de Bruijn graph, those data
* must be wrapped into a class that inherits from the abstract class CCDBG_Data_t.
* To associate data of type "MyUnitigData" to unitigs, class MyUnitigData must be declared as follows:
* \code{.cpp}
* class MyUnitigData : public CCDBG_Data_t<MyUnitigData> { ... };
* ...
* ColoredCDBG<MyUnitigData> ccdbg;
* \endcode
* An object of type MyUnitigData represents an instanciation of user data associated to one unitig of the graph.
* CCDBG_Data_t has one template parameter: the type of unitig data ("MyUnitigData").
* Because CCDBG_Data_t is an abstract class, all the methods of the base class (CCDBG_Data_t) must be implemented in
* your wrapper (the derived class, aka MyUnitigData in this example). IMPORTANT: If you do not implement those methods,
* default ones that have no effects will be applied. Do not forget to implement copy and move constructors/destructors
* as well as copy and move assignment operators.
* An example of using such a structure is shown in snippets/test.cpp.
*/
template<typename Unitig_data_t> //Curiously Recurring Template Pattern (CRTP)
class CCDBG_Data_t {

    typedef Unitig_data_t U;

    public:

        /**
        * Clear the data associated with a unitig.
        * @param um_dest is a UnitigColorMap object representing a unitig (the reference sequence of um_dest) for which
        * the data must be cleared. The object calling this function represents the data associated with the reference unitig
        * of um_dest.
        */
        void clear(const UnitigColorMap<U>& um_dest){}

        /**
        * Join data of two unitigs which are going to be concatenated. Specifically, if A is the reference unitig
        * of the UnitigColorMap um_dest and B is the reference unitig of the UnitigColorMap um_src, then after the call to
        * this function, unitigs A and B will be removed and a unitig C = AB will be added to the graph. The object calling
        * this function represents the data to associate with the new unitig C = AB. If um_dest.strand = false,
        * then the reverse-complement of A is going to be used in the concatenation. Reciprocally, if um_src.strand = false,
        * then the reverse-complement of B is going to be used in the concatenation. The data of each unitig can be accessed
        * through the method UnitigColorMap::getData(). The two unitigs A and B are guaranteed to be from the same graph.
        * @param um_dest is a UnitigColorMap object representing a unitig (the reference sequence of the mapping) to which
        * another unitig is going to be appended. The object calling this function represents the data associated with the
        * reference unitig of um_dest.
        * @param um_src is a UnitigColorMap object representing a unitig (and its data) that will be appended at the end of
        * the unitig represented by parameter um_dest.
        */
        void concat(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src){}

        /**
        * Merge the data of a sub-unitig B to the data of a sub-unitig A.
        * The object calling this function represents the data associated with the reference unitig of um_dest.
        * The two unitigs A and B are NOT guaranteed to be from the same graph. The data of each unitig can be accessed
        * through the UnitigMap::getData.
        * @param um_dest is a UnitigColorMap object representing a sub-unitig (the mapped sequence of the mapping) A. The object
        * calling this function represents the data associated with the reference unitig of um_dest.
        * @param um_src is a UnitigColorMap object representing a sub-unitig (the mapped sequence of the mapping) for which the
        * data must be merged with the data of sub-unitig B (given by parameter um_dest).
        */
        void merge(const UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src){}

        /** Extract data corresponding to a sub-unitig of a unitig A. The extracted sub-unitig, called B in the following, is defined
        * as a mapping to A given by the input UnitigColorMap object um_src. Hence, B = A[um_src.dist, um_src.dist + um_src.len + k - 1]
        * or B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1]) if um_src.strand == false (B is reverse-complemented). After this
        * function returns, unitig A is deleted from the graph and B is inserted in the graph (along with their data and colors) IF the
        * input parameter last_extraction == true. The object calling this function represents the data associated with sub-unitig B.
        * @param uc_dest is a constant pointer to a UnitigColors object representing the k-mer color sets of colored unitig B.
        * If uc_dest == nullptr, no UnitigColors object will be associated with unitig B.
        * @param um_src is a UnitigColorMap object representing the mapping to a colored unitig A from which a new colored unitig B will be
        * extracted, i.e, B = A[um_src.dist, um_src.dist + um_src.len + k - 1] or B = rev(A[um_src.dist, um_src.dist + um_src.len + k - 1])
        * if um_src.strand == false.
        * @param last_extraction is a boolean indicating if this is the last call to this function on the reference unitig A used for the
        * mapping given by um_src. If last_extraction is true, the reference unitig A of um_src will be removed from the graph right after
        * this function returns. Also, all unitigs B extracted from the reference unitig A, along with their data and colors, will be
        * inserted in the graph.
        */
        void extract(const UnitigColors* uc_dest, const UnitigColorMap<U>& um_src, const bool last_extraction){}

        /** Serialize the data to a GFA-formatted string. This function is used when the graph is written to disk in GFA format.
        * If the returned string is not empty, the string is appended as an optional field to the Segment line matching the unitig to which
        * this data is associated. Note that it is your responsability to add GFA-compatible tags matching your data in the string.
        * @param um_src is a const_UnitigColorMap object representing the (reference) unitig to which the data to serialize is
        * associated.
        * @return a string which is the serialization of the data.
        */
        string serialize(const const_UnitigColorMap<U>& um_src) const {

            return string();
        }
};

/** @class ColoredCDBG
* @brief Represent a Colored and Compacted de Bruijn graph. The class inherits from CompactedDBG
* which means that all public functions available with CompactedDBG are also available with ColoredCDBG.
* \code{.cpp}
* ColoredCDBG<> ccdbg_1; // No unitig data
* ColoredCDBG<void> ccdbg_2; // Equivalent to previous notation
* ColoredCDBG<MyUnitigData> ccdbg_3; // An object of type MyUnitigData will be associated with each unitig (along the colors)
* \endcode
* If data are to be associated with the unitigs, these data must be wrapped into a class that inherits from the abstract class
* CCDBG_Data_t, such as in:
* \code{.cpp}
* class MyUnitigData : public CCDBG_Data_t<MyUnitigData> { ... };
* CompactedDBG<MyUnitigData> cdbg;
* \endcode
* Because CCDBG_Data_t is an abstract class, all the methods of the base class (CCDBG_Data_t) must be implemented in
* your wrapper (the derived class, aka MyUnitigData in this example). IMPORTANT: If you do not implement those methods, default
* ones that have no effects will be applied.
*/
template<typename Unitig_data_t = void>
class ColoredCDBG : public CompactedDBG<DataAccessor<Unitig_data_t>, DataStorage<Unitig_data_t>> {

    static_assert(is_void<Unitig_data_t>::value || is_base_of<CCDBG_Data_t<Unitig_data_t>, Unitig_data_t>::value,
                  "Type Unitig_data_t of data associated with vertices of class ColoredCDBG<Unitig_data_t> must "
                  " be void (no data) or a class extending class CCDBG_Data_t");

    typedef Unitig_data_t U;

    template<typename U> friend class DataAccessor;

    public:

        /** Constructor (set up an empty colored cdBG).
        * @param kmer_length is the length k of k-mers used in the graph (each unitig is of length at least k).
        * @param minimizer_length is the length g of minimizers (g < k) used in the graph.
        */
        ColoredCDBG(int kmer_length = DEFAULT_K, int minimizer_length = -1);

        /** Copy constructor (copy a colored cdBG).
        * This function is expensive in terms of time and memory as the content of a colored and compacted
        * de Bruijn graph is copied. After the call to this function, the same graph exists twice in memory.
        * @param o is a constant reference to the colored and compacted de Bruijn graph to copy.
        */
        ColoredCDBG(const ColoredCDBG& o);

        /** Move constructor (move a colored cdBG).
        * The content of o is moved ("transfered") to a new colored and compacted de Bruijn graph.
        * The colored and compacted de Bruijn graph referenced by o will be empty after the call to this constructor.
        * @param o is a reference on a reference to the colored and compacted de Bruijn graph to move.
        */
        ColoredCDBG(ColoredCDBG&& o);

        /** Copy assignment operator (copy a colored cdBG).
        * This function is expensive in terms of time and memory as the content of a colored and compacted
        * de Bruijn graph is copied.  After the call to this function, the same graph exists twice in memory.
        * @param o is a constant reference to the colored and compacted de Bruijn graph to copy.
        * @return a reference to the colored and compacted de Bruijn which is the copy.
        */
        ColoredCDBG& operator=(const ColoredCDBG& o);

        /** Move assignment operator (move a colored cdBG).
        * The content of o is moved ("transfered") to a new colored and compacted de Bruijn graph.
        * The colored and compacted de Bruijn graph referenced by o will be empty after the call to this operator.
        * @param o is a reference on a reference to the colored and compacted de Bruijn graph to move.
        * @return a reference to the colored and compacted de Bruijn which has (and owns) the content of o.
        */
        ColoredCDBG& operator=(ColoredCDBG&& o);

        /** Equality operator.
        * @return a boolean indicating if two compacted de Bruijn graphs have the same colored unitigs (does not
        * compare the data associated with the unitigs).
        */
        bool operator==(const ColoredCDBG& o) const;

        /** Inequality operator.
        * @return a boolean indicating if two compacted de Bruijn graphs have different colored unitigs (does not
        * compare the data associated with the unitigs).
        */
        inline bool operator!=(const ColoredCDBG& o) const;

        /** Addition assignment operator (merge a colored cdBG).
        * After merging, all unitigs and colors of o have been added to and compacted with the current colored
        * and compacted de Bruijn graph (this). If the unitigs of o had data of type "MyUnitigData" associated,
        * they have been added to the current colored and compacted de Bruijn graph using the functions of the
        * class MyUnitigData which are in base class CCDBG_Data_t<MyUnitigData>. This function is similar to
        * ColoredCDBG::merge except that it uses only one thread while ColoredCDBG::merge can work with multiple
        * threads (number of threads provided as a parameter). Note that if multiple colored and compacted
        * de Bruijn graphs have to be merged, it is more efficient to call ColoredCDBG::merge with a vector of
        * ColoredCDBG as input.
        * @param o is a constant reference to the colored and compacted de Bruijn graph to merge.
        * @return a reference to the current colored and compacted de Bruijn after merging.
        */
        ColoredCDBG& operator+=(const ColoredCDBG& o);

        /** Clear the graph: remove unitigs, user data and colors + reset its parameters.
        */
        void clear();

        /** Build the Colored and compacted de Bruijn graph (only the unitigs).
        * A call to ColoredCDBG::mapColors is required afterwards to map colors to unitigs.
        * @param opt is a structure from which the members are parameters of this function. See CCDBG_Build_opt.
        * @return boolean indicating if the graph has been built successfully.
        */
        bool buildGraph(const CCDBG_Build_opt& opt);

        /** Map the colors to the unitigs. This is done by reading the input files and querying the graph.
        * If a color filename is provided in opt.filename_colors_in, colors are loaded from that file instead.
        * @param opt is a structure from which the members are parameters of this function. See CCDBG_Build_opt.
        * @return boolean indicating if the colors have been mapped successfully.
        */
        bool buildColors(const CCDBG_Build_opt& opt);

        /** Write a colored and compacted de Bruijn graph to disk.
        * @param prefix_output_filename is a string which is the prefix of the filename for the two files that are
        * going to be written to disk. Assuming the prefix is "XXX", two files "XXX.gfa" and "XXX.bfg_colors" will
        * be written to disk.
        * @param nb_threads is the number of threads that can be used to write the graph to disk.
        * @param verbose is a boolean indicating if information message are printed during writing (true) or not (false).
        * @return a boolean indicating if the graph was successfully written.
        */
        bool write(const string& prefix_output_filename, const size_t nb_threads = 1, const bool verbose = false) const;

        /** Read a colored and compacted de Bruijn graph from disk. The graph (in GFA format) must have been produced
        * by Bifrost.
        * @param prefix_input_filename is a string which is the prefix of the filename for the two files that are
        * going to be read from disk. Assuming the prefix is "XXX", two files "XXX.gfa" and "XXX.bfg_colors" will
        * be read from disk.
        * @param nb_threads is the number of threads that can be used to read the graph and its colors from disk.
        * @param verbose is a boolean indicating if information messages are printed during reading (true) or not (false).
        * @return a boolean indicating if the graph was successfully read.
        */
        bool read(const string& input_graph_filename, const string& input_colors_filename, const size_t nb_threads = 1, const bool verbose = false);

        /** Merge a colored and compacted de Bruijn graph.
        * After merging, all unitigs and colors of the input graph have been added to and compacted with the current
        * colored and compacted de Bruijn graph (this). If the unitigs of the input graph had data of type "MyUnitigData"
        * associated, they have been added to the current colored and compacted de Bruijn graph using the functions of
        * the class MyUnitigData which are also present in its base class CCDBG_Data_t<MyUnitigData>.
        * Note that if multiple colored and compacted de Bruijn graphs have to be merged, it is more efficient to call
        * ColoredCDBG::merge with a vector of ColoredCDBG as input.
        * @param o is a constant reference to the colored and compacted de Bruijn graph to merge.
        * @param nb_threads is an integer indicating how many threads can be used during the merging.
        * @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
        * @return a boolean indicating if the graph has been successfully merged.
        */
        bool merge(const ColoredCDBG& o, const size_t nb_threads = 1, const bool verbose = false);

        /** Merge and clear a colored and compacted de Bruijn graph.
        * After merging, all unitigs and colors of the input graph have been added to and compacted with the current colored
        * and compacted de Bruijn graph (this). The input graph is cleared before the function returns. If the unitigs of the
        * input graph had data of type "MyUnitigData" associated, they have been added to the current colored and compacted
        * de Bruijn graph using the functions of the class MyUnitigData which are also present in its base class
        * CCDBG_Data_t<MyUnitigData>.
        * Note that if multiple colored and compacted de Bruijn graphs have to be merged, it is more efficient to call
        * ColoredCDBG::merge with a vector of ColoredCDBG as input.
        * @param o is a reference on a reference to the colored and compacted de Bruijn graph to merge. It can be obtained using
        * std::move(). After merging, the graph pointed by o is cleared.
        * @param nb_threads is an integer indicating how many threads can be used during the merging.
        * @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
        * @return a boolean indicating if the graph has been successfully merged.
        */
        bool merge(ColoredCDBG&& o, const size_t nb_threads = 1, const bool verbose = false);

        /** Merge multiple colored and compacted de Bruijn graphs.
        * After merging, all unitigs and colors of the input colored and compacted de Bruijn graphs have been added to and
        * compacted with the current colored and compacted de Bruijn graph (this). If the unitigs had data of type "MyUnitigData"
        * associated, they have been added to the current colored and compacted de Bruijn graph using the functions of the
        * class MyUnitigData which are also present in its base class CCDBG_Data_t<MyUnitigData>.
        * @param v is a constant reference to a vector of colored and compacted de Bruijn graphs to merge.
        * @param nb_threads is an integer indicating how many threads can be used during the merging.
        * @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
        * @return a boolean indicating if the graphs have been successfully merged.
        */
        bool merge(const vector<ColoredCDBG>& v, const size_t nb_threads = 1, const bool verbose = false);

        /** Merge and clear multiple colored and compacted de Bruijn graphs.
        * After merging, all unitigs and colors of the input colored and compacted de Bruijn graphs have been added to and
        * compacted with the current colored and compacted de Bruijn graph (this). The input graphs are cleared before the
        * function returns. If the input unitigs had data of type "MyUnitigData" associated, they have been added to the
        * current colored and compacted de Bruijn graph using the functions of the class MyUnitigData which are also present
        * in its base class CCDBG_Data_t<MyUnitigData>.
        * @param v is a reference on a reference to a vector of colored and compacted de Bruijn graphs to merge.  It can be
        * obtained using std::move(). After merging, the graphs in v are cleared.
        * @param nb_threads is an integer indicating how many threads can be used during the merging.
        * @param verbose is a boolean indicating if information messages must be printed during the execution of the function.
        * @return a boolean indicating if the graphs have been successfully merged.
        */
        bool merge(vector<ColoredCDBG>&& v, const size_t nb_threads = 1, const bool verbose = false);

        /** Get the name of a color. As colors match the input files, the color names match the input filenames.
        * @return a string which is either a color name or an empty string if the color ID is invalid or if the
        * colors have not yet been mapped to the unitigs.
        */
        string getColorName (const size_t color_id) const;

        /** Get the names of all colors. As colors match the input files, the color names match the input filenames.
        * @return a vector of strings for which each string is either a color name or an empty string if the color ID is
        * invalid or if the colors have not yet been mapped to the unitigs.
        */
        vector<string> getColorNames() const;

        /** Get the number of colors in the graph.
        * @return the number of colors in the graph.
        */
        inline size_t getNbColors() const { return this->getData()->getNbColors(); }

        bool search(const vector<string>& query_filenames, const string& out_filename_prefix,
                    const double ratio_kmers, const bool inexact_search, const size_t nb_threads,
                    const bool verbose = false) const;

    private:

        void checkColors(const vector<string>& filename_seq_in) const;

        void initUnitigColors(const CCDBG_Build_opt& opt, const size_t max_nb_hash = 31);
        void buildUnitigColors(const size_t nb_threads);
        //void buildUnitigColors2(const size_t nb_threads);

        void resizeDataUC(const size_t sz, const size_t nb_threads = 1, const size_t max_nb_hash = 31);

        bool invalid;
};

#include "ColoredCDBG.tcc"

#endif
