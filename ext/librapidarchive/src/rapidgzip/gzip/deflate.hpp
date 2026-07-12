/**
 * @file
 *
 * - Note that this implementation avoids C++ exceptions because invalid data is assumed to happen rather often,
 *   which is the case when searching for deflate blocks without knowing the exact offsets! Exceptions are too
 *   slow for that!
 * - In the same manner as exceptions, it turns out that using std::array with a maximum possible size instead of
 *   dynamically sized std::vectors improves speed for checking and decoding a lot by avoiding heap allocations.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <core/Error.hpp>
#include <rapidgzip/DecodedDataView.hpp>
#include <rapidgzip/MarkerReplacement.hpp>
//#include <huffman/HuffmanCodingLinearSearch.hpp>
//#include <huffman/HuffmanCodingSymbolsPerLength.hpp>
#include <rapidgzip/huffman/HuffmanCodingReversedBitsCachedCompressed.hpp>
#include <rapidgzip/huffman/HuffmanCodingReversedBitsCached.hpp>
//#include <rapidgzip/huffman/HuffmanCodingReversedCodesPerLength.hpp>

//#define WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER
//#define WITH_MULTI_CACHED_HUFFMAN_DECODER

#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    //#include <rapidgzip/huffman/HuffmanCodingDistanceISAL.hpp>
    #include <rapidgzip/huffman/HuffmanCodingISAL.hpp>
#elif defined( WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER )
    #include <rapidgzip/huffman/HuffmanCodingShortBitsCachedDeflate.hpp>
#elif defined( WITH_MULTI_CACHED_HUFFMAN_DECODER )
    #include <rapidgzip/huffman/HuffmanCodingShortBitsMultiCached.hpp>
#else
    //#include <rapidgzip/huffman/HuffmanCodingDoubleLiteralCached.hpp>
    #include <huffman/HuffmanCodingShortBitsCached.hpp>
#endif

#include "definitions.hpp"
#include "RFCTables.hpp"


namespace rapidgzip::deflate
{
/**
 * @verbatim
 * function benchmarkRapidgzipParallel()
 * {
 *     m rapidgzip &>/dev/null && for (( i=0; i<10; ++i)); do
 *         src/tools/rapidgzip -v -d -o /dev/null "$1" 2>&1 | sed -nr 's|.*Decompressed in total.* -> ([0-9.]+) .*|\1|p'
 *     done
 * }
 * function benchmarkRapidgzipParallelFiles()
 * {
 *     for file in test-files/silesia/20xsilesia.tar.gz test-files/fastq/10xSRR22403185_2.fastq.gz 4GiB-base64.gz; do
 *         echo "$file"
 *         uncertainValue $( benchmarkRapidgzipParallel "$file" )
 *     done
 * }
 *
 * Decompressed in total 4239155200 B from 20xsilesia.tar.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 4810 | 5024 +- 10 | 5127
 *     HuffmanCodingDoubleLiteralCached                    : 3072 | 3123 +-  4 | 3178
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 3425 | 3505 +-  4 | 3564
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 3849 | 3953 +-  6 | 4025
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 3752 | 3927 +-  8 | 4017
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 3736 | 3880 +-  6 | 3953
 *
 * Decompressed in total 3618153020 B from 10xSRR22403185_2.fastq.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 2701 | 2871 +- 10 | 3056
 *     HuffmanCodingDoubleLiteralCached                    : 2431 | 2600 +- 10 | 2719
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 2719 | 2815 +-  8 | 3000
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 2742 | 2868 +-  7 | 2945
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 2809 | 2938 +-  8 | 3046
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 2734 | 2803 +-  5 | 2888
 *
 * Decompressed in total 4294967296 B from 4GiB-base64.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 6794 | 6973 +- 9 | 7081
 *     HuffmanCodingDoubleLiteralCached                    : 3537 | 3591 +- 3 | 3650
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 3977 | 4038 +- 4 | 4096
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 3876 | 3964 +- 6 | 4065
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 3926 | 4035 +- 6 | 4096
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 3924 | 4024 +- 5 | 4079
 * @endverbatim
 *
 * BEWARE: These timings are HIGHLY dependent on something that I cannot fully reproduce.
 *         It might be RAM usage, maybe owing to my Frankensystem not being able to use full dual-channel speed
 *         on the whole addressable range.
 *            2x16GiB DIMM DDR4 Synchronous Unbuffered (Unregistered) 3600 MHz (0.3 ns)
 *            2x32GiB DIMM DDR4 Synchronous Unbuffered (Unregistered) 3600 MHz (0.3 ns)
 *         It might even be the single Youtube video running in the background, which is GPU-accellerated and
 *         does not cause any CPU utilization (0-1% in total) but still might result in context switches and/or
 *         cache interference.
 *         Interestingly, the least susceptible is base64, which is fairly constant and the most susceptible
 *         is FASTQ, which yields 2.5 GB/s for LUT_BITS_COUNT=11 for one test and 3.3 GB/s after freeing 10 GB
 *         of RAM, so ~20% variation! Silesia only changes by ~10%. These tests have been repeated 10 times
 *         during which the results are fairly stable. They only vary over longer time spans.
 *
 * -> Even though HuffmanCodingShortBitsCached is fairly simple and does not even cache longer codes and
 *    instead falls back to >bit-wise< code reading, it still outperforms the previous contender:
 *    HuffmanCodingDoubleLiteralCached. All of the test cases are faster with HuffmanCodingShortBitsCached!
 *    The highest improvements are achieved for silesia.tar.gz.
 *    We are still far away from HuffmanCodingISAL for base64.gz.
 *    For FASTQ, we are actually even with HuffmanCodingISAL!
 *    This shows how much the Huffman table creation bottle-necked the decoding.
 *    @todo Future improvements on this should also cache some of the length and distance codes
 *          following non-literal symbols and/or double-cache symbols.
 *
 * Redo non-parallelized to reduce contributions of memory bandwidth and CPU utilization etc.
 *
 * @verbatim
 * function benchmarkRapidgzipSequential()
 * {
 *     m rapidgzip &>/dev/null && for (( i=0; i<10; ++i)); do
 *         src/tools/rapidgzip -P 1 -v -d -o /dev/null "$1" 2>&1 |
 *             sed -nr 's|.*Decompressed in total.* -> ([0-9.]+) .*|\1|p'
 *     done
 * }
 * function benchmarkRapidgzipSequentialFiles()
 * {
 *     for file in test-files/silesia/silesia.tar.gz test-files/fastq/SRR22403185_2.fastq.gz base64-512MiB.gz; do
 *         echo "$file"
 *         uncertainValue $( benchmarkRapidgzipSequential "$file" )
 *     done
 * }
 *
 * Decompressed in total  B from silesia.tar.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 703.8 | 720.5 +- 1.8 | 770.6
 *     HuffmanCodingDoubleLiteralCached                    : 247.34 | 252.48 +- 0.19 | 254.12
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 269.2 | 273.2 +- 0.3 | 280.9
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 322.3 | 330.4 +- 0.4 | 335.9
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 320.1 | 327.6 +- 0.5 | 338.9
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 323.5 | 327.7 +- 0.3 | 332.5
 *     HuffmanCodingShortBitsCachedDeflate with 11 Bits    : 307.4 | 312.8 +- 0.4 | 317.7
 *     HuffmanCodingShortBitsMultiCached with 11 Bits      : 329.3 | 337.3 +- 0.8 | 356.5
 *
 * Decompressed in total  B from 10xSRR22403185_2.fastq.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 857.8 | 879.1 +- 1.2 | 896.5
 *     HuffmanCodingDoubleLiteralCached                    : 334.3 | 342.3 +- 0.4 | 351.0
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 350.67 | 356.18 +- 0.27 | 361.27
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 358.3 | 366.5 +- 0.4 | 371.2
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 356.4 | 366.8 +- 0.4 | 371.4
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 360.9 | 365.8 +- 0.3 | 371.2
 *     HuffmanCodingShortBitsCachedDeflate with 11 Bits    : 335.6 | 349.7 +- 0.6 | 357.8
 *     HuffmanCodingShortBitsMultiCached with 11 Bits      : 363.9 | 376.3 +- 0.9 | 393.2
 *
 * Decompressed in total  B from 4GiB-base64.gz in MB/s:
 *     HuffmanCodingISAL with LIBRAPIDARCHIVE_WITH_ISAL=ON : 527.2 | 538.8 +- 0.7 | 545.6
 *     HuffmanCodingDoubleLiteralCached                    : 252.9 | 254.95 +- 0.19 | 258.83
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=8  : 219.4 | 244.4 +- 1.9 | 272.6
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=10 : 210.4 | 234.6 +- 1.7 | 264.9
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=11 : 213.9 | 238.3 +- 2.0 | 262.6
 *     HuffmanCodingShortBitsCached with LUT_BITS_COUNT=12 : 209.2 | 221.2 +- 1.1 | 240.0
 *     HuffmanCodingShortBitsCachedDeflate with 11 Bits    : 201.1 | 243.3 +- 2.0 | 260.5
 *     HuffmanCodingShortBitsMultiCached with 11 Bits      : 214.9 | 229.6 +- 1.0 | 242.1
 * @endverbatim
 *
 * It really is insane how much these benchmarks differ from the multi-threaded ones.
 * While base64 is the fastest multi-threaded test case, it is the slowest using a single-thread.
 * Similarly, HuffmanCodingDoubleLiteralCached is faster with a single-thread but slower with multiple
 * threads, probably because the LUT becomes too large for the caches when two hardware threads use
 * the same core.
 */
#ifdef LIBRAPIDARCHIVE_WITH_ISAL
using LiteralOrLengthHuffmanCoding = HuffmanCodingISAL;
#elif defined( WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER )
using LiteralOrLengthHuffmanCoding = HuffmanCodingShortBitsCachedDeflate</* LUT_BITS_COUNT */ 11>;
#elif defined( WITH_MULTI_CACHED_HUFFMAN_DECODER )
using LiteralOrLengthHuffmanCoding = HuffmanCodingShortBitsMultiCached</* LUT_BITS_COUNT */ 11>;
#else
//using LiteralOrLengthHuffmanCoding =
//    HuffmanCodingDoubleLiteralCached<uint16_t, MAX_CODE_LENGTH, uint16_t, MAX_LITERAL_HUFFMAN_CODE_COUNT>;
using LiteralOrLengthHuffmanCoding = HuffmanCodingShortBitsCached<
    uint16_t, MAX_CODE_LENGTH, uint16_t, MAX_LITERAL_HUFFMAN_CODE_COUNT,
    /* LUT_BITS_COUNT */ 11, /* REVERSE_BITS */ true, /* CHECK_OPTIMALITY */ true>;
#endif

/**
 * Because the fixed Huffman coding is used by different threads it HAS TO BE immutable. It is constant anyway
 * but it also MUST NOT have mutable members. This means that HuffmanCodingDoubleLiteralCached does NOT work
 * because it internally saves the second symbol.
 * @todo Make it such that the implementations can handle the case that the construction might result in
 *       larger symbol values than are allowed to appear in the first place! I.e., cut-off construction there.
 *       Note that changing this from 286 to 512, lead to an increase of the runtime! We need to reduce it again! */
using FixedHuffmanCoding =
    HuffmanCodingReversedBitsCached<uint16_t, MAX_CODE_LENGTH, uint16_t, MAX_LITERAL_OR_LENGTH_SYMBOLS + 2>;



/**
 * [findDeflateBlocksRapidgzipLUT with 13 bits, Walk Tree LUT] ( 42.1 <= 42.6 +- 0.3 <= 43 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 14 bits, Walk Tree LUT] ( 42.14 <= 42.56 +- 0.22 <= 42.9 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 15 bits, Walk Tree LUT] ( 41.9 <= 42.4 +- 0.3 <= 42.9 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 16 bits, Walk Tree LUT] ( 40.7 <= 42 +- 0.6 <= 42.8 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 17 bits, Walk Tree LUT] ( 41.7 <= 42.5 +- 0.5 <= 43.2 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 18 bits, Walk Tree LUT] ( 40.57 <= 41.01 +- 0.28 <= 41.41 ) MB/s
 *
 * Cumulative time spent during tests with deflate::block::readDynamicHuffmanCoding:
 *     readDynamicHuffmanCoding : 10.368 s
 *     Read precode             : 0.468309 s
 *     Create precode HC        : 2.72764 s
 *     Apply precode HC         : 6.68265 s
 *     Create distance HC       : 0.138837 s
 *     Create literal HC        : 0.0533746 s
 * -> @todo The precode check is much more lax for this!
 *          The block finder performance is only comparable to others because checkPrecode is called beforehand!
 */
// using PrecodeHuffmanCoding = HuffmanCodingLinearSearch<uint8_t, uint8_t>;

/**
 * [findDeflateBlocksRapidgzipLUT with 13 bits, Walk Tree LUT] ( 50.1 <= 50.6 +- 0.4 <= 51.4 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 14 bits, Walk Tree LUT] ( 49.4 <= 50.6 +- 0.7 <= 51.4 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 15 bits, Walk Tree LUT] ( 50.04 <= 50.42 +- 0.3 <= 50.96 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 16 bits, Walk Tree LUT] ( 48.2 <= 50.1 +- 0.8 <= 51 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 17 bits, Walk Tree LUT] ( 49.3 <= 50.1 +- 0.4 <= 50.6 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 18 bits, Walk Tree LUT] ( 45.9 <= 47.8 +- 1 <= 48.8 ) MB/s
 *
 * Cumulative time spent during tests with deflate::block::readDynamicHuffmanCoding:
 *     readDynamicHuffmanCoding : 10.233 s
 *     Read precode             : 0.466273 s
 *     Create precode HC        : 2.61544 s
 *     Apply precode HC         : 6.6635 s
 *     Create distance HC       : 0.138498 s
 *     Create literal HC        : 0.0532511 s
 * -> @todo The precode check is much more lax for this!
 *          The block finder performance is only comparable to others because checkPrecode is called beforehand!
 */
// using PrecodeHuffmanCoding = HuffmanCodingSymbolsPerLength<uint8_t, MAX_PRECODE_LENGTH, uint8_t, MAX_PRECODE_COUNT>;

/**
 * [findDeflateBlocksRapidgzipLUT with 13 bits, Walk Tree LUT] ( 48.4 <= 49 +- 0.4 <= 49.8 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 14 bits, Walk Tree LUT] ( 48.3 <= 49.5 +- 0.7 <= 50.4 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 15 bits, Walk Tree LUT] ( 47.4 <= 49 +- 0.7 <= 49.6 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 16 bits, Walk Tree LUT] ( 48.9 <= 49.8 +- 0.3 <= 50.2 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 17 bits, Walk Tree LUT] ( 48.4 <= 49.6 +- 0.5 <= 50.2 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 18 bits, Walk Tree LUT] ( 46.3 <= 47.2 +- 0.5 <= 47.8 ) MB/s
 *
 * Cumulative time spent during tests with deflate::block::readDynamicHuffmanCoding:
 *     readDynamicHuffmanCoding : 1.8663 s
 *     Read precode             : 0.417445 s
 *     Create precode HC        : 1.06137 s
 *     Apply precode HC         : 0.0399441 s
 *     Create distance HC       : 0.00760224 s
 *     Create literal HC        : 0.0447448 s
 */
// using PrecodeHuffmanCoding = HuffmanCodingReversedCodesPerLength<uint8_t, MAX_PRECODE_LENGTH,
//                                                                  uint8_t, MAX_PRECODE_COUNT>;

/**
 * [findDeflateBlocksRapidgzipLUT with 13 bits, Walk Tree LUT] ( 51.2 <= 52.5 +- 0.8 <= 53.6 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 14 bits, Walk Tree LUT] ( 52.2 <= 53 +- 0.3 <= 53.3 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 15 bits, Walk Tree LUT] ( 51 <= 52 +- 0.4 <= 52.4 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 16 bits, Walk Tree LUT] ( 52.1 <= 52.6 +- 0.3 <= 53.1 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 17 bits, Walk Tree LUT] ( 51.1 <= 52.2 +- 0.5 <= 52.7 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 18 bits, Walk Tree LUT] ( 48.9 <= 50.4 +- 0.7 <= 50.9 ) MB/s
 *
 * Cumulative time spent during tests with deflate::block::readDynamicHuffmanCoding:
 *     readDynamicHuffmanCoding : 2.01957 s
 *     Read precode             : 0.409709 s
 *     Create precode HC        : 1.23679 s
 *     Apply precode HC         : 0.0188343 s
 *     Create distance HC       : 0.00916463 s
 *     Create literal HC        : 0.0444881 s
 * -> This sums up to 1.71898603 s < 2.01957 s. I can only surmise that the missing time is spend in the calls
 *    to now() and the increments for failed counters and other statistics.
 */
// using PrecodeHuffmanCoding = HuffmanCodingReversedBitsCached<uint8_t, MAX_PRECODE_LENGTH, uint8_t, MAX_PRECODE_COUNT>;

/**
 * [findDeflateBlocksRapidgzipLUT with 13 bits, Walk Tree LUT] ( 52.2 <= 52.9 +- 0.4 <= 53.7 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 14 bits, Walk Tree LUT] ( 52.6 <= 53.6 +- 0.4 <= 53.9 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 15 bits, Walk Tree LUT] ( 52.14 <= 52.42 +- 0.13 <= 52.57 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 16 bits, Walk Tree LUT] ( 52.3 <= 52.7 +- 0.3 <= 53.1 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 17 bits, Walk Tree LUT] ( 51.5 <= 53.3 +- 0.7 <= 53.8 ) MB/s
 * [findDeflateBlocksRapidgzipLUT with 18 bits, Walk Tree LUT] ( 50.2 <= 50.9 +- 0.5 <= 51.6 ) MB/s
 *
 * Cumulative time spent during tests with deflate::block::readDynamicHuffmanCoding:
 *     readDynamicHuffmanCoding : 1.84971 s
 *     Read precode             : 0.417705 s
 *     Create precode HC        : 1.06757 s
 *     Apply precode HC         : 0.0182615 s
 *     Create distance HC       : 0.00743017 s
 *     Create literal HC        : 0.0440123 s
 */
using PrecodeHuffmanCoding = HuffmanCodingReversedBitsCachedCompressed<uint8_t, MAX_PRECODE_LENGTH,
                                                                       uint8_t, MAX_PRECODE_COUNT>;

/**
 * HuffmanCodingReversedBitsCached is definitely faster for silesia.tar.gz which has more back-references than
 * base64.gz for which the difference in changing this Huffman coding is negligible. Note that we can't use
 * double caching for this because that would mean merging the cache with ne next literal/length Huffman code!
 *
 * HuffmanCodingDistanceISAL:
 *
 *     m rapidgzip && src/tools/rapidgzip -d -o /dev/null test-files/fastq/10xSRR22403185_2.fastq.gz
 *     Decompressed in total 3618153020 B in:
 *         1.60722 s -> 2251.18 MB/s
 *         1.63562 s -> 2212.1 MB/s
 *         1.63213 s -> 2216.83 MB/s
 *
 *     m rapidgzip && src/tools/rapidgzip -d -o /dev/null test-files/silesia/20xsilesia.tar.gz
 *     Decompressed in total 4239155200 B in:
 *         1.10059 s -> 3851.72 MB/s
 *         1.13037 s -> 3750.25 MB/s
 *         1.16631 s -> 3634.66 MB/s
 *         1.12481 s -> 3768.77 MB/s
 *
 * HuffmanCodingReversedBitsCached:
 *
 *     m rapidgzip && src/tools/rapidgzip -d -o /dev/null test-files/fastq/10xSRR22403185_2.fastq.gz
 *     Decompressed in total 3618153020 B in:
 *         1.61128 s -> 2245.52 MB/s
 *         1.61067 s -> 2246.36 MB/s
 *         1.65374 s -> 2187.86 MB/s
 *         1.60478 s -> 2254.61 MB/s
 *
 *     m rapidgzip && src/tools/rapidgzip -d -o /dev/null test-files/silesia/20xsilesia.tar.gz
 *     Decompressed in total 4239155200 B in:
 *         1.11193 s -> 3812.43 MB/s
 *         1.0941 s -> 3874.56 MB/s
 *         1.0993 s -> 3856.23 MB/s
 *
 * -> ISA-l is actually slightly (~1-2%) slower than my own simple distance Huffman decoder.
 *    Probably because the table is small enough that short/long caching hinders performance more than it helps.
 **/
using DistanceHuffmanCoding = HuffmanCodingReversedBitsCached<uint16_t, MAX_CODE_LENGTH,
                                                              uint8_t, MAX_DISTANCE_SYMBOL_COUNT>;
//using DistanceHuffmanCoding = HuffmanCodingDistanceISAL;

/* Include 256 safety buffer so that we can avoid branches while filling. */
using LiteralAndDistanceCLBuffer = std::array<uint8_t, MAX_LITERAL_OR_LENGTH_SYMBOLS + MAX_DISTANCE_SYMBOL_COUNT + 256>;


[[nodiscard]] constexpr FixedHuffmanCoding
createFixedHC()
{
    std::array<uint8_t, MAX_LITERAL_OR_LENGTH_SYMBOLS + 2> encodedFixedHuffmanTree{};
    for ( size_t i = 0; i < encodedFixedHuffmanTree.size(); ++i ) {
        // NOLINTBEGIN(bugprone-branch-clone)
        if ( i < 144 ) {
            encodedFixedHuffmanTree[i] = 8;
        } else if ( i < 256 ) {
            encodedFixedHuffmanTree[i] = 9;
        } else if ( i < 280 ) {
            encodedFixedHuffmanTree[i] = 7;
        } else {
            encodedFixedHuffmanTree[i] = 8;
        }
        // NOLINTEND(bugprone-branch-clone)
    }

    FixedHuffmanCoding result;
    const auto error = result.initializeFromLengths( { encodedFixedHuffmanTree.data(),
                                                       encodedFixedHuffmanTree.size() } );
    if ( error != Error::NONE ) {
        throw std::logic_error( "Fixed Huffman Tree could not be created!" );
    }

    return result;
}


namespace
{
/**
 * @note Initially this was a static member of Block but when compiling with the manylinux_2_28 container,
 *       which contains "g++ (GCC) 11.2.1 20220127 (Red Hat 11.2.1-9)", I'd get redefinition errors for this
 *       even though it compiles fine with "g++ (Ubuntu 11.2.0-19ubuntu1) 11.2.0".
 *       I suspect that it does not like static methods inside templated classes. Multiple instantiations
 *       of it might be mangled into the same name?
 *       > rapidgzip.cpp:638:1: error: redefinition of ‘const char
 *       > _ZTSZN7rapidgzip7deflate5BlockILb0ELb0EE33readDistanceAndLiteralCodeLengthsERSt5arrayIhLm572EER9BitReader
 *       > ILb0EmERKNS_41HuffmanCodingReversedBitsCachedCompressedIhLh7EhLm19EEEmRKSt8functionIFhhEEEd_UlhE_ []’
 */
[[nodiscard]] forceinline Error
readDistanceAndLiteralCodeLengths( LiteralAndDistanceCLBuffer&              literalCL,
                                   gzip::BitReader&                         bitReader,
                                   const PrecodeHuffmanCoding&              precodeCoding,
                                   const size_t                             literalCLSize,
                                   const std::function<uint8_t( uint8_t )>& translateSymbol
                                   = [] ( uint8_t symbol ) { return symbol; } )
{
    size_t i = 0;
    for ( ; i < literalCLSize; ) {
        const auto decoded = precodeCoding.decode( bitReader );
        if ( !decoded ) {
            return Error::INVALID_HUFFMAN_CODE;
        }
        const auto code = translateSymbol( *decoded );

        /* Note that this interpretation of the alphabet results in the maximum code length being 15! */
        if ( code <= 15 ) {
            literalCL[i] = code;
            ++i;
        } else if ( code == 16 ) {
            if ( i == 0 ) {
                return Error::INVALID_CL_BACKREFERENCE;
            }
            const auto lastValue = literalCL[i - 1];

            /* Unroll 3U + 0b11U = 6 times to avoid branches. Do it manually to be portable. */
            literalCL[i + 0] = lastValue;
            literalCL[i + 1] = lastValue;
            literalCL[i + 2] = lastValue;
            literalCL[i + 3] = lastValue;
            literalCL[i + 4] = lastValue;
            literalCL[i + 5] = lastValue;

            i += bitReader.read<2>() + 3;
        } else if ( code == 17 ) {
            /* Unroll 3U + 0b111U = 10 times to avoid branches. Do it manually to be portable. */
            literalCL[i + 0] = 0;
            literalCL[i + 1] = 0;
            literalCL[i + 2] = 0;
            literalCL[i + 3] = 0;
            literalCL[i + 4] = 0;
            literalCL[i + 5] = 0;
            literalCL[i + 6] = 0;
            literalCL[i + 7] = 0;
            literalCL[i + 8] = 0;
            literalCL[i + 9] = 0;

            i += bitReader.read<3>() + 3;
        } else if ( code == 18 ) {
            /* Decode fixed number of zeros. The vector is initialized to zeros, so we can simply skip these. */
            #if defined( __GNUC__ )
                #pragma GCC unroll 16
            #endif
            for ( size_t j = 0; j < 11U + ( 1U << 7U ) - 1U; ++j ) {
                literalCL[i + j] = 0;
            }
            i += bitReader.read<7>() + 11;
        }
    }

    return i == literalCLSize ? Error::NONE : Error::EXCEEDED_LITERAL_RANGE;
}
}


/**
 * It should be fine to have these data members even when not needed.
 * It's not like they are expensive to initialize and deflate::block shouldn't be created in quick successions
 * anyway, it can and should be reused!
 */
class BlockStatistics
{
public:
    uint64_t failedPrecodeInit{ 0 };
    uint64_t failedDistanceInit{ 0 };
    uint64_t failedLiteralInit{ 0 };
    uint64_t failedPrecodeApply{ 0 };
    uint64_t missingEOBSymbol{ 0 };

    std::array<uint64_t, /* codeLengthCount - 4 is 4 bits = 16 possible values */ 16> precodeCLHistogram{};

    struct {
        uint32_t precode{ 0 };
        uint32_t distance{ 0 };
        uint32_t literal{ 0 };  // Minimum value is 257!
    } codeCounts;

    /* These are cumulative counters but they can be manually reset before calls to readHeader. */
    struct {
        uint64_t literal{ 0 };
        uint64_t backreference{ 0 };
        uint64_t copies{ 0 };
    } symbolTypes;

    /* These are cumulative counters but they can be manually reset before calls to readHeader. */
    struct {
        double readDynamicHeader{ 0 };
        double readPrecode{ 0 };
        double createPrecodeHC{ 0 };
        double applyPrecodeHC{ 0 };
        double createDistanceHC{ 0 };
        double createLiteralHC{ 0 };
        double readData{ 0 };
    } durations;

    /** These are time points used to calculate the durations and are necessary to constexpr hide calls to now(). */
    struct {
        using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

        TimePoint readDynamicStart;
        TimePoint readPrecode;
        TimePoint createdPrecodeHC;
        TimePoint appliedPrecodeHC;
        TimePoint createdDistanceHC;

        TimePoint readDataStart;
    } times;
};


/**
 * @todo Silesia is ~70% slower when writing back and calculating CRC32.
 * When only writing the result and not calculating CRC32, then it is ~60% slower.
 * Both, LZ77 back-references and CRC32 calculation can still be improved upon by a lot, I think.
 * Silesia contains a lot of 258 length back-references with distance 1, which could be replaced with memset
 * with the last byte.
 */
template<bool ENABLE_STATISTICS = false>
class Block :  // NOLINT(cppcoreguidelines-pro-type-member-init)
    public BlockStatistics
{
public:
    using CompressionType = deflate::CompressionType;

    struct Backreference{
        uint16_t distance{ 0 };
        uint16_t length{ 0 };
    };

public:
    [[nodiscard]] bool
    eob() const noexcept
    {
        return m_atEndOfBlock;
    }

    [[nodiscard]] constexpr bool
    eos() const noexcept
    {
        return m_atEndOfBlock && m_isLastBlock;
    }

    [[nodiscard]] constexpr bool
    eof() const noexcept
    {
        return m_atEndOfFile;
    }

    [[nodiscard]] constexpr bool
    isLastBlock() const noexcept
    {
        return m_isLastBlock;
    }

    [[nodiscard]] constexpr CompressionType
    compressionType() const noexcept
    {
        return m_compressionType;
    }

    [[nodiscard]] constexpr uint8_t
    padding() const noexcept
    {
        return m_padding;
    }

    /**
     * @param nMaxToDecode Maximum bytes to decode. It might decode less even when there is enough data.
     *                     It will only decode as much as fits into the internal buffer.
     *                     It might decode more when it is an uncompressed block.
     *                     Check for @ref eob to test for the end of the block instead of testing the read byte count.
     */
    [[nodiscard]] std::pair<DecodedDataView, Error>
    read( gzip::BitReader& bitReader,
          size_t           nMaxToDecode = std::numeric_limits<size_t>::max() );

    /**
     * @tparam treatLastBlockAsError This parameter is intended when using readHeader for finding valid headers.
     *         Ignoring last headers, filters candidates by 25% and filtering them sooner avoids reading the
     *         Huffman codings, which saves almost 50% of time!
     */
    template<bool treatLastBlockAsError = false>
    [[nodiscard]] Error
    readHeader( gzip::BitReader& bitReader );

    /**
     * Reads the dynamic Huffman code. This is called by @ref readHeader after reading the first three header bits
     * and determining that it is a dynamic Huffman encoded block.
     * @note Normally, you want to call @ref readHeader instead. This is only for very specific edge use cases!
     */
    [[nodiscard]] Error
    readDynamicHuffmanCoding( gzip::BitReader& bitReader );

    [[nodiscard]] constexpr size_t
    uncompressedSize() const noexcept
    {
        return m_compressionType == CompressionType::UNCOMPRESSED ? m_uncompressedSize : 0;
    }

    /**
     * Primes the deflate decoder with a window to be used for the LZ77 backreferences.
     * There are two use cases for this function:
     *  - To set a window before decoding in order to resume decoding and for seeking in the gzip stream.
     *  - To replace marker bytes with real data in post.
     * The only real use case for the latter would be huge deflate blocks. In practice, all gzip implementations
     * I encountered produced deflate blocks not larger than 64 KiB. In that case, it would be simpler to create
     * a new deflate::Block object on the next block and then set the initial window before decoding with the
     * data from the last read calls whose markers will have to be replaced using @ref replaceMarkerBytes.
     * This method does not much more but has to account for wrap-around, too.
     */
    void
    setInitialWindow( VectorView<uint8_t> const& initialWindow = {} );

    [[nodiscard]] constexpr bool
    isValid() const noexcept
    {
        switch ( m_compressionType )
        {
        case CompressionType::RESERVED:
            return false;

        case CompressionType::UNCOMPRESSED:
            return true;

        case CompressionType::FIXED_HUFFMAN:
        #ifdef _MSC_VER
            return true;
        #else
            return m_fixedHC.isValid();
        #endif

        case CompressionType::DYNAMIC_HUFFMAN:
            return m_literalHC.isValid();
        }

        return false;
    }

    [[nodiscard]] constexpr const auto&
    precodeCL() const noexcept
    {
        return m_precodeCL;
    }

    [[nodiscard]] constexpr const auto&
    distanceAndLiteralCL() const noexcept
    {
        return m_literalCL;
    }

    void
    setTrackBackreferences( bool enable ) noexcept
    {
        m_trackBackreferences = enable;
    }

    [[nodiscard]] bool
    trackBackreferences() const noexcept
    {
        return m_trackBackreferences;
    }

    [[nodiscard]] const std::vector<Backreference>&
    backreferences() const noexcept
    {
        return m_backreferences;
    }

    /**
     * Reinitializes this block to behave basically as if default-constructed,
     * This avoids a generic reinitialization, e.g., by copying a default-constructed Block to it
     * because it might be more expensive than necessary for multi-stream gzips because it would zero the whole
     * 128 KiB decode buffer and all the 64 KiB DistanceHuffmanCoding buffer even though that is unnecessary.
     */
    void
    reset( const std::optional<VectorView<uint8_t> > initialWindow = {} )
    {
        m_uncompressedSize = 0;

        m_atEndOfBlock = false;
        m_atEndOfFile = false;

        m_isLastBlock = false;
        m_compressionType = CompressionType::RESERVED;
        m_padding = 0;

        m_windowPosition = 0;
        m_containsMarkerBytes = true;
        m_decodedBytes = 0;

        m_distanceToLastMarkerByte = 0;

        m_trackBackreferences = false;
        m_decodedBytesAtBlockStart = 0;
        m_backreferences.clear();

        if ( initialWindow ) {
            setInitialWindow( *initialWindow );
        } else {
            m_window16 = initializeMarkedWindowBuffer();
        }
    }

private:
    template<typename Window>
    forceinline void
    appendToWindow( Window&                     window,
                    typename Window::value_type decodedSymbol );

    template<typename Window>
    forceinline void
    appendToWindowUnsafe( Window&                     window,
                          typename Window::value_type decodedSymbol );

    template<typename Window>
    forceinline void
    resolveBackreference( Window&  window,
                          uint16_t distance,
                          uint16_t length,
                          size_t   nBytesRead );

    template<typename Window>
    [[nodiscard]] std::pair<size_t, Error>
    readInternal( gzip::BitReader& bitReader,
                  size_t           nMaxToDecode,
                  Window&          window );

    template<typename Window>
    [[nodiscard]] std::pair<size_t, Error>
    readInternalUncompressed( gzip::BitReader& bitReader,
                              Window&          window );

    template<typename Window,
             typename HuffmanCoding>
    [[nodiscard]] std::pair<size_t, Error>
    readInternalCompressed( gzip::BitReader&     bitReader,
                            size_t               nMaxToDecode,
                            Window&              window,
                            const HuffmanCoding& coding );

#if defined( LIBRAPIDARCHIVE_WITH_ISAL ) || defined( WITH_MULTI_CACHED_HUFFMAN_DECODER )
    template<typename Window>
    [[nodiscard]] std::pair<size_t, Error>
    readInternalCompressedMultiCached( gzip::BitReader&                    bitReader,
                                       size_t                              nMaxToDecode,
                                       Window&                             window,
                                       const LiteralOrLengthHuffmanCoding& coding );

#elif defined( WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER )

    template<typename Window>
    [[nodiscard]] std::pair<size_t, Error>
    readInternalCompressedSpecialized( gzip::BitReader&                    bitReader,
                                       size_t                              nMaxToDecode,
                                       Window&                             window,
                                       const LiteralOrLengthHuffmanCoding& coding );

#endif

    [[nodiscard]] forceinline std::pair<uint16_t, Error>
    getDistance( gzip::BitReader& bitReader ) const;

    /**
     * @param position The position in the window where the next byte would be appended. Similar to std::end();
     * @param size How many of the bits before @ref position are requested. @ref position - @ref size is begin.
     * @return the areas last written in the circular window buffer. Because of the circularity, two VectorViews
     *         are returned and both are non-empty in the case of the last written data wrapping around.
     */
    template<typename Window,
             typename Symbol = typename Window::value_type,
             typename View = VectorView<Symbol> >
    [[nodiscard]] static std::array<View, 2>
    lastBuffers( const Window& window,
                 size_t        position,
                 size_t        size )
    {
        if ( size > window.size() ) {
            throw std::invalid_argument( "Requested more bytes than fit in the buffer. Data is missing!" );
        }

        std::array<View, 2> result{};
        if ( size == 0 ) {
            return result;
        }

        /* Calculate wrapped around begin without unsigned underflow during the difference. */
        const auto begin = ( position + window.size() - ( size % window.size() ) ) % window.size();
        if ( begin < position ) {
            result[0] = View( window.data() + begin, position - begin );
            return result;
        }

        result[0] = View( window.data() + begin, window.size() - begin );  // up to end of window
        result[1] = View( window.data(), position );  // wrapped around part at start of window
        return result;
    }

private:
    /**
     * Size is max back-reference distance + max back-reference length to avoid the case of "in-place" updating
     * (overlapping input and output). Round up to power of two in the hope of making modulo faster...
     * Note that this buffer may be used for 16-bit half-decompressed data for when the initial window buffer is
     * unknown as well as for the case of the window buffer being known which only requires uint8_t.
     * For the former we need twice the size!
     * @note The buffer size should probably be a power of two or else I observed a slowdown probably because the
     *       circular buffer index modulo operation cannot be executed by a simple bitwise 'and' anymore.
     * @note 128 KiB is quite a lot of stack pressure. It actually leads to stack overflows on MacOS when creating
     *       multiple Block objects in the function call hierarchy such as in getUsedWindowSymbols!
     */
    using PreDecodedBuffer = std::array<uint16_t, 2 * MAX_WINDOW_SIZE>;
    using DecodedBuffer = WeakArray<std::uint8_t, PreDecodedBuffer().size() * sizeof( uint16_t ) / sizeof( uint8_t )>;

    /* The marker byte buffer does not have to fit an uncompressed block because larger uncompressed blocks will
     * trigger a conversion from PreDecodedBuffer to DecodedBuffer anyway. */
    static_assert( PreDecodedBuffer().size() * sizeof( uint16_t ) / sizeof( uint8_t ) >= MAX_UNCOMPRESSED_SIZE,
                   "Buffer should at least be able to fit one uncompressed block." );
    static_assert( std::min( PreDecodedBuffer().size(),
                             PreDecodedBuffer().size() * sizeof( uint16_t ) / sizeof( uint8_t ) )
                   >= MAX_WINDOW_SIZE + MAX_RUN_LENGTH,
                   "Buffers should at least be able to fit the back-reference window plus the maximum match length." );

private:
    /**
     * @note Making this constexpr or an immediately evaluated lambda expression to initialize the buffer,
     * increases release build compile time from 30s to 135s with GCC 11, and no measurable change (~80s) with
     * Clang 15! But it is important to avoid frequent recomputations for streams of many small fixed huffman
     * blocks, so use a function-local static variable, the best kind of static besides static constexpr where
     * possible.
     *
     * E.g. create a test file with this:
     * @verbatim
     * echo foo | gzip > many-small-streams.gz
     * for (( i=0; i < 1000; ++i )); do cat many-small-streams.gz >> 1000-24B-fixed-huffman-streams.gz; done
     * for (( i=0; i < 1000; ++i )); do cat 1000-24B-fixed-huffman-streams.gz >> 1M-24B-fixed-huffman-streams.gz; done
     * @endverbatim
     *
     * Comparison with default tools:
     * @verbatim
     * time igzip -d -c 1M-24B-fixed-huffman-streams.gz > /dev/null
     *   -> 0.471s 0.526s 0.471s 0.482s 0.518s
     * time gzip -d -c 1M-24B-fixed-huffman-streams.gz > /dev/null
     *   -> 5.566s 5.409s 5.431s 5.514s 5.470s
     * @endverbatim
     *
     * And benchmark with:
     * @verbatim
     * cmake -DLIBRAPIDARCHIVE_WITH_ISAL=OFF .. && cmake --build -- rapidgzip
     * for (( i=0; i < 5; ++i )); do
     *     src/tools/rapidgzip -v -d -o /dev/null 1M-24B-fixed-huffman-streams.gz 2>&1 | grep "Decompressed"
     * done
     * @endverbatim
     * Before:
     * @verbatim
     * Decompressed in total 4000000 B in:
     *     8.06588 s -> 0.495916 MB/s
     *     8.00692 s -> 0.499568 MB/s
     *     8.17901 s -> 0.489057 MB/s
     *     8.12103 s -> 0.492548 MB/s
     *     8.13215 s -> 0.491875 MB/s
     * @endverbatim
     * After introducing the static function-local variable and returning a reference:
     * @verbatim
     * Decompressed in total 4000000 B in:
     *     2.63232 s -> 1.51957 MB/s
     *     2.66496 s -> 1.50096 MB/s
     *     2.54554 s -> 1.57138 MB/s
     *     2.60184 s -> 1.53737 MB/s
     *     2.64525 s -> 1.51214 MB/s
     * @endverbatim
     * After avoiding unnecessary calls to this function during resetting of the block:
     * @verbatim
     * Decompressed in total 4000000 B in:
     *     0.550918 s -> 7.2606 MB/s
     *     0.561135 s -> 7.12841 MB/s
     *     0.572299 s -> 6.98935 MB/s
     *     0.592707 s -> 6.7487 MB/s
     *     0.606601 s -> 6.59412 MB/s
     * @endverbatim
     */
    [[nodiscard]] static const PreDecodedBuffer&
    initializeMarkedWindowBuffer()
    {
        static const PreDecodedBuffer markers =
            [] ()
            {
                PreDecodedBuffer result{};
                for ( size_t i = 0; i < MAX_WINDOW_SIZE; ++i ) {
                    result[result.size() - MAX_WINDOW_SIZE + i] = i + MAX_WINDOW_SIZE;
                }
                return result;
            } ();
        return markers;
    }

    [[nodiscard]] constexpr DecodedBuffer
    getWindow() noexcept
    {
        return DecodedBuffer{ reinterpret_cast<std::uint8_t*>( m_window16.data() ) };
    }

private:
    uint16_t m_uncompressedSize{ 0 };

private:
    /* These flags might get triggered by the read function. */
    mutable bool m_atEndOfBlock{ false };
    mutable bool m_atEndOfFile{ false };

    bool m_isLastBlock{ false };
    CompressionType m_compressionType{ CompressionType::RESERVED };
    /**
     * For UNCOMPRESSED blocks, this will hold the encountered padding, which probably is 0
     * but we might want to check that.
     */
    uint8_t m_padding{ 0 };

#ifndef _MSC_VER
    /**
     * Initializing m_fixedHC statically without constexpr leads to very weird problems when compiling with ASAN.
     * The code might be too complex and run into the static initialization order fiasco.
     * But having this static const is very important to get a 10-100x speedup for finding deflate blocks!
     * MSVC chokes, i.e., crashes, on this line. That's why on MSVC, m_fixedHC is a function-local static variable.
     * > fatal error C1001: An internal error has occurred in the compiler.
     */
    static constexpr FixedHuffmanCoding m_fixedHC = createFixedHC();
#endif
    LiteralOrLengthHuffmanCoding m_literalHC;

    DistanceHuffmanCoding m_distanceHC;

    alignas( 64 ) PreDecodedBuffer m_window16{ initializeMarkedWindowBuffer() };

public:
    /**
     * Points to the index of the next code to be written in @ref getWindow. I.e., can also be interpreted as
     * the size of getWindow (in the beginning as long as it does not wrap).
     */
    size_t m_windowPosition{ 0 };
    /**
     * If true, then @ref m_window16 should be used, else @ref getWindow!
     * When @ref m_distanceToLastMarkerByte reaches a sufficient threshold, @ref m_window16 will be converted
     * to @ref getWindow and this variable will be set to true.
     */
    bool m_containsMarkerBytes{ true };
    /**
     * Sum of decoded bytes over all read calls. Also will be set when calling setInitialWindow.
     * It is used to determine whether a backreference references valid data.
     */
    size_t m_decodedBytes{ 0 };

    /**
     * This is incremented whenever a symbol could be fully decoded and it gets reset when a marker byte is
     * encountered. It is used to determine when the last window buffer has been fully decoded.
     * The exact value does not matter and is undefined when @ref m_containsMarkerBytes is false.
     */
    size_t m_distanceToLastMarkerByte{ 0 };

    bool m_trackBackreferences{ false };
    size_t m_decodedBytesAtBlockStart{ 0 };
    std::vector<Backreference> m_backreferences;

    /* Large buffers required only temporarily inside readHeader. */
    alignas( 64 ) std::array<uint8_t, MAX_PRECODE_COUNT> m_precodeCL{};
    alignas( 64 ) PrecodeHuffmanCoding m_precodeHC{};
    alignas( 64 ) LiteralAndDistanceCLBuffer m_literalCL{};
};


template<bool ENABLE_STATISTICS>
template<bool treatLastBlockAsError>
Error
Block<ENABLE_STATISTICS>::readHeader( gzip::BitReader& bitReader )
{
    try {
        m_isLastBlock = bitReader.read<1>();
    } catch ( const gzip::BitReader::EndOfFileReached& ) {
        return Error::END_OF_FILE;
    }

    if constexpr ( treatLastBlockAsError ) {
        if ( m_isLastBlock ) {
            return Error::UNEXPECTED_LAST_BLOCK;
        }
    }
    m_compressionType = static_cast<CompressionType>( bitReader.read<2>() );

    Error error = Error::NONE;

    switch ( m_compressionType )
    {
    case CompressionType::UNCOMPRESSED:
    {
        /* @todo There is no mention what the padding is. But there is mention for the flags, that the reserved ones
         *       ones should be zero. Could I also check for the padding to be zero? I just don't want to believe,
         *       that anyone would store random data here ... Although it might be good for stenography :D */
        if ( bitReader.tell() % BYTE_SIZE != 0 ) {
            m_padding = static_cast<uint8_t>( bitReader.read( BYTE_SIZE - ( bitReader.tell() % BYTE_SIZE ) ) );
            if ( m_padding != 0 ) {
                return Error::NON_ZERO_PADDING;
            }
        }

        m_uncompressedSize = bitReader.read<2 * BYTE_SIZE>();
        const auto negatedLength = bitReader.read<2 * BYTE_SIZE>();
        if ( m_uncompressedSize != static_cast<uint16_t>( ~negatedLength ) ) {
            return Error::LENGTH_CHECKSUM_MISMATCH;
        }
        break;
    }

    case CompressionType::FIXED_HUFFMAN:
        break;

    case CompressionType::DYNAMIC_HUFFMAN:
        error = readDynamicHuffmanCoding( bitReader );
        break;

    case CompressionType::RESERVED:
        return Error::INVALID_COMPRESSION;
    };

    m_atEndOfBlock = false;
    m_decodedBytesAtBlockStart = m_decodedBytes;
    m_backreferences.clear();

    return error;
}


template<bool ENABLE_STATISTICS>
Error
Block<ENABLE_STATISTICS>::readDynamicHuffmanCoding( gzip::BitReader& bitReader )
{
    if constexpr ( ENABLE_STATISTICS ) {
        times.readDynamicStart = now();
    }

    /**
     * Huffman codings map variable length (bit) codes to symbols.
     * Huffman codings are given a as a tuple of code lengths, i.e., number of bits for Huffman code to use.
     * The elements of the tuple correspond to the elements of the ordered set of symbols, i.e., the alphabet.
     * For reading the block header it is important to understand that there are three different Huffmann condings
     * and also alphabets:
     *  - Alphabet L: the mixed alphabet containing 286 literals and lengths / instructions.
     *  - Alphabet D: contains distances in 30 different symbols / instructions.
     *  - Alphabet P: contains 19 different symbols / instructions for reconstructing the code length tuples.
     *                It is also called Precode and used to encode L and D! It itself is "encoded" a sequence of
     *                3-bit numbers for the bit lengths.
     *                This means, there can be no longer Huffman code than 7 for this, i.e., fits into a char.
     */

    const auto literalCodeCount = 257 + bitReader.read<5>();
    if ( literalCodeCount > MAX_LITERAL_OR_LENGTH_SYMBOLS ) {
        durations.readDynamicHeader += duration( times.readDynamicStart );
        return Error::EXCEEDED_LITERAL_RANGE;
    }
    const auto distanceCodeCount = 1 + bitReader.read<5>();
    if ( distanceCodeCount > MAX_DISTANCE_SYMBOL_COUNT ) {
        durations.readDynamicHeader += duration( times.readDynamicStart );
        return Error::EXCEEDED_DISTANCE_RANGE;
    }
    const auto codeLengthCount = 4 + bitReader.read<4>();

    if constexpr ( ENABLE_STATISTICS ) {
        this->precodeCLHistogram[codeLengthCount - 4]++;
        this->codeCounts.precode = codeLengthCount;
        this->codeCounts.distance = distanceCodeCount;
        this->codeCounts.literal = literalCodeCount;
    }

    /* Get code lengths (CL) for alphabet P. */
    std::memset( m_precodeCL.data(), 0, m_precodeCL.size() * sizeof( m_precodeCL[0] ) );
    for ( size_t i = 0; i < codeLengthCount; ++i ) {
        m_precodeCL[PRECODE_ALPHABET[i]] = bitReader.read<PRECODE_BITS>();
    }

    if constexpr ( ENABLE_STATISTICS ) {
        times.readPrecode = now();
        durations.readPrecode += duration( times.readDynamicStart, times.readPrecode );
    }

    auto error = m_precodeHC.initializeFromLengths( VectorView<uint8_t>( m_precodeCL.data(), m_precodeCL.size() ) );

    if constexpr ( ENABLE_STATISTICS ) {
        times.createdPrecodeHC = now();
        this->durations.createPrecodeHC += duration( times.readPrecode, times.createdPrecodeHC );
    }

    if ( error != Error::NONE ) {
        if constexpr ( ENABLE_STATISTICS ) {
            this->failedPrecodeInit++;
            durations.readDynamicHeader += duration( times.readDynamicStart );
        }
        return error;
    }

    /* Decode the code lengths for the literal/length and distance alphabets. */
    auto precodeApplyError = readDistanceAndLiteralCodeLengths(
        m_literalCL, bitReader, m_precodeHC, literalCodeCount + distanceCodeCount );

    if constexpr ( ENABLE_STATISTICS ) {
        times.appliedPrecodeHC = now();
        durations.applyPrecodeHC += duration( times.createdPrecodeHC, times.appliedPrecodeHC );
    }

    if ( precodeApplyError != Error::NONE ) {
        if constexpr ( ENABLE_STATISTICS ) {
            this->failedPrecodeApply++;
            durations.readDynamicHeader += duration( times.readDynamicStart );
        }
        return precodeApplyError;
    }

    /* Check for end-of-block symbol to have a non-zero code length. */
    if ( m_literalCL[deflate::END_OF_BLOCK_SYMBOL] == 0 ) {
        if constexpr ( ENABLE_STATISTICS ) {
            durations.readDynamicHeader += duration( times.readDynamicStart );
            this->missingEOBSymbol++;
        }
        return Error::INVALID_CODE_LENGTHS;
    }

    /* Create distance HC
     * When encoding base64-encoded random-data, I encountered a length of 9, so uint16_t is necessary! */
    error = m_distanceHC.initializeFromLengths(
        VectorView<uint8_t>( m_literalCL.data() + literalCodeCount, distanceCodeCount ) );

    if constexpr ( ENABLE_STATISTICS ) {
        times.createdDistanceHC = now();
        durations.createDistanceHC += duration( times.appliedPrecodeHC, times.createdDistanceHC );
    }

    if ( error != Error::NONE ) {
        if constexpr ( ENABLE_STATISTICS ) {
            durations.readDynamicHeader += duration( times.readDynamicStart );
            this->failedDistanceInit++;
        }
        return error;
    }

    /* Create literal HC */
#ifdef WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER
    error = m_literalHC.initializeFromLengths( VectorView<uint8_t>( m_literalCL.data(), literalCodeCount ),
                                               m_distanceHC );
#else
    error = m_literalHC.initializeFromLengths( VectorView<uint8_t>( m_literalCL.data(), literalCodeCount ) );
#endif
    if ( error != Error::NONE ) {
        if constexpr ( ENABLE_STATISTICS ) {
            this->failedLiteralInit++;
        }
    }

    if constexpr ( ENABLE_STATISTICS ) {
        const auto tFinish = now();
        durations.createLiteralHC += duration( times.createdDistanceHC, tFinish );
        durations.readDynamicHeader += duration( times.readDynamicStart, tFinish );
    }

    return error;
}


template<bool ENABLE_STATISTICS>
forceinline
std::pair<uint16_t, Error>
Block<ENABLE_STATISTICS>::getDistance( gzip::BitReader& bitReader ) const
{
    uint16_t distance = 0;
    if ( m_compressionType == CompressionType::FIXED_HUFFMAN ) {
        distance = reverseBits( static_cast<uint8_t>( bitReader.read<5>() ) ) >> 3U;
        if ( UNLIKELY( distance >= MAX_DISTANCE_SYMBOL_COUNT ) ) [[unlikely]] {
            return { 0, Error::EXCEEDED_DISTANCE_RANGE };
        }
    } else {
        const auto decodedDistance = m_distanceHC.decode( bitReader );
        if ( UNLIKELY( !decodedDistance ) ) [[unlikely]] {
            return { 0, Error::INVALID_HUFFMAN_CODE };
        }
        distance = static_cast<uint16_t>( *decodedDistance );
    }

    if ( distance <= 3U ) {
        distance += 1U;
    } else if ( distance <= 29U ) {
        const auto extraBitsCount = ( distance - 2U ) / 2U;
        const auto extraBits = bitReader.read( extraBitsCount );
        distance = distanceLUT[distance] + extraBits;
    } else {
        throw std::logic_error( "Invalid distance codes encountered!" );
    }

    return { distance, Error::NONE };
}


template<bool ENABLE_STATISTICS>
std::pair<DecodedDataView, Error>
Block<ENABLE_STATISTICS>::read( gzip::BitReader& bitReader,
                                size_t           nMaxToDecode )
{
    if ( eob() ) {
        return { {}, Error::NONE };
    }

    if ( m_compressionType == CompressionType::RESERVED ) {
        throw std::domain_error( "Invalid deflate compression type!" );
    }

    if constexpr ( ENABLE_STATISTICS ) {
        times.readDataStart = now();
    }

    DecodedDataView result;
    const auto window = getWindow();

    if ( m_compressionType == CompressionType::UNCOMPRESSED ) {
        std::optional<size_t> nBytesRead;
        if ( m_uncompressedSize >= MAX_WINDOW_SIZE ) {
            /* Special case for uncompressed blocks larger or equal than the window size.
             * Because, in this case, we can simply copy the uncompressed block to the beginning of the window
             * without worrying about wrap-around and also now we know that there are no markers remaining! */
            m_windowPosition = m_uncompressedSize;
            nBytesRead = bitReader.read( reinterpret_cast<char*>( window.data() ), m_uncompressedSize );
        } else if ( m_containsMarkerBytes && ( m_distanceToLastMarkerByte + m_uncompressedSize >= MAX_WINDOW_SIZE ) ) {
            /* Special case for when the new uncompressed data plus some fully-decoded data from the window buffer
             * together exceed the maximum backreference distance. */
            assert( m_distanceToLastMarkerByte <= m_decodedBytes );

            /* Copy and at the same time downcast enough data for the window from the 16-bit element buffer. */
            std::vector<uint8_t> remainingData( MAX_WINDOW_SIZE - m_uncompressedSize );
            size_t downcastedSize{ 0 };
            for ( const auto buffer : lastBuffers( m_window16, m_windowPosition, remainingData.size() ) ) {
                if ( std::any_of( buffer.begin(), buffer.end(), [] ( const auto symbol ) { return symbol > 255; } ) ) {
                    throw std::logic_error( "Encountered marker byte even though there shouldn't be one!" );
                }

                std::transform( buffer.begin(), buffer.end(), remainingData.data() + downcastedSize,
                                [] ( const auto symbol ) { return static_cast<uint8_t>( symbol ); } );
                downcastedSize += buffer.size();
            }

            m_windowPosition = MAX_WINDOW_SIZE;

            std::memcpy( window.data(), remainingData.data(), remainingData.size() );
            nBytesRead = bitReader.read( reinterpret_cast<char*>( window.data() + remainingData.size() ),
                                         m_uncompressedSize );
        } else if ( !m_containsMarkerBytes ) {
            /* When there are no markers, it means we can simply memcpy into the uint8_t window.
             * This speeds things up from ~400 MB/s to ~ 6 GB/s compared to calling appendToWindow for each byte!
             * We can use @ref lastBuffers, which are also returned, to determine where to copy to. */
            m_windowPosition = ( m_windowPosition + m_uncompressedSize ) % window.size();
            size_t totalBytesRead{ 0 };
            auto buffers =
                lastBuffers<DecodedBuffer, uint8_t, WeakVector<uint8_t> >(
                    window, m_windowPosition, m_uncompressedSize );
            for ( auto& buffer : buffers ) {
                totalBytesRead += bitReader.read( reinterpret_cast<char*>( buffer.data() ), buffer.size() );
            }
            nBytesRead = totalBytesRead;
        }

        if ( nBytesRead ) {
            m_containsMarkerBytes = false;
            m_atEndOfBlock = true;
            m_decodedBytes += *nBytesRead;

            result.data = lastBuffers( window, m_windowPosition, *nBytesRead );

            if constexpr ( ENABLE_STATISTICS ) {
                durations.readData += duration( times.readDataStart );
            }

            return { result, *nBytesRead == m_uncompressedSize ? Error::NONE : Error::EOF_UNCOMPRESSED };
        }
    }

    size_t nBytesRead{ 0 };
    auto error = Error::NONE;
    if ( m_containsMarkerBytes ) {
        /* This is the only case that should increment or reset m_distanceToLastMarkerByte. */
        std::tie( nBytesRead, error ) = readInternal( bitReader, nMaxToDecode, m_window16 );

        /* Theoretically, it would be enough if m_distanceToLastMarkerByte >= MAX_WINDOW_SIZE but that complicates
         * things because we can only convert up to m_distanceToLastMarkerByte of data even though we might need
         * to return up to nBytesRead of data! Furthermore, the wrap-around, again, would be more complicated. */
        if ( ( m_distanceToLastMarkerByte >= m_window16.size() )
             || ( ( m_distanceToLastMarkerByte >= MAX_WINDOW_SIZE )
                  && ( m_distanceToLastMarkerByte == m_decodedBytes ) ) ) {
            setInitialWindow();
            result.data = lastBuffers( window, m_windowPosition, nBytesRead );
        } else {
            result.dataWithMarkers = lastBuffers( m_window16, m_windowPosition, nBytesRead );
        }
    } else {
        std::tie( nBytesRead, error ) = readInternal( bitReader, nMaxToDecode, window );
        result.data = lastBuffers( window, m_windowPosition, nBytesRead );
    }

    if constexpr ( ENABLE_STATISTICS ) {
        durations.readData += duration( times.readDataStart );
    }

    return { result, error };
}


template<bool ENABLE_STATISTICS>
template<typename Window>
inline void
Block<ENABLE_STATISTICS>::appendToWindow( Window&                     window,
                                          typename Window::value_type decodedSymbol )
{
    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<typename Window::value_type>, uint16_t>;

    if constexpr ( containsMarkerBytes ) {
        if ( decodedSymbol > std::numeric_limits<uint8_t>::max() ) {
            m_distanceToLastMarkerByte = 0;
        } else {
            ++m_distanceToLastMarkerByte;
        }
    }

    window[m_windowPosition] = decodedSymbol;
    m_windowPosition++;
    m_windowPosition %= window.size();
}


template<bool ENABLE_STATISTICS>
template<typename Window>
inline void
Block<ENABLE_STATISTICS>::appendToWindowUnsafe( Window&                     window,
                                                typename Window::value_type decodedSymbol )
{
    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<typename Window::value_type>, uint16_t>;

    if constexpr ( containsMarkerBytes ) {
        if ( decodedSymbol > std::numeric_limits<uint8_t>::max() ) {
            m_distanceToLastMarkerByte = 0;
        } else {
            ++m_distanceToLastMarkerByte;
        }
    }

    window[m_windowPosition] = decodedSymbol;
    m_windowPosition++;
}


template<bool ENABLE_STATISTICS>
template<typename Window>
inline void
Block<ENABLE_STATISTICS>::resolveBackreference( Window&        window,
                                                const uint16_t distance,
                                                const uint16_t length,
                                                const size_t   nBytesRead )
{
    if ( m_trackBackreferences ) {
        if ( m_decodedBytes < m_decodedBytesAtBlockStart ) {
            throw std::logic_error( "Somehow the decoded bytes counter seems to have shrunk!" );
        }
        const auto decodedBytesInBlock = m_decodedBytes - m_decodedBytesAtBlockStart + nBytesRead;
        if ( distance > decodedBytesInBlock ) {
            m_backreferences.emplace_back( Backreference{
                static_cast<uint16_t>( distance - decodedBytesInBlock ),
                std::min( length, distance )
            } );
        }
    }

    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<decltype( *window.data() ) >, uint16_t>;

    const auto offset = ( m_windowPosition + window.size() - distance ) % window.size();
    const auto nToCopyPerRepeat = std::min( distance, length );
    assert( nToCopyPerRepeat != 0 );
    /* Note: NOT "<= window.size()" but only "<" because for equality I would have to
     *       compute modulo window.size() instead of simply: m_windowPosition += length. */
    if ( LIKELY( m_windowPosition + length < window.size() ) ) [[likely]] {
        if ( LIKELY( ( length <= distance ) && ( distance <= m_windowPosition ) ) ) [[likely]] {
            std::memcpy( &window[m_windowPosition], &window[offset], length * sizeof( window.front() ) );
            m_windowPosition += length;

            if constexpr ( containsMarkerBytes ) {
                size_t distanceToLastMarkerByte{ 0 };
                for ( ; distanceToLastMarkerByte < length; ++distanceToLastMarkerByte ) {
                    if ( window[m_windowPosition - 1 - distanceToLastMarkerByte]
                         > std::numeric_limits<uint8_t>::max() ) {
                        m_distanceToLastMarkerByte = distanceToLastMarkerByte;
                        return;
                    }
                }
                m_distanceToLastMarkerByte += length;
            }
            return;
        }

        if constexpr ( !containsMarkerBytes ) {
            if ( UNLIKELY( nToCopyPerRepeat == 1 ) ) [[unlikely]] {
                std::memset( &window[m_windowPosition], window[offset], length );
                m_windowPosition += length;
                return;
            }
        }

        for ( size_t nCopied = 0; nCopied < length; ) {
            for ( size_t position = offset;
                  ( position < offset + nToCopyPerRepeat ) && ( nCopied < length );
                  ++position, ++nCopied )
            {
                const auto copiedSymbol = window[position % window.size()];
                appendToWindowUnsafe( window, copiedSymbol );
            }
        }
        return;
    }

    for ( size_t nCopied = 0; nCopied < length; ) {
        for ( size_t position = offset;
              ( position < offset + nToCopyPerRepeat ) && ( nCopied < length );
              ++position, ++nCopied )
        {
            const auto copiedSymbol = window[position % window.size()];
            appendToWindow( window, copiedSymbol );
        }
    }
}


template<bool ENABLE_STATISTICS>
template<typename Window>
std::pair<size_t, Error>
Block<ENABLE_STATISTICS>::readInternal( gzip::BitReader& bitReader,
                                        size_t           nMaxToDecode,
                                        Window&          window )
{
    if ( m_compressionType == CompressionType::UNCOMPRESSED ) {
        /* This does not take into account nMaxToDecode to avoid additional state to keep track off. */
        return readInternalUncompressed( bitReader, window );
    }

    if ( m_compressionType == CompressionType::FIXED_HUFFMAN ) {
    #ifdef _MSC_VER
        /**
         * Initialization of a local static variable is thread-safe and happens on first pass as opposed to
         * the static initialization ordering fiasco for global or class-scope static variables.
         * @see https://stackoverflow.com/questions/8102125/is-local-static-variable-initialization-thread-safe-in-c11
         */
        static const auto fixedHC = createFixedHC();
        return readInternalCompressed( bitReader, nMaxToDecode, window, fixedHC );
    #else
        return readInternalCompressed( bitReader, nMaxToDecode, window, m_fixedHC );
    #endif
    }

#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    if constexpr ( std::is_same_v<LiteralOrLengthHuffmanCoding, HuffmanCodingISAL> ) {
        return readInternalCompressedMultiCached( bitReader, nMaxToDecode, window, m_literalHC );
    } else {
        return readInternalCompressed( bitReader, nMaxToDecode, window, m_literalHC );
    }
#elif defined( WITH_MULTI_CACHED_HUFFMAN_DECODER )
    return readInternalCompressedMultiCached( bitReader, nMaxToDecode, window, m_literalHC );
#elif defined( WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER )
    return readInternalCompressedSpecialized( bitReader, nMaxToDecode, window, m_literalHC );
#else
    return readInternalCompressed( bitReader, nMaxToDecode, window, m_literalHC );
#endif
}


template<bool ENABLE_STATISTICS>
template<typename Window>
std::pair<size_t, Error>
Block<ENABLE_STATISTICS>::readInternalUncompressed( gzip::BitReader& bitReader,
                                                    Window&          window )
{
    /**
     * Because the non-compressed deflate block size is 16-bit, the uncompressed data is limited to 65535 B!
     * The buffer can hold MAX_WINDOW_SIZE 16-bit values (for markers) or twice the amount of decoded bytes.
     * Therefore, this routine is safe to call in respect of "buffer overflows" before returning the view to
     * the buffer.
     *
     * Timings for different buffer sizes in MB/s for 2GiB-random.gz:
     * @verbatim
     *    8 B : 398.55  411.779 409.841
     *   16 B : 386.543 385.621 385.567
     *   32 B : 412.783 407.354 402.352 402.129
     *   64 B : 397.71  412.952 413.265 416.339
     *  128 B : 379.629 380.691 387.439 377.22
     *  256 B : 380.17  389.722 387.635 405.699
     *  512 B : 382.466 379.642 390.317 381.801
     * 1024 B : 384.92  386.544 381.748 388.71
     * 2048 B : 378.362 394.002 391.357 389.728
     * 4096 B : 380.87  379.09  386.711 395.955
     * @endverbatim
     */
    uint32_t totalBytesRead{ 0 };
    std::array<uint8_t, 64> buffer{};
    for ( ; totalBytesRead + buffer.size() <= m_uncompressedSize; totalBytesRead += buffer.size() ) {
        const auto nBytesRead = bitReader.read( reinterpret_cast<char*>( buffer.data() ), buffer.size() );
        for ( size_t i = 0; i < nBytesRead; ++i ) {
            appendToWindow( window, buffer[i] );
        }
    }
    for ( ; totalBytesRead < m_uncompressedSize; ++totalBytesRead ) {
        appendToWindow( window, static_cast<uint8_t>( bitReader.read<BYTE_SIZE>() ) );
    }
    m_atEndOfBlock = true;
    m_decodedBytes += m_uncompressedSize;
    return { m_uncompressedSize, Error::NONE };
}


template<bool ENABLE_STATISTICS>
template<typename Window,
         typename HuffmanCoding>
std::pair<size_t, Error>
Block<ENABLE_STATISTICS>::readInternalCompressed( gzip::BitReader&     bitReader,
                                                  size_t               nMaxToDecode,
                                                  Window&              window,
                                                  const HuffmanCoding& coding )
{
    if ( !coding.isValid() ) {
        throw std::invalid_argument( "No Huffman coding loaded! Call readHeader first!" );
    }

    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<decltype( *window.data() ) >, uint16_t>;

    nMaxToDecode = std::min( nMaxToDecode, window.size() - MAX_RUN_LENGTH );

    size_t nBytesRead{ 0 };
    for ( nBytesRead = 0; nBytesRead < nMaxToDecode; )
    {
        const auto decoded = coding.decode( bitReader );
        if ( !decoded ) {
            return { nBytesRead, Error::INVALID_HUFFMAN_CODE };
        }
        auto code = *decoded;

        if ( code <= 255 ) {
            if constexpr ( ENABLE_STATISTICS ) {
                symbolTypes.literal++;
            }

            appendToWindow( window, code );
            ++nBytesRead;
            continue;
        }

        if ( UNLIKELY( code == END_OF_BLOCK_SYMBOL /* 256 */ ) ) [[unlikely]] {
            m_atEndOfBlock = true;
            break;
        }

        if ( UNLIKELY( code > 285 ) ) [[unlikely]] {
            return { nBytesRead, Error::INVALID_HUFFMAN_CODE };
        }

        if constexpr ( ENABLE_STATISTICS ) {
            symbolTypes.backreference++;
        }

        const auto length = getLength( code, bitReader );
        if ( length != 0 ) {
            if constexpr ( ENABLE_STATISTICS ) {
                symbolTypes.copies += length;
            }
            const auto [distance, error] = getDistance( bitReader );
            if ( error != Error::NONE ) {
                return { nBytesRead, error };
            }

            if constexpr ( !containsMarkerBytes ) {
                if ( distance > m_decodedBytes + nBytesRead ) {
                    return { nBytesRead, Error::EXCEEDED_WINDOW_RANGE };
                }
            }

            resolveBackreference( window, distance, length, nBytesRead );
            nBytesRead += length;
        }
    }

    m_decodedBytes += nBytesRead;
    return { nBytesRead, Error::NONE };
}


#if defined( LIBRAPIDARCHIVE_WITH_ISAL ) || defined( WITH_MULTI_CACHED_HUFFMAN_DECODER )
template<bool ENABLE_STATISTICS>
template<typename Window>
std::pair<size_t, Error>
Block<ENABLE_STATISTICS>::readInternalCompressedMultiCached
(
    gzip::BitReader&                    bitReader,
    size_t                              nMaxToDecode,
    Window&                             window,
    const LiteralOrLengthHuffmanCoding& coding )
{
    if ( !coding.isValid() ) {
        throw std::invalid_argument( "No Huffman coding loaded! Call readHeader first!" );
    }

    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<decltype( *window.data() ) >, uint16_t>;

    nMaxToDecode = std::min( nMaxToDecode, window.size() - MAX_RUN_LENGTH );

    size_t nBytesRead{ 0 };
    for ( nBytesRead = 0; nBytesRead < nMaxToDecode; )
    {
        auto [symbol, symbolCount] = coding.decode( bitReader );
        if ( symbolCount == 0 ) {
            return { nBytesRead, Error::INVALID_HUFFMAN_CODE };
        }

        for ( ; symbolCount > 0; symbolCount--, symbol >>= 8U ) {
            const auto code = static_cast<uint16_t>( symbol & 0xFFFFU );

            if ( ( code <= 255 ) || ( symbolCount > 1 ) ) {
                if constexpr ( ENABLE_STATISTICS ) {
                    symbolTypes.literal++;
                }

                appendToWindow( window, static_cast<uint8_t>( code ) );
                ++nBytesRead;
                continue;
            }

            if ( UNLIKELY( code == END_OF_BLOCK_SYMBOL /* 256 */ ) ) [[unlikely]] {
                m_atEndOfBlock = true;
                m_decodedBytes += nBytesRead;
                return { nBytesRead, Error::NONE };
            }

            static constexpr auto MAX_LIT_LEN_SYM = 512U;
            if ( UNLIKELY( code > MAX_LIT_LEN_SYM ) ) [[unlikely]] {
                return { nBytesRead, Error::INVALID_HUFFMAN_CODE };
            }

            if constexpr ( ENABLE_STATISTICS ) {
                symbolTypes.backreference++;
            }

            /* If the next symbol is a repeat length, read in the length extra bits, the distance code, the distance
             * extra bits. Then write out the corresponding data and update the state data accordingly. */
            const auto length = symbol - 254U;
            if ( length != 0 ) {
                if constexpr ( ENABLE_STATISTICS ) {
                    symbolTypes.copies += length;
                }
                const auto [distance, error] = getDistance( bitReader );
                if ( error != Error::NONE ) {
                    return { nBytesRead, error };
                }

                if constexpr ( !containsMarkerBytes ) {
                    if ( distance > m_decodedBytes + nBytesRead ) {
                        return { nBytesRead, Error::EXCEEDED_WINDOW_RANGE };
                    }
                }

                resolveBackreference( window, distance, length, nBytesRead );
                nBytesRead += length;
            }
        }
    }

    m_decodedBytes += nBytesRead;
    return { nBytesRead, Error::NONE };
}

#elif defined( WITH_DEFLATE_SPECIFIC_HUFFMAN_DECODER )

template<bool ENABLE_STATISTICS>
template<typename Window>
std::pair<size_t, Error>
Block<ENABLE_STATISTICS>::readInternalCompressedSpecialized
(
    gzip::BitReader&                    bitReader,
    size_t                              nMaxToDecode,
    Window&                             window,
    const LiteralOrLengthHuffmanCoding& coding )
{
    if ( !coding.isValid() ) {
        throw std::invalid_argument( "No Huffman coding loaded! Call readHeader first!" );
    }

    constexpr bool containsMarkerBytes = std::is_same_v<std::decay_t<decltype( *window.data() ) >, uint16_t>;

    nMaxToDecode = std::min( nMaxToDecode, window.size() - MAX_RUN_LENGTH );

    size_t nBytesRead{ 0 };
    LiteralOrLengthHuffmanCoding::CacheEntry cacheEntry;
    for ( nBytesRead = 0; nBytesRead < nMaxToDecode; )
    {
        try {
            cacheEntry = coding.decode( bitReader, m_distanceHC );
        } catch ( const Error& errorCode ) {
            return { nBytesRead, errorCode };
        }

        switch ( cacheEntry.distance )
        {
        case 0xFFFFU:
            m_atEndOfBlock = true;
            m_decodedBytes += nBytesRead;
            return { nBytesRead, Error::NONE };

        case 0U:
            if constexpr ( ENABLE_STATISTICS ) {
                symbolTypes.literal++;
            }
            appendToWindow( window, cacheEntry.symbolOrLength );
            ++nBytesRead;
            break;

        default:
        {
            const auto length = cacheEntry.symbolOrLength + 3U;
            if constexpr ( ENABLE_STATISTICS ) {
                symbolTypes.backreference++;
                symbolTypes.copies += length;
            }

            if constexpr ( !containsMarkerBytes ) {
                if ( cacheEntry.distance > m_decodedBytes + nBytesRead ) {
                    return { nBytesRead, Error::EXCEEDED_WINDOW_RANGE };
                }
            }

            resolveBackreference( window, cacheEntry.distance, length, nBytesRead );
            nBytesRead += length;
            break;
        }
        }
    }

    m_decodedBytes += nBytesRead;
    return { nBytesRead, Error::NONE };
}
#endif  // ifdef LIBRAPIDARCHIVE_WITH_ISAL


template<bool ENABLE_STATISTICS>
void
Block<ENABLE_STATISTICS>::setInitialWindow( VectorView<uint8_t> const& initialWindow )
{
    if ( !m_containsMarkerBytes ) {
        return;
    }

    const auto window = getWindow();

    /* Set an initial window before decoding has started. */
    if ( ( m_decodedBytes == 0 ) && ( m_windowPosition == 0 ) ) {
        if ( !initialWindow.empty() ) {
            std::memcpy( window.data(), initialWindow.data(), initialWindow.size() );
            m_windowPosition = initialWindow.size();
            m_decodedBytes = initialWindow.size();
        }
        m_containsMarkerBytes = false;
        return;
    }

    /* The buffer is initialized with markers! We need to take care that we do not try to replace those. */
    for ( size_t i = 0; m_decodedBytes + i < m_window16.size(); ++i ) {
        m_window16[( m_windowPosition + i ) % m_window16.size()] = 0;
    }
    replaceMarkerBytes( { m_window16.data(), m_window16.size() }, initialWindow );

    /* We cannot simply move each byte to m_window because it has twice as many elements as m_windows16
     * and simply filling it from left to right will result in wrapping not working because the right half
     * is empty. It would only work if there is no wrapping necessary because it is a contiguous block!
     * To achieve that, map i -> i' such that m_windowPosition is m_window.size() - 1.
     * This way all back-references will not wrap around on the left border. */
    std::array<uint8_t, decltype( m_window16 )().size()> conflatedBuffer{};

    for ( size_t i = 0; i < m_window16.size(); ++i ) {
        conflatedBuffer[i] = m_window16[( i + m_windowPosition ) % m_window16.size()];
    }

    std::memcpy( window.data() + ( window.size() - conflatedBuffer.size() ),
                 conflatedBuffer.data(),
                 conflatedBuffer.size() );

    m_windowPosition = 0;

    m_containsMarkerBytes = false;
}


[[nodiscard]] inline bool
verifySparseWindow( gzip::BitReader&          bitReader,
                    const std::vector<bool>&  windowByteIsRequired,
                    const VectorView<uint8_t> expectedOutput )
{
    Block</* ENABLE_STATISTICS */ false> block;

    /* Check that the created sparse window is correct by setting all sparse bytes to some arbitrary canary token.
     * If not (should not happen), show a warning and fall back to a non-sparse window. */
    std::vector<uint8_t> initialWindow( MAX_WINDOW_SIZE, 0 );
    for ( size_t i = 0; i < windowByteIsRequired.size(); ++i ) {
        if ( !windowByteIsRequired[i] ) {
            initialWindow[i] = 1;
        }
    }
    block.setInitialWindow( { initialWindow.data(), initialWindow.size() } );

    for ( size_t nBytesRead = 0; nBytesRead < MAX_WINDOW_SIZE; ) {
        const auto headerError = block.readHeader( bitReader );
        if ( headerError == Error::END_OF_FILE ) {
            break;
        }
        if ( headerError != Error::NONE ) {
            throw std::invalid_argument( "Failed to decode the deflate block header! " + toString( headerError ) );
        }

        size_t nBytesReadFromBlock{ 0 };
        while ( ( nBytesRead + nBytesReadFromBlock < MAX_WINDOW_SIZE ) && !block.eob() ) {
            const auto [view, readError] = block.read( bitReader, MAX_WINDOW_SIZE - nBytesRead );
            if ( readError != Error::NONE ) {
                throw std::invalid_argument( "Failed to read deflate block data! " + toString( readError ) );
            }

            if ( view.dataWithMarkersSize() > 0 ) {
                throw std::logic_error( "Result should not contain markers because we have set a window!" );
            }
            for ( const auto& buffer : view.data ) {
                const auto sizeToCompare = std::min( expectedOutput.size() - ( nBytesRead + nBytesReadFromBlock ),
                                                     buffer.size() );
                if ( !std::equal( buffer.data(), buffer.data() + sizeToCompare,
                                  expectedOutput.data() + nBytesRead + nBytesReadFromBlock ) )
                {
                    return false;
                }
                nBytesReadFromBlock += buffer.size();
            }
        }

        nBytesRead += nBytesReadFromBlock;
        if ( block.eos() ) {
            break;
        }
    }

    return true;
}


[[nodiscard]] inline std::vector<bool>
getUsedWindowSymbols( gzip::BitReader& bitReader )
{
    std::vector<bool> window( MAX_WINDOW_SIZE, false );

    static constexpr bool CHECK_CORRECTNESS{ true };
    /* Store the decompressed data to check for sparsity correctness. Initialize to 0 for correctness and also
     * in order to set a dummy window with only zeros so that we do not get any marker bytes! This simplifies
     * the correctness check and assuming that the sparse bytes are cleared to 0, then using zeros for the
     * dummy window is perfect because it will not give mismatches in the case some sparsity is wrong but the
     * data to be sparsed out is 0 anyway! */
    [[maybe_unused]] std::vector<uint8_t> decompressed;
    [[maybe_unused]] const auto oldOffset = bitReader.tell();
    size_t nBytesRead{ 0 };

    /* Anonymous namespace to ensure that the 208kB deflate::Block is cleared from the stack before
     * allocating another such block inside verifySparseWindow! */
    {
        /* Size of deflate::Block is 207616. Allocation on stack did result in a SIGBUS
         * because of a stack overflow on MacOS:
         * #0 rapidgzip::deflate::Block<false>::setInitialWindow
         * #1 rapidgzip::deflate::getUsedWindowSymbols
         * #2 rapidgzip::GzipChunk<rapidgzip::ChunkData>::determineUsedWindowSymbols
         * #3 rapidgzip::GzipChunk<rapidgzip::ChunkData>::appendDeflateBlockBoundary
         * #4 rapidgzip::GzipChunk<rapidgzip::ChunkData>::decodeChunkWithRapidgzip
         * #5 rapidgzip::GzipChunk<rapidgzip::ChunkData>::decodeChunk
         * #6 rapidgzip::GzipChunkFetcher<FetchingStrategy::FetchMultiStream, rapidgzip::ChunkData>::decodeBlock
         * #7 rapidgzip::GzipChunkFetcher<FetchingStrategy::FetchMultiStream, rapidgzip::ChunkData>::decodeBlock
         * #8 BlockFetcher<...>::decodeAndMeasureBlock
         * #9 BlockFetcher<...>::submitOnDemandTask
         * #10 decltype(std::declval<BlockFetcher<rapidgzip::GzipBlockFinder, rapidgzip::ChunkData, FetchingStrategy::FetchMultiStream>::submitOnDemandTask
         * #11 std::__1::__packaged_task_func<BlockFetcher<rapidgzip::GzipBlockFinder, rapidgzip::ChunkData, FetchingStrategy::FetchMultiStream>::submitOnDemandTask
         * #12 std::__1::__packaged_task_function<rapidgzip::ChunkData ()>::operator()() const future:1880
         * #13 std::__1::packaged_task<rapidgzip::ChunkData ()>::operator()() future:1957
         * #14 ThreadPool::PackagedTaskWrapper::SpecializedFunctor<...>::operator()()
         * #15 ThreadPool::PackagedTaskWrapper::operator()()
         * #16 ThreadPool::workerMain(unsigned long)
         * [...]
         * HOWEVER, allocation on the heap with std::make_unique leads to some weird memory leak with
         * benchmarkIndexCompression. I cannot debug it with heaptrack because it shows only 1 GB memory usage,
         * not the 50+ GB observerd. It seems that POSIX brk is used and this leads to some fragmentation for
         * because it basically is only a linear allocator and can only free in LIFO order. Why exactly this
         * happens, I have no idea. It may be a tiny memory leak or it may be some weird allocation ordering
         * problem, especially with the heap-allocated return result ... The memory rises really fast with Clang 15,
         * but also seems to leak, albeit much much slower, with GCC 11.
         * It is made even harder to debug because this memory leak disappears, when enabling -fsanitize=leak, or
         * -fsanitize=address, and even when using valgrind --tool=massif! Utter fucking garbage!
         */
        Block</* ENABLE_STATISTICS */ false> block;
        block.setTrackBackreferences( true );

        if constexpr ( CHECK_CORRECTNESS ) {
            decompressed.assign( MAX_WINDOW_SIZE, 0 );
            block.setInitialWindow( { decompressed.data(), decompressed.size() } );
        }

        for ( ; nBytesRead < MAX_WINDOW_SIZE; ) {
            /* Block::readHeader also clears the backreferences. */
            const auto headerError = block.readHeader( bitReader );
            if ( headerError == Error::END_OF_FILE ) {
                break;
            }
            if ( headerError != Error::NONE ) {
                throw std::invalid_argument( "Failed to decode the deflate block header! " + toString( headerError ) );
            }

            size_t nBytesReadFromBlock{ 0 };
            while ( ( nBytesRead + nBytesReadFromBlock < MAX_WINDOW_SIZE ) && !block.eob() ) {
                const auto [view, readError] = block.read( bitReader, MAX_WINDOW_SIZE - nBytesRead );
                if ( readError != Error::NONE ) {
                    throw std::invalid_argument( "Failed to read deflate block data! " + toString( readError ) );
                }

                if constexpr ( CHECK_CORRECTNESS ) {
                    if ( view.dataWithMarkersSize() > 0 ) {
                        throw std::logic_error( "Result should not contain markers because we have set a window!" );
                    }
                    for ( const auto& buffer : view.data ) {
                        const auto sizeToCopy = std::min( decompressed.size() - ( nBytesRead + nBytesReadFromBlock ),
                                                          buffer.size() );
                        if ( sizeToCopy == 0 ) {
                            continue;
                        }
                        std::memcpy( decompressed.data() + nBytesRead + nBytesReadFromBlock, buffer.data(),
                                     sizeToCopy );
                        nBytesReadFromBlock += sizeToCopy;
                    }
                } else {
                    nBytesReadFromBlock += view.size();
                }
            }

            const auto& backreferences = block.backreferences();
            for ( const auto& reference : backreferences ) {
                /* The back-references are relative to the current block, so we need to subtract nBytesRead
                 * from the relative distance to get the distance relative to the first block start.
                 * If the result would become negative, then nothing from the window is needed and we can skip it. */
                if ( reference.distance < nBytesRead ) {
                    continue;
                }

                const auto distanceFromEnd = reference.distance - nBytesRead;
                if ( distanceFromEnd > window.size() ) {
                    std::stringstream message;
                    message << "The back-reference distance should not exceed MAX_WINDOW_SIZE ("
                            << formatBytes( MAX_WINDOW_SIZE ) << ") but got: " << formatBytes( distanceFromEnd ) << "!";
                    throw std::logic_error( std::move( message ).str() );
                }
                if ( reference.length == 0 ) {
                    continue;
                }
                const auto startOffset = window.size() - distanceFromEnd;

                for ( size_t i = 0; ( i < reference.length ) && ( startOffset + i < window.size() ); ++i ) {
                    window[startOffset + i] = true;
                }
            }

            nBytesRead += nBytesReadFromBlock;
            if ( block.eos() ) {
                break;
            }
        }
    }

    if constexpr ( CHECK_CORRECTNESS ) {
        bitReader.seekTo( oldOffset );
        if ( !verifySparseWindow( bitReader, window, { decompressed.data(), nBytesRead } ) ) {
            std::stringstream message;
            message << "[Warning] Sparse window detection failed at offset " << formatBits( oldOffset )
                    << ". Will fall back to full window\n";
        #ifdef RAPIDGZIP_FATAL_PERFORMANCE_WARNINGS
            throw std::logic_error( std::move( message ).str() );
        #else
            std::cerr << std::move( message ).str();
        #endif
            window.assign( MAX_WINDOW_SIZE, true );
            return window;
        }
    }

    return window;
}


template<typename Container>
[[nodiscard]] std::vector<uint8_t>
getSparseWindow( gzip::BitReader& bitReader,
                 const Container& window )
{
    const auto usedSymbols = getUsedWindowSymbols( bitReader );
    std::vector<uint8_t> sparseWindow( std::min<size_t>( 32_Ki, window.size() ), 0 );
    for ( size_t i = 0; i < sparseWindow.size(); ++i ) {
        if ( usedSymbols[i + ( usedSymbols.size() - sparseWindow.size() )] ) {
            sparseWindow[i] = window[i + ( window.size() - sparseWindow.size() )];
        }
    };
    return sparseWindow;
}
}  // namespace rapidgzip::deflate
