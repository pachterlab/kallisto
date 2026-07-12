#pragma once

/**
 * @file Contains some compile-time computation of magic constants for the precode huffman codings.
 *       This whole file is only used in tests, benchmarks, and SingleLUT.hpp, which is now also only used in
 *       tests and benchmarks. I.e., this file should not increase normal compile-times or binary size.
 */

#include <algorithm>
#include <array>
#include <cstdint>

#include <rapidgzip/huffman/HuffmanCodingReversedBitsCachedCompressed.hpp>

#include "definitions.hpp"


namespace rapidgzip::deflate::precode
{
constexpr auto MAX_DEPTH = MAX_PRECODE_LENGTH;

/* Contains how often a the code lengths [1,7] appear. */
using Histogram = std::array<uint8_t, MAX_DEPTH>;


/* operator== is only constexpr for std::array since C++20 */
[[nodiscard]] constexpr bool
operator==( const Histogram& a,
            const Histogram& b )
{
    for ( size_t i = 0; i < a.size(); ++i ) {
        if ( a[i] != b[i] ) {
            return false;
        }
    }
    return true;
}


/**
 * This alternative version tries to reduce the number of instructions required for creation so that it can be
 * used with MSVC, which is the worst in evaluating @ref createPrecodeFrequenciesValidLUT constexpr and runs into
 * internal compiler errors or out-of-heap-space errors.
 * This alternative implementation uses the fact only very few of the LUT entries are actually valid, meaning
 * we can reduce instructions by initializing to invalid and then iterating only over the valid possibilities.
 * @param depth A depth of 1 means that we should iterate over 1-bit codes, which can only be 0,1,2.
 * @param freeBits This can be calculated from the histogram but it saves constexpr instructions when
 *        the caller updates this value outside.
 */
template<typename Result,
         typename Functor,
         uint8_t  DEPTH = 1>
constexpr void
iterateValidPrecodeHistograms( Result&        result,
                               const Functor& processValidHistogram,
                               const uint32_t remainingCount = MAX_PRECODE_COUNT,
                               Histogram      histogram = Histogram{},
                               const uint32_t freeBits = 2 )
{
    static_assert( DEPTH <= MAX_DEPTH, "Cannot descend deeper than the frequency counts!" );

    /* The for loop maximum is given by the invalid Huffman code check, i.e.,
     * when there are more code lengths on a tree level than there are nodes. */
    for ( uint32_t count = 0; count <= std::min( remainingCount, freeBits ); ++count ) {
        histogram[DEPTH - 1] = count;
        const auto newFreeBits = ( freeBits - count ) * 2;

        /* The first layer may not be fully filled or even empty. This does not fit any of the general tests. */
        if constexpr ( DEPTH == 1 ) {
            if ( count == 1 ) {
                processValidHistogram( histogram );
            }
        }

        if constexpr ( DEPTH == MAX_DEPTH ) {
            if ( newFreeBits == 0 ) {
                processValidHistogram( histogram );
            }
        } else {
            if ( count == freeBits ) {
                processValidHistogram( histogram );
            } else {
                const auto newRemainingCount = remainingCount - count;
                iterateValidPrecodeHistograms<Result, Functor, DEPTH + 1>( result, processValidHistogram,
                                                                           newRemainingCount, histogram, newFreeBits );
            }
        }
    }
}


constexpr auto VALID_HISTOGRAMS_COUNT =
    [] ()
    {
        size_t validCount{ 0 };
        iterateValidPrecodeHistograms( validCount, [&validCount] ( const auto& ) { ++validCount; } );
        return validCount;
    }();

static_assert( VALID_HISTOGRAMS_COUNT == 1526 );


/* Size: sizeof( std::array<uint8_t, MAX_DEPTH = 7> ) * 1526 = 10.682 kB */
static constexpr auto VALID_HISTOGRAMS =
    [] ()
    {
        size_t validCount{ 0 };
        std::array<Histogram, VALID_HISTOGRAMS_COUNT> validHistograms{};

        iterateValidPrecodeHistograms(
            validHistograms,
            [&validCount, &validHistograms] ( const Histogram& histogram ) {
                validHistograms[validCount++] = histogram;
            } );

        return validHistograms;
    }();

static_assert( VALID_HISTOGRAMS.back() == Histogram{ { /* code length 1 */ 2, } } );
}  // namespace rapidgzip::deflate::precode
