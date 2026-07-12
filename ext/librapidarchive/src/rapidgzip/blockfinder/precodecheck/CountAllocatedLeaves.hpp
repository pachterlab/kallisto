/**
 * @file This method counts the number of allocated leaves for each code length. E.g., a code length equal to the
 *       maximum code length takes up 1 leave, while a code length 1 shorter takes up 2 leaves, and a code length
 *       2 shorter takes up 4 leaves, and so on. Except for the case with exactly 1 code length of 1, the full
 *       tree should be occupied in terms of projected leaves.
 *       This version is so easy and does not require huge lookup tables that it is embarrassing that I did not think
 *       of this approach much sooner. I was too occupied to incrementally improve the histogram based approach.
 */

#pragma once

#include <array>
#include <cstdint>

#include <core/BitManipulation.hpp>
#include <core/Error.hpp>
#include <core/common.hpp>
#include <rapidgzip/gzip/definitions.hpp>


namespace rapidgzip::PrecodeCheck::CountAllocatedLeaves
{
using LeafCount = uint16_t;


[[nodiscard]] constexpr LeafCount
getVirtualLeafCount( const uint64_t codeLength )
{
    /* A code length of 1 takes up  64 leaves at maximum depth of 7.
     * A code length of 2 takes up  32 leaves at maximum depth of 7.
     * A code length of 7 takes up   1 leaves at maximum depth of 7. */
    return codeLength > 0 ? 1U << ( deflate::MAX_PRECODE_LENGTH - codeLength ) : 0U;
}


[[nodiscard]] constexpr LeafCount
getVirtualLeafCount( const uint64_t precodeBits,
                     const size_t   codeLengthCount )
{
    size_t virtualLeafCount{ 0 };
    for ( size_t i = 0; i < codeLengthCount; ++i ) {
        virtualLeafCount += getVirtualLeafCount( ( precodeBits >> ( i * deflate::PRECODE_BITS ) )
                                                 & nLowestBitsSet<uint64_t, deflate::PRECODE_BITS>() );
    }
    return virtualLeafCount;
}


template<uint8_t VALUE_BITS,
         uint8_t VALUE_COUNT>
[[nodiscard]] constexpr LeafCount
computeLeafCount( uint64_t values )
{
    LeafCount result{ 0 };
    for ( uint8_t i = 0; i < VALUE_COUNT; ++i ) {
        result += getVirtualLeafCount(
            ( values >> ( i * VALUE_BITS ) ) & nLowestBitsSet<LeafCount, VALUE_BITS>() );
    }
    return result;
}


template<uint8_t PRECODE_CHUNK_SIZE>
static constexpr auto PRECODE_TO_LEAF_COUNT_LUT =
    [] ()
    {
        std::array<LeafCount, 1ULL << uint16_t( PRECODE_CHUNK_SIZE * deflate::PRECODE_BITS )> result{};
        for ( size_t i = 0; i < result.size(); ++i ) {
            result[i] = computeLeafCount<deflate::PRECODE_BITS, PRECODE_CHUNK_SIZE>( i );
        }
        return result;
    }();


template<uint8_t PRECODE_CHUNK_SIZE>
[[nodiscard]] constexpr forceinline LeafCount
precodesToLeafCount( uint64_t precodeBits )
{
    constexpr auto CACHED_BITS = deflate::PRECODE_BITS * PRECODE_CHUNK_SIZE;
    constexpr auto CHUNK_COUNT = ceilDiv( deflate::MAX_PRECODE_COUNT, PRECODE_CHUNK_SIZE );
    LeafCount histogram{ 0 };
    for ( size_t chunk = 0; chunk < CHUNK_COUNT; ++chunk ) {
        auto precodeChunk = precodeBits >> ( chunk * CACHED_BITS );
        /* The last requires no bit masking because @ref next57Bits is already sufficiently masked.
         * This branch will hopefully get unrolled, else it could hinder performance. */
        if ( chunk != CHUNK_COUNT - 1 ) {
            precodeChunk &= nLowestBitsSet<uint64_t, CACHED_BITS>();
        }
        histogram += PRECODE_TO_LEAF_COUNT_LUT<PRECODE_CHUNK_SIZE>[precodeChunk];
    }
    return histogram;
}


[[nodiscard]] constexpr Error
checkPrecode( const uint64_t next4Bits,
              const uint64_t next57Bits )
{
    const auto codeLengthCount = 4 + next4Bits;
    constexpr auto PRECODE_CHUNK_SIZE = 4U;
#if 0
    const auto precodeBits = next57Bits & nLowestBitsSet<uint64_t>( codeLengthCount * deflate::PRECODE_BITS );

    constexpr auto CACHED_BITS = deflate::PRECODE_BITS * PRECODE_CHUNK_SIZE;
    const auto CHUNK_COUNT = ceilDiv( codeLengthCount, PRECODE_CHUNK_SIZE );
    LeafCount virtualLeafCount{ 0 };
    /* Non-fixed CHUNK_COUNT based on codeLengthCount makes it much slower: 95 -> 67 MB/s */
    for ( size_t chunk = 0; chunk < CHUNK_COUNT; ++chunk ) {
        auto precodeChunk = precodeBits >> ( chunk * CACHED_BITS );
        /* The last requires no bit masking because @ref next57Bits is already sufficiently masked.
         * This branch will hopefully get unrolled, else it could hinder performance. */
        if ( chunk != CHUNK_COUNT - 1 ) {
            precodeChunk &= nLowestBitsSet<uint64_t, CACHED_BITS>();
        }
        virtualLeafCount += PRECODE_TO_LEAF_COUNT_LUT<PRECODE_CHUNK_SIZE>[precodeChunk];
    }
#elif 1
    /* Try manual Duff's device loop unrolling making use of the fact that 4 <= codeLengthCount <= 19 */
    constexpr auto& LUT = PRECODE_TO_LEAF_COUNT_LUT<PRECODE_CHUNK_SIZE>;
    constexpr auto CACHED_BITS = deflate::PRECODE_BITS * PRECODE_CHUNK_SIZE;
    const auto precodeBits = next57Bits & nLowestBitsSet<uint64_t>( codeLengthCount * deflate::PRECODE_BITS );

    LeafCount virtualLeafCount{ 0 };
    virtualLeafCount += LUT[precodeBits & nLowestBitsSet<uint64_t, CACHED_BITS>()];
    virtualLeafCount += LUT[( precodeBits >> CACHED_BITS ) & nLowestBitsSet<uint64_t, CACHED_BITS>()];
    /* Adding if ( codeLengthCount <= 12 ) here makes it slower: 95 -> 80 MB/s */
    virtualLeafCount += LUT[( precodeBits >> ( 2U * CACHED_BITS ) ) & nLowestBitsSet<uint64_t, CACHED_BITS>()];
    virtualLeafCount += LUT[( precodeBits >> ( 3U * CACHED_BITS ) ) & nLowestBitsSet<uint64_t, CACHED_BITS>()];
   /* The last requires no bit masking because BitReader::read already has masked to the lowest 57 bits
    * and this shifts 48 bits to the right, leaving only 9 (<12) bits set anyways. */
    virtualLeafCount += LUT[precodeBits >> ( 4U * CACHED_BITS )];

#elif 0
    /**
     * [13 bits, Single LUT] ( 90.2 <= 94.8 +- 2.2 <= 97.0 ) MB/s (candidates: 495)
     * [14 bits, Single LUT] ( 89.5 <= 91.0 +- 0.7 <= 91.8 ) MB/s (candidates: 495)
     * [15 bits, Single LUT] ( 86.5 <= 87.6 +- 0.9 <= 88.8 ) MB/s (candidates: 495)
     * [16 bits, Single LUT] ( 88.3 <= 89.3 +- 0.8 <= 90.5 ) MB/s (candidates: 495)
     * [17 bits, Single LUT] ( 87.5 <= 88.4 +- 0.7 <= 89.0 ) MB/s (candidates: 495)
     * [18 bits, Single LUT] ( 82.0 <= 84.5 +- 1.1 <= 85.3 ) MB/s (candidates: 495)
     *
     * [13 bits, this with chunk count 1] ( 69.0 <= 71.0 +- 1.3 <= 73.2 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 1] ( 70.4 <= 71.2 +- 0.6 <= 72.4 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 1] ( 67.9 <= 71.2 +- 1.6 <= 73.3 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 1] ( 66.3 <= 67.7 +- 0.9 <= 69.3 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 1] ( 65.8 <= 66.6 +- 0.7 <= 67.9 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 1] ( 62.3 <= 66.6 +- 1.8 <= 68.3 ) MB/s (candidates: 495)
     *  -> still faster than the non-LUT version even though there is no chunking at all!
     *
     * [13 bits, this with chunk count 2] ( 81.1 <= 84.9 +- 1.7 <= 86.4 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 2] ( 85.7 <= 87.4 +- 1.0 <= 89.2 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 2] ( 85.5 <= 86.7 +- 0.9 <= 88.2 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 2] ( 82.9 <= 84.0 +- 0.9 <= 85.6 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 2] ( 82.0 <= 84.0 +- 1.2 <= 85.2 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 2] ( 81.1 <= 83.6 +- 1.5 <= 85.2 ) MB/s (candidates: 495)
     *
     * [13 bits, this with chunk count 3] ( 86.0 <= 90.7 +- 1.8 <= 91.9 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 3] ( 94.3 <= 95.1 +- 0.6 <= 96.0 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 3] ( 92.7 <= 94.4 +- 0.9 <= 95.2 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 3] ( 84.6 <= 88.2 +- 1.5 <= 89.7 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 3] ( 87.8 <= 90.2 +- 1.3 <= 91.4 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 3] ( 87.8 <= 89.5 +- 0.9 <= 90.5 ) MB/s (candidates: 495)
     *
     * [13 bits, this with chunk count 4] ( 90.7 <= 95.3 +- 1.9 <= 96.7 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 4] ( 93.0 <= 96.1 +- 1.8 <= 97.5 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 4] ( 90.2 <= 92.8 +- 1.5 <= 94.7 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 4] ( 89.8 <= 90.8 +- 0.8 <= 91.9 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 4] ( 89.1 <= 91.0 +- 0.8 <= 91.8 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 4] ( 88.5 <= 90.5 +- 1.0 <= 91.4 ) MB/s (candidates: 495)
     *  -> This seems to be ideal. Chunk count of 5 is faster for 13 bits skip LUT, but this might
     *     just be a fluke because all others are slower.
     *
     * [13 bits, this with chunk count 5] ( 92.2 <= 95.6 +- 1.5 <= 96.9 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 5] ( 90.9 <= 91.7 +- 0.6 <= 92.8 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 5] ( 86.6 <= 88.2 +- 0.9 <= 89.3 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 5] ( 80.6 <= 85.6 +- 2.1 <= 88.5 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 5] ( 87.2 <= 88.6 +- 0.8 <= 89.3 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 5] ( 85.8 <= 86.6 +- 0.5 <= 87.4 ) MB/s (candidates: 495)
     *
     * [13 bits, this with chunk count 6] ( 88.7 <= 92.8 +- 1.7 <= 94.5 ) MB/s (candidates: 495)
     * [14 bits, this with chunk count 6] ( 88.9 <= 90.0 +- 0.7 <= 90.9 ) MB/s (candidates: 495)
     * [15 bits, this with chunk count 6] ( 83.1 <= 86.1 +- 1.4 <= 87.9 ) MB/s (candidates: 495)
     * [16 bits, this with chunk count 6] ( 81.3 <= 83.7 +- 1.8 <= 86.0 ) MB/s (candidates: 495)
     * [17 bits, this with chunk count 6] ( 80.7 <= 82.6 +- 1.2 <= 84.3 ) MB/s (candidates: 495)
     * [18 bits, this with chunk count 6] ( 76.5 <= 78.8 +- 1.4 <= 80.5 ) MB/s (candidates: 495)
     *  -> "Finally" some measurable change, probably caused by cache spills from the huge LUT.
     *     LUT is expected to be ( 2 << ( 6 * 3 ) ) * sizeof( uint16_t ) = 512 KiB
     */
    const auto precodeBits = next57Bits & nLowestBitsSet<uint64_t>( codeLengthCount * deflate::PRECODE_BITS );
    const auto virtualLeafCount = precodesToLeafCount<PRECODE_CHUNK_SIZE>( precodeBits );
#else
    /* This version without any LUTs already reaches 57 MB/s! */
    const auto virtualLeafCount = getVirtualLeafCount( next57Bits, codeLengthCount );
#endif

    /* Is allowed for a single code length of 1. It results in some additional false positives unfortunately,
     * but this does not seem to hinder performance. */
    if ( virtualLeafCount == 64 ) {
        return Error::NONE;
    }
    return virtualLeafCount == 128 ? Error::NONE : Error::INVALID_CODE_LENGTHS;

    /* The nice property of this approach is that we can return exact error codes.
     * It is not even slower, probably because the != Error::NONE check outside lets the compmiler optimize away
     * some of these calls. */
    if ( virtualLeafCount > 128 ) {
        return Error::INVALID_CODE_LENGTHS;
    }
    if ( virtualLeafCount < 128 ) {
        return Error::BLOATING_HUFFMAN_CODING;
    }
    return Error::NONE;
}
}  // namespace rapidgzip::PrecodeCheck::CountAllocatedLeaves
