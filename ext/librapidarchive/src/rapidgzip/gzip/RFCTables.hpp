#pragma once

#include <array>
#include <cstdint>

#include "definitions.hpp"


namespace rapidgzip::deflate
{
/* Distance Code Table */

/**
 * Currently, only used for tests.
 */
[[nodiscard]] constexpr uint16_t
calculateDistance( uint16_t distance,
                   uint8_t  extraBitsCount,
                   uint16_t extraBits ) noexcept
{
    assert( distance >= 4 );
    return 1U + ( 1U << ( extraBitsCount + 1U ) ) + ( ( distance % 2U ) << extraBitsCount ) + extraBits;
};


[[nodiscard]] constexpr uint16_t
calculateDistanceExtraBits( uint16_t distance ) noexcept
{
    return distance <= 3U ? 0U : ( distance - 2U ) / 2U;
}


/**
 * This only makes sense to use for LUT creation because, else, calculating the extra bits count
 * would be work done twice.
 * @return A kind of intermediary distance. In order to get the real distance,
 *         extraBits need to be added to the returned value.
 */
[[nodiscard]] constexpr uint16_t
calculateDistance( uint16_t distance ) noexcept
{
    assert( distance >= 4 );
    const auto extraBitsCount = calculateDistanceExtraBits( distance );
    return 1U + ( 1U << ( extraBitsCount + 1U ) ) + ( ( distance % 2U ) << extraBitsCount );
};


using DistanceLUT = std::array<uint16_t, 30>;

[[nodiscard]] constexpr DistanceLUT
createDistanceLUT() noexcept
{
    DistanceLUT result{};
    for ( uint16_t i = 0; i < 4; ++i ) {
        result[i] = i + 1;
    }
    for ( uint16_t i = 4; i < static_cast<uint16_t>( result.size() ); ++i ) {
        result[i] = calculateDistance( i );
    }
    return result;
}


alignas( 8 ) static constexpr DistanceLUT
distanceLUT = createDistanceLUT();


/* Length Code Table */


[[nodiscard]] constexpr uint16_t
calculateLength( uint16_t code ) noexcept
{
    assert( code < 285 - 261 );
    const auto extraBits = code / 4U;
    return 3U + ( 1U << ( extraBits + 2U ) ) + ( ( code % 4U ) << extraBits );
};


using LengthLUT = std::array<uint16_t, 285 - 261>;

[[nodiscard]] constexpr LengthLUT
createLengthLUT() noexcept
{
    LengthLUT result{};
    for ( uint16_t i = 0; i < static_cast<uint16_t>( result.size() ); ++i ) {
        result[i] = calculateLength( i );
    }
    return result;
}


alignas( 8 ) static constexpr LengthLUT
lengthLUT = createLengthLUT();


[[nodiscard]] inline uint16_t
getLength( uint16_t         code,
           gzip::BitReader& bitReader )
{
    if ( code <= 264 ) {
        return code - 257U + 3U;
    }
    if ( code < 285 ) {
        code -= 261;
        const auto extraBits = code / 4;
        return calculateLength( code ) + bitReader.read( extraBits );
    }
    if ( code == 285 ) {
        return 258;
    }

    throw std::invalid_argument( "Invalid Code!" );
}


[[nodiscard]] inline uint8_t
getLengthMinus3( uint16_t         code,
                 gzip::BitReader& bitReader )
{
    /* The largest length is 258 and the smallest is 3, so length-3 fits into uint8_t. */
    return static_cast<uint8_t>( getLength( code, bitReader ) - 3U );
}


template<typename DistanceHuffmanCoding>
inline std::pair<uint16_t, Error>
getDistance( CompressionType              compressionType,
             const DistanceHuffmanCoding& distanceHC,
             gzip::BitReader&             bitReader )
{
    uint16_t distance = 0;
    if ( compressionType == CompressionType::FIXED_HUFFMAN ) {
        distance = reverseBits( static_cast<uint8_t>( bitReader.read<5>() ) ) >> 3U;
        if ( UNLIKELY( distance >= MAX_DISTANCE_SYMBOL_COUNT ) ) [[unlikely]] {
            return { 0, Error::EXCEEDED_DISTANCE_RANGE };
        }
    } else {
        const auto decodedDistance = distanceHC.decode( bitReader );
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
}  // namespace rapidgzip::deflate
