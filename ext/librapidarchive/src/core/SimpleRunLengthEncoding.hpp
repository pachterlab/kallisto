#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include <core/VectorView.hpp>


namespace rapidgzip::SimpleRunLengthEncoding
{
/**
 * For ~0, then last encoded varint byte will be mostly zeros because 7 * 9  = 63, i.e., the last bit
 * will have its own byte. It would have been beautiful in its own way to move that last bit into the
 * previous bytes' high bit. This would be size-optimal, but there would be no redundancy anymore and
 * we would have to complicate the algorithm to count in order to treat the last byte differently.
 * See the example for -2: https://protobuf.dev/programming-guides/encoding/#varints
 * @verbatim
 * 11111110 11111111 11111111 11111111 11111111
 * 11111111 11111111 11111111 11111111 00000001
 * @endverbatim
 */
void
writeVarInt( std::vector<uint8_t>& target,
             uint64_t              value )
{
    do {
        target.push_back( ( ( value >> 7U ) > 0 ) ? ( value | 0b1000'0000ULL ) : ( value & 0b0111'1111ULL ) );
        value >>= 7U;
    } while ( value > 0 );
}


/**
 * Reads a varint and returns its value as well as the number of bytes read.
 * The number of read bytes is 0 on error.
 */
template<typename Container>
[[nodiscard]] constexpr std::pair<uint64_t, uint8_t>
readVarInt( const Container& source,
            size_t           offset = 0 )
{
    uint64_t value{ 0 };
    uint8_t nBytesRead{ 0 };
    for ( size_t i = offset; i < source.size(); ++i ) {
        const auto byte = source[i];
        /* The last varint byte contains the 64-th bit and nothing more! */
        if ( ( nBytesRead == 9 ) && ( byte > 1 ) ) {
            return { 0, 0 };
        }
        value += static_cast<uint64_t>( byte & 0b0111'1111ULL ) << ( 7U * nBytesRead );
        ++nBytesRead;
        if ( ( byte & 0b1000'0000U ) == 0 ) {
            return { value, nBytesRead };
        }
    }
    /* Should only happen when reaching the end of the input container while expecting further varint bytes. */
    return { 0, 0 };
}


/**
 * @return The offset and length of a repeated symbol run of specified minimum length.
 */
template<typename Container>
[[nodiscard]] constexpr std::pair<size_t, size_t>
findRun( const Container& data,
         size_t           offset = 0,
         size_t           minLength = 4 )
{
    for ( ; offset < data.size(); ++offset ) {
        size_t length = 1;
        while ( ( offset + length < data.size() ) && ( data[offset + length] == data[offset] ) ) {
            ++length;
        }
        if ( length >= minLength ) {
            return { offset, length };
        }
    }
    return { offset, 0 };
}


/**
 * @verbatim
 * (0 for literals, 1 for repeated symbols, (>1 varint could be backward reference, but not supported yet))
 * (varint number of literals/repetitions/match length)
 * (literals or the symbol to repeat)
 * (varint number for following literals) l.. literals ...
 * (varint number for repeated next symbol) symbol
 * (varint number for following literals) l.. literals ...
 * ...
 * @endverbatim
 * Number can be zero if there are no literals or something to repeat, although a literal could also
 * be stored as "number for repeated next symbol" = 1.
 * Varint is a vector of up to 10 bytes similar to UTF-8, the next bytes' bits are only read if the highest
 * bit of the previous one is set.
 * See the Protobuf varint specification. I like this format because it stores numbers space-efficiently
 * and at the same time avoids LSB / MSB issues because the number is read byte-wise anyway.
 * Subsequent bytes are for higher-order bits, i.e., they are in little-endian order, but one would not be able
 * to call a simple load / memcpy to it anyway.
 * I don't understand the LZ4 format, which ADDs not bit-shifts subsequent bytes. There is no benefit to that!?
 * Heck, if it is easier on the compressor, for streaming, it could also always write out a fixed 10-bit value
 * and simply set the higher bits to 0! This way, one could stream out gigabytes of literals and later update the
 * number.
 * @see https://protobuf.dev/programming-guides/encoding/#varints
 *
 * @note If one wanted, one could probably do this compression constexpr in-place in the input std::array.
 */
[[nodiscard]] std::vector<uint8_t>
simpleRunLengthEncode( const VectorView<uint8_t> data )
{
    std::vector<uint8_t> encoded;

    if ( data.empty() ) {
        return encoded;
    }

    /* The base case is us trying to store everything as literals.
     * While doing so, we look for repeated counts that would make us switch to a compressed write.
     * The first of the symbols to repeat would still be encoded in the literals, followed by two varints,
     * for the back-reference (offset=1, length). This means that any repeat count <= 3 would not compress!
     * If we have to start a new literal section only for a single symbol for the repeated symbols, then
     * this would require a run length >= 6 to be compressing! */
    size_t i{ 0 };
    while ( i < data.size() ) {
        /* Find literals end point by searching for repeated symbols. */
        const auto [runOffset, runLength] = findRun( data, i, 6 );

        /* Add literals until runOffset + 1 (first symbol to repeat must be included as literal). */
        writeVarInt( encoded, 0 /* literal operation */ );
        const auto literalCount = std::min( data.size() - i, runOffset + 1 - i );
        writeVarInt( encoded, literalCount );
        for ( size_t j = 0; j < literalCount; ++j ) {
            encoded.push_back( data[i + j] );
        }
        i += literalCount;

        if ( i >= data.size() ) {
            break;
        }
        /* Should not really happen. The above case should be true if the run length were zero. */
        if ( runLength <= 1 ) {
            continue;
        }

        writeVarInt( encoded, 1 /* repeat last symbol */ );
        writeVarInt( encoded, runLength - 1 );
        i += runLength - 1;
    }

    return encoded;
}


/**
 * Constexpr-compatible decode function that works with std::array or any container reference.
 * Returns the number of bytes decoded, or 0 on error.
 * The output buffer must be large enough to hold the decoded data.
 * @note Unfortunately cannot give the container as argument in the constexpr version:
 *       error: modification of ‘decoded’ (via output) from outside current evaluation is not a constant expression
 */
template<typename OutputContainer,
         typename InputContainer>
[[nodiscard]] constexpr OutputContainer
simpleRunLengthDecode( const InputContainer& data,
                       size_t                decompressedSize )
{
    OutputContainer output{};

    if constexpr ( std::is_same_v<OutputContainer, std::vector<uint8_t> > ) {
        output.resize( decompressedSize, 0 );
    } else {
        if ( ( decompressedSize > 0 ) && ( output.size() != decompressedSize ) ) {
            throw std::logic_error( "Requested decompressed size does not match container!" );
        }
    }

    size_t decodedSize = 0;
    size_t i = 0;

    while ( i < data.size() ) {
        auto [backwardReference, nBytesRead] = readVarInt( data, i );
        if ( nBytesRead == 0 ) {
            throw std::domain_error( "Partial varint read for operation type!" );
        }
        i += nBytesRead;

        /* Repeat the last decoded symbol. */
        if ( backwardReference > decodedSize ) {
            throw std::domain_error( "Backreference points past the file start!" );
        }

        /* Read varint length (literals count or match length). */
        auto [length, nBytesRead2] = readVarInt( data, i );
        if ( nBytesRead2 == 0 ) {
            throw std::domain_error( "Partial varint read for literal count/match length!" );
        }
        i += nBytesRead2;

        switch ( backwardReference )
        {
        case 0:
        {
            if ( i + length > data.size() ) {
                throw std::domain_error( "Literal count points past the end!" );
            }

            /* Copy literals. */
            for ( size_t j = 0; ( j < length ) && ( decodedSize + j < output.size() ); ++j ) {
                output[decodedSize + j] = data[i + j];
            }
            i += length;
            decodedSize += length;

            break;
        }
        case 1:
        {
            /* When output is already full, then the actual symbol does not matter. We still keep decoding
             * to count the decoded size so that it could be used to constexpr-compatibly allocate std::array! */
            const auto symbol = decodedSize - backwardReference < output.size()
                                ? output[decodedSize - backwardReference]
                                : uint8_t( 0 );
            /* Add repeated symbols to output. Skipping this for zeros is a compile time optimization!
             * It reduces instructions, making use of the zero-initialized containers! */
            if ( symbol != 0 ) {
                for ( size_t j = 0; ( j < length ) && ( decodedSize + j < output.size() ); ++j ) {
                    output[decodedSize + j] = symbol;
                }
            }
            decodedSize += length;

            break;
        }
        default:
            throw std::domain_error( "Unsupported backward reference!" );
        }
    }

    if ( decodedSize != output.size() ) {
        throw std::logic_error( "Decompressed size (" + std::to_string( decodedSize )
                                + ") does not match container (" + std::to_string( output.size() ) + ")!" );
    }

    return output;
}
}  // namespace rapidgzip::SimpleRunLengthEncoding
