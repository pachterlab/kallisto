#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <core/common.hpp>
#include <core/Error.hpp>
#include <core/VectorView.hpp>


namespace rapidgzip
{
/**
 * This was the first and is the most straight-forward implementation.
 * - During initialization, it stores all codes and their lengths in two vectors.
 * - During decoding, it reads the bits one by one and for each intermediary, checks whether
 *   there is a matching code with the current length in the vectors.
 * Note that reading the bits one by one is also necessary to reverse the codes.
 */
template<typename T_HuffmanCode,
         typename T_Symbol>
class HuffmanCodingLinearSearch
{
public:
    using HuffmanCode = T_HuffmanCode;
    using Symbol = T_Symbol;
    using BitCount = uint8_t;

    static constexpr auto MAX_CODE_LENGTH = std::numeric_limits<HuffmanCode>::digits;

    [[nodiscard]] bool
    isValid() const
    {
        return !m_codes.empty() && ( m_minCodeLength > 0 );
    }

    [[nodiscard]] std::vector<HuffmanCode>
    countFrequencies( const std::vector<uint8_t>& values )
    {
        std::vector<HuffmanCode> frequencies( std::numeric_limits<uint8_t>::max(), 0 );
        for ( const auto value : values ) {
            ++frequencies[value];
        }
        return frequencies;
    }

public:
    /**
     * Creates a Huffman coding for the alphabet (0, ..., codeLengths.size() - 1) using the respective
     * bit lengths for the Huffman codes.
     * @note Reusing this struct by calling this method multiple times is allowed. All members will be reinitialized.
     */
    Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        m_codeLengths.assign( codeLengths.begin(), codeLengths.end() );
        m_minCodeLength = getMinPositive( m_codeLengths );
        m_maxCodeLength = getMax( m_codeLengths );

        static_assert( std::is_unsigned_v<HuffmanCode>, "Huffman code type must be unsigned" );

        /* Huffman tree give variable-length binary codes for each leaf node representing a symbol
         * Maximum code length: ???
         * Further restrictions on Huffman tree by Deflate:
         *  - Shorter codes lexicographically precede longer codes.
         *  - All codes of a given bit length have lexicographically consecutive values,
         *    in the same order as the symbols they represent
         * With this, the Huffman tree can be exactly specified as a vector, whose elements represent
         * the ordered alphabet's symbol's code lengths, e.g, for alphabet ABCD: (2,1,3,3) -> 0b10, 0b0, 0b110, 0b111
         */

        const auto& bitLengths = m_codeLengths;
        auto bitLengthFrequencies = countFrequencies( bitLengths );
        if ( codeLengths.size() > std::numeric_limits<std::decay_t<decltype( bitLengthFrequencies[0] )> >::max() ) {
            throw std::logic_error( "The frequency count type must fit the count even if all code lengths are equal!" );
        }

        /* It's important for the check to know the highest non-zero value. */
        const auto lastNonZero = std::find_if( bitLengthFrequencies.rbegin(), bitLengthFrequencies.rend(),
                                               [] ( const auto value ) { return value != 0; } );
        bitLengthFrequencies.resize( bitLengthFrequencies.rend() - lastNonZero );

        if ( bitLengthFrequencies.empty() ) {
            return Error::EMPTY_INPUT;
        }

        /* Bit lengths of zero mean, that the particular "letter" is not used at all and should be skipped
         * for calculating the codes. However, the value must be zero for the initial condition of the next
         * algorithms to work. */
        bitLengthFrequencies[0] = 0;

        /* Check whether these values make sense. A bit length can only represent finite amount of codes:
         * Bit length 1 -> 2 if no longer bit length, 1 if there are longer bit lengths, 0 if any of the next levels
         * contain more than 2^(level-1) leafs/symbols! It is easier to check this in reverse! If the last level is
         * completely filled with symbols, then the level before it cannot have any leaf nodes. Or more generically,
         * half the amount of possible values of the previous level are filled up for being parents to the leaf nodes of
         * the next level. */
        auto maxSymbolsPossible = uint8_t( 1 ) << bitLengthFrequencies.size();
        for ( int bitLength = bitLengthFrequencies.size() - 1; bitLength > 0; --bitLength ) {
            const auto frequency = bitLengthFrequencies.begin() + bitLength;
            if ( *frequency > maxSymbolsPossible ) {
                return Error::EXCEEDED_CL_LIMIT;
            }

            /* E.g.:
             * fixed tree has 255 - 144 + 1 = 112 symbols of length 9
             * fixed tree has 143 - 0 + 1 + 287 - 280 + 1 = 152 symbols of length 8
             *   -> can store 256 values but because of the previous one only 256 - 112/2 = 200
             *      -> not even utilized fully ...?
             * fixed tree has 279 - 256 + 1 = 24 symbols of length 7
             *   -> can store 128 - 152 / 2 = 52
             *      -> also not fully utilized ... Why use such a suboptimal fixed tree?! */
            maxSymbolsPossible = ( uint8_t( 1 ) << ( bitLength - 1 ) )
                                 - ceilDiv( *frequency, 2 );  // ( *frequency >> 2U ) + ( frequency & 1 );
        }

        /**
         * Find the minimum code value for each bitLength to start assigning values from.
         * This is basically the above test in reverse plus assuming that lower tree levels
         * are filled from the left, i.e., the highest code values on a level are all non-leaf nodes.
         * @verbatim
         *          O  -> 0         level / bit length 1
         *         / \
         *  01 <- A   O -> 10           1
         *           / \
         *   100 <- B   C -> 101        2
         * @endverbatim
         * This would be encoded as (ABC) -> (1,3,3) => frequencies: (1,2,3) -> (1,0,2).
         * The smallest code length value on level 2 (or l in general), is the smallest
         * non-leaf node on level l-1 with a 0 appended in binary, i.e., multiplied by 2.
         * The smallest non-leaf node on a level is given by the smallest code value on
         * the previous level plus the number of leaf nodes on that level, i.e., symbols
         * with the previous code length.
         */
        std::vector<HuffmanCode> minimumCodeValuesPerLevel( bitLengthFrequencies.size() + 1 );
        minimumCodeValuesPerLevel[0] = 0;
        for ( size_t bits = 1; bits <= bitLengthFrequencies.size(); ++bits ) {
            minimumCodeValuesPerLevel[bits] = ( minimumCodeValuesPerLevel[bits - 1]
                                                + bitLengthFrequencies[bits - 1] ) << 1U;
        }

        /* Now, begin assigning the alphabet consecutively to the codes starting from the minimumCodes.
         * @todo It might be possible to merge this with the minimumCodeValuesPerLevel creation loop to avoid
         * that intermediary vector and hopefully shave off time.
         * Note that this for loop destructibly works on minimumCodeValuesPerLevel anyway.
         * @todo The minimum codes alone might even be enough to decode Huffman codes! */
        std::vector<HuffmanCode> huffmanCodes( bitLengths.size() );
        for ( size_t i = 0; i < bitLengths.size(); ++i ) {
            const auto bitLength = bitLengths[i];
            if ( bitLength != 0 ) {
                huffmanCodes[i] = minimumCodeValuesPerLevel[bitLength]++;
            }
        }

        m_codes = std::move( huffmanCodes );

        #if 1
        if ( m_codeLengths.size() > std::numeric_limits<Symbol>::max() ) {
            throw std::invalid_argument( "The range of the symbol type cannot represent the implied alphabet!" );
        }

        for ( const auto codeLength : m_codeLengths ) {
            if ( codeLength > std::numeric_limits<HuffmanCode>::digits ) {
                std::stringstream message;
                message << "The range of the code type cannot represent the given code lengths!\n"
                        << "Got length " << static_cast<int>( codeLength ) << " but code type width is "
                        << std::numeric_limits<HuffmanCode>::digits << "\n";
                throw std::invalid_argument( std::move( message ).str() );
            }
        }
        #endif

        return Error::NONE;
    }

    template<typename BitReader>
    [[nodiscard]] forceinline constexpr std::optional<Symbol>
    decode( BitReader& bitReader ) const
    {
        HuffmanCode code = 0;

        /** Read the first n bytes. Note that we can't call the bitReader with argument > 1 because the bit order
         * would be inversed. @todo Reverse the Huffman codes and prepend bits instead of appending, so that this
         * first step can be conflated and still have the correct order for comparison! */
        for ( BitCount i = 0; i < m_minCodeLength; ++i ) {
            code = ( code << 1U ) | ( bitReader.template read<1>() );
        }

        for ( BitCount bitLength = m_minCodeLength; bitLength <= m_maxCodeLength; ++bitLength ) {
            /* Look up whether it is a huffman code.
             * @todo faster lookup than simple linear search, e.g., binning by length and binary search between.
             *       -> then again binary search is not even necessary because the values in each bin are given
             *          sequentially. Might be possible to simply add the lookup value to the minimumValueCode or so! */
            for ( size_t j = 0; j < m_codeLengths.size(); ++j ) {
                if ( ( m_codeLengths[j] == bitLength ) && ( m_codes[j] == code ) ) {
                    return static_cast<Symbol>( j );
                }
            }

            code <<= 1;
            code |= bitReader.template read<1>();
        }

        return std::nullopt;
    }

private:
    std::vector<BitCount> m_codeLengths;
    std::vector<HuffmanCode> m_codes;
    BitCount m_minCodeLength{ 0 };
    BitCount m_maxCodeLength{ 0 };
};
}  // namespace rapidgzip
