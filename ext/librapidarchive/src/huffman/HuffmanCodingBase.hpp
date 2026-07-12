#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include <core/Error.hpp>
#include <core/VectorView.hpp>


namespace rapidgzip
{
template<typename T_HuffmanCode,
         uint8_t  T_MAX_CODE_LENGTH,
         typename T_Symbol,
         size_t   T_MAX_SYMBOL_COUNT,
         bool     CHECK_OPTIMALITY = true>
class HuffmanCodingBase
{
public:
    using HuffmanCode = T_HuffmanCode;
    using Symbol = T_Symbol;

    static constexpr auto MAX_CODE_LENGTH = T_MAX_CODE_LENGTH;
    static constexpr auto MAX_SYMBOL_COUNT = T_MAX_SYMBOL_COUNT;

    static_assert( MAX_CODE_LENGTH <= std::numeric_limits<HuffmanCode>::digits,
                   "The huffman code type must fit the max code length!" );
    static_assert( MAX_SYMBOL_COUNT <= std::numeric_limits<Symbol>::max(),
                   "The symbol type must fit the highest possible symbol!" );

    using BitCount = uint8_t;
    using Frequency = HuffmanCode;
    using CodeLengthFrequencies = std::array<Frequency, MAX_CODE_LENGTH + 1>;

    [[nodiscard]] bool
    isValid() const
    {
        return m_minCodeLength <= m_maxCodeLength;
    }

protected:
    constexpr Error
    initializeMinMaxCodeLengths( const VectorView<BitCount>& codeLengths )
    {
        static_assert( std::is_unsigned_v<HuffmanCode>, "Huffman code type must be unsigned" );

        if ( UNLIKELY( codeLengths.empty() ) ) [[unlikely]] {
            return Error::EMPTY_ALPHABET;
        }

        if ( UNLIKELY( codeLengths.size() > MAX_SYMBOL_COUNT ) ) [[unlikely]] {
            throw std::invalid_argument( "The range of the symbol type cannot represent the implied alphabet!" );
        }

        /* A maximum code length of 0 is valid! It happens when encoding this with pigz:
         * python3 -c 'import sys; sys.stdout.buffer.write(bytes(range(256)))' | pigz > 0CL.pigz */
        m_maxCodeLength = getMax( codeLengths );

        m_minCodeLength = getMinPositive( codeLengths );
        if ( UNLIKELY( m_maxCodeLength > MAX_CODE_LENGTH ) ) [[unlikely]] {
            throw std::invalid_argument( "The range of the code type cannot represent the given code lengths!" );
        }

        return Error::NONE;
    }

    constexpr Error
    checkCodeLengthFrequencies( const CodeLengthFrequencies& bitLengthFrequencies,
                                size_t                       codeLengthsSize )
    {
        /* E.g.:
         * fixed tree has 255 - 144 + 1 = 112 symbols of length 9
         * fixed tree has (143 - 0 + 1) + (287 - 280 + 1) = 152 symbols of length 8
         *   -> can store 256 values but because of the previous one only 256 - 112/2 = 200
         *      -> Note that the rest, i.e., 200-152 = 48 = 24*2 values cannot be used
         *         because the tree ends in the lower level for code length 7!
         * fixed tree has 279 - 256 + 1 = 24 symbols of length 7
         *   -> can store 128 - 152 / 2 - 112 / 4 = 24
         *      -> Exactly the amount actually being used!
         * -> The tree is full utilized on all levels! I think this makes it possible to do stronger checks!
         *    The only reasonable case for unused tree nodes is if there is only one letter.
         *    E.g., consider a single letter, which only uses half the tree with one bit.
         *    But consider 3 letters. Making all of CL 2 makes no sense, because the third could have been
         *    one bit shorter without any problems.
         *      -> this basically corresponds to the check that the number of letters with the longest code length
         *         should be even!
         * I think with this, it is always possible to exactly use up all leaf nodes for any alphabet,
         * and we could test for that! */
        const auto nonZeroCount = codeLengthsSize - bitLengthFrequencies[0];
        HuffmanCode unusedSymbolCount = HuffmanCode( 1 ) << m_minCodeLength;
        for ( int bitLength = m_minCodeLength; bitLength <= m_maxCodeLength; ++bitLength ) {
            const auto frequency = bitLengthFrequencies[bitLength];
            if ( frequency > unusedSymbolCount ) {
                return Error::INVALID_CODE_LENGTHS;
            }
            unusedSymbolCount -= frequency;
            unusedSymbolCount *= 2;  /* Because we go down one more level for all unused tree nodes! */
        }

        if constexpr ( CHECK_OPTIMALITY ) {
            if ( ( ( nonZeroCount == 1 ) && ( unusedSymbolCount != ( 1U << m_maxCodeLength ) ) ) ||
                 ( ( nonZeroCount >  1 ) && ( unusedSymbolCount != 0 ) ) ) {
                return Error::BLOATING_HUFFMAN_CODING;
            }
        }

        return Error::NONE;
    }

    /**
     * @note Reusing this struct by calling this method multiple times is allowed. All members will be reinitialized.
     */
    constexpr void
    initializeMinimumCodeValues( CodeLengthFrequencies& bitLengthFrequencies )
    {
        /* Bit lengths of zero mean, that the particular "letter" is not used at all and should be skipped
         * for calculating the codes. However, the value must be zero for the initial condition of the next
         * algorithms to work. */
        bitLengthFrequencies[0] = 0;

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
        HuffmanCode minCode = 0;
        /* minCodeLength might be zero for empty deflate blocks as can happen when compressing an empty file! */
        for ( size_t bits = std::max<size_t>( 1U, m_minCodeLength ); bits <= m_maxCodeLength; ++bits ) {
            minCode = static_cast<HuffmanCode>( HuffmanCode( minCode + bitLengthFrequencies[bits - 1U] ) << 1U );
            m_minimumCodeValuesPerLevel[bits - m_minCodeLength] = minCode;
        }
    }

public:
    /**
     * Creates a Huffman coding for the implicit alphabet (0, ..., codeLengths.size() - 1)
     * using the respectively given code lengths for the Huffman codes.
     *
     * A Huffman tree gives a variable-length binary code for each leaf node representing a symbol.
     * Further restrictions on Huffman tree by deflate:
     *  - Maximum code length: 15 because the deflate code length alphabet only goes from 0 to 15.
     *  - Shorter codes lexicographically precede longer codes.
     *  - All codes of a given bit length have lexicographically consecutive values,
     *    in the same order as the symbols they represent
     * With this, the Huffman tree can be exactly specified as a vector, whose elements represent
     * the ordered alphabet's symbol's code lengths, e.g, for alphabet ABCD: (2,1,3,3) -> 0b10, 0b0, 0b110, 0b111
     */
    constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        if ( const auto errorCode = initializeMinMaxCodeLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        CodeLengthFrequencies bitLengthFrequencies = {};
        for ( const auto value : codeLengths ) {
            ++bitLengthFrequencies[value];
        }

        if ( const auto errorCode = checkCodeLengthFrequencies( bitLengthFrequencies, codeLengths.size() );
             errorCode != Error::NONE ) {
            return errorCode;
        }

        initializeMinimumCodeValues( bitLengthFrequencies );  // resets bitLengthFrequencies[0] to 0!

        return Error::NONE;
    }

    [[nodiscard]] constexpr BitCount
    minCodeLength() const noexcept
    {
        return m_minCodeLength;
    }

    [[nodiscard]] constexpr BitCount
    maxCodeLength() const noexcept
    {
        return m_maxCodeLength;
    }

    [[nodiscard]] constexpr auto const&
    minimumCodeValuesPerLevel() const noexcept
    {
        return m_minimumCodeValuesPerLevel;
    }

protected:
    BitCount m_minCodeLength{ std::numeric_limits<BitCount>::max() };
    BitCount m_maxCodeLength{ std::numeric_limits<BitCount>::min() };

    /** Only indexes [0, m_maxCodeLength - m_minCodeLength) contain valid data! */
    std::array<HuffmanCode, MAX_CODE_LENGTH + 1> m_minimumCodeValuesPerLevel{};
};


template<uint8_t  MAX_CODE_LENGTH,
         typename CodeLengths>
[[nodiscard]] constexpr bool
checkHuffmanCodeLengths( const CodeLengths& codeLengths )
{
    size_t virtualLeafCount{ 0 };
    for ( const auto codeLength : codeLengths ) {
        virtualLeafCount += codeLength > 0 ? 1U << ( MAX_CODE_LENGTH - codeLength ) : 0U;
    }

    if ( virtualLeafCount == ( 1U << ( MAX_CODE_LENGTH - 1U ) ) ) {
        size_t greaterThanOne{ 0 };
        for ( const auto codeLength : codeLengths ) {
            if ( codeLength > 1 ) {
                ++greaterThanOne;
            }
        }
        return greaterThanOne == 0;
    }
    return virtualLeafCount == ( 1U << MAX_CODE_LENGTH );
}
}  // namespace rapidgzip
