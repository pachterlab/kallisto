#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include "HuffmanCodingBase.hpp"


namespace rapidgzip
{
/**
 * This is an iterative improvement over HuffmanCodingLinearSearch.
 * - During initialization, it stores all symbols (for each code) sorted by length in an array and also stores
 *   offsets to jump to subarrays of symbols with a given code length. The subarray size is given by the next offset.
 *   This avoids going over all elements all the time and also already implements usage of maximum-sized and
 *   manually managed memory chunks by using std::array to avoid heap allocations.
 * - During decoding, it reads the bits one by one and for each intermediary, checks whether
 *   there is a matching code, which is calculated ad-hoc from m_minimumCodeValuesPerLevel.
 *   If the code matches, it gets the symbol for current code length from the corresponding subarray.
 * Note that reading the bits one by one is also necessary to reverse the codes.
 */
template<typename HuffmanCode,
         uint8_t  MAX_CODE_LENGTH,
         typename Symbol,
         size_t   MAX_SYMBOL_COUNT,
         bool     CHECK_OPTIMALITY = true>
class HuffmanCodingSymbolsPerLength :
    public HuffmanCodingBase<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT, CHECK_OPTIMALITY>
{
public:
    using BaseType = HuffmanCodingBase<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT, CHECK_OPTIMALITY>;
    using BitCount = typename BaseType::BitCount;
    using CodeLengthFrequencies = typename BaseType::CodeLengthFrequencies;

protected:
    /**
     * @note Reusing this struct by calling this method multiple times is allowed. All members will be reinitialized.
     */
    constexpr void
    initializeSymbolsPerLength( const VectorView<BitCount>&  codeLengths,
                                const CodeLengthFrequencies& bitLengthFrequencies )
    {
        /* Calculate cumulative frequency sums to be used as offsets for each code-length
         * for the code-length-sorted alphabet vector. */
        size_t sum = 0;
        for ( uint8_t bitLength = this->m_minCodeLength; bitLength <= this->m_maxCodeLength; ++bitLength ) {
            m_offsets[bitLength - this->m_minCodeLength] = static_cast<uint16_t>( sum );
            sum += bitLengthFrequencies[bitLength];
        }
        m_offsets[this->m_maxCodeLength - this->m_minCodeLength + 1] = static_cast<uint16_t>( sum );

        /* The codeLengths.size() check above should implicitly check this already. */
        assert( sum <= m_symbolsPerLength.size() && "Specified max symbol range exceeded!" );
        assert( sum <= std::numeric_limits<uint16_t>::max() && "Symbol count limited to 16-bit!" );

        /* Fill the code-length-sorted alphabet vector. */
        auto sizes = m_offsets;
        for ( size_t symbol = 0; symbol < codeLengths.size(); ++symbol ) {
            if ( codeLengths[symbol] != 0 ) {
                const auto k = codeLengths[symbol] - this->m_minCodeLength;
                m_symbolsPerLength[sizes[k]++] = static_cast<Symbol>( symbol );
            }
        }
    }

public:
    [[nodiscard]] constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        if ( const auto errorCode = BaseType::initializeMinMaxCodeLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        CodeLengthFrequencies bitLengthFrequencies = {};
        for ( const auto value : codeLengths ) {
            ++bitLengthFrequencies[value];
        }

        if ( const auto errorCode = BaseType::checkCodeLengthFrequencies( bitLengthFrequencies, codeLengths.size() );
             errorCode != Error::NONE ) {
            return errorCode;
        }

        BaseType::initializeMinimumCodeValues( bitLengthFrequencies );  // resets bitLengthFrequencies[0] to 0!

        initializeSymbolsPerLength( codeLengths, bitLengthFrequencies );

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
        for ( BitCount i = 0; i < this->m_minCodeLength; ++i ) {
            code = static_cast<HuffmanCode>( HuffmanCode( code << 1U ) | bitReader.template read<1>() );
        }

        for ( BitCount k = 0; k <= this->m_maxCodeLength - this->m_minCodeLength; ++k ) {
            const auto minCode = this->m_minimumCodeValuesPerLevel[k];
            if ( minCode <= code ) {
                const auto subIndex = m_offsets[k] + static_cast<size_t>( code - minCode );
                if ( subIndex < m_offsets[k + 1] ) {
                    return m_symbolsPerLength[subIndex];
                }
            }

            code <<= 1;
            code |= bitReader.template read<1>();
        }

        return std::nullopt;
    }

protected:
    /**
     * Contains the alphabet, first sorted by code length, then by given alphabet order. E.g., it could look like this:
     * @verbatim
     * +-------+-----+---+
     * | B D E | A F | C |
     * +-------+-----+---+
     *   CL=3   CL=4  CL=5
     * @endverbatim
     * The starting index for a given code length (CL) can be queried with m_offsets.
     */
    alignas( 64 ) std::array<Symbol, MAX_SYMBOL_COUNT> m_symbolsPerLength{};
    /* +1 because it stores the size at the end as well as 0 in the first element, which might be redundant but fast. */
    alignas( 64 ) std::array<uint16_t, MAX_CODE_LENGTH + 1> m_offsets{};
    static_assert( MAX_SYMBOL_COUNT + MAX_CODE_LENGTH <= std::numeric_limits<uint16_t>::max(),
                   "Offset type must be able to point at all symbols!" );
};
}  // namespace rapidgzip
