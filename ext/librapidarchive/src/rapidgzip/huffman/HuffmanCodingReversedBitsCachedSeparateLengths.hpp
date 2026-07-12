#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>

#include <huffman/HuffmanCodingSymbolsPerLength.hpp>
#include <rapidgzip/gzip/definitions.hpp>


namespace rapidgzip
{
/**
 * Same as HuffmanCodingReversedBitsCached but the lengths are stored separately, requiring an additional lookup
 * into a much smaller lookup table.
 */
template<typename HuffmanCode,
         uint8_t  MAX_CODE_LENGTH,
         typename Symbol,
         size_t   MAX_SYMBOL_COUNT>
class HuffmanCodingReversedBitsCachedSeparateLengths :
    public HuffmanCodingSymbolsPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>
{
public:
    using BaseType = HuffmanCodingSymbolsPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>;
    using BitCount = typename BaseType::BitCount;
    using CodeLengthFrequencies = typename BaseType::CodeLengthFrequencies;

public:
    [[nodiscard]] constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        if ( const auto errorCode = BaseType::initializeFromLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        if ( codeLengths.size() > MAX_SYMBOL_COUNT ) {
            throw std::invalid_argument( "Number of code lengths exceeds the maximum symbol count!" );
        }

        if ( m_needsToBeZeroed ) {
            // Works constexpr
            for ( size_t symbol = 0; symbol < ( 1ULL << this->m_maxCodeLength ); ++symbol ) {
                m_codeCache[symbol] = 0;
            }
        }

        auto codeValues = this->m_minimumCodeValuesPerLevel;
        for ( size_t symbol = 0; static_cast<size_t>( symbol ) < codeLengths.size(); ++symbol ) {
            const auto length = codeLengths[symbol];
            m_codeLengths[symbol + 1] = length;
            if ( length == 0 ) {
                continue;
            }

            const auto k = length - this->m_minCodeLength;
            const auto code = codeValues[k]++;
            const auto reversedCode = reverseBits( code, length );

            const auto fillerBitCount = this->m_maxCodeLength - length;
            const auto maximumPaddedCode = static_cast<HuffmanCode>(
                reversedCode | ( nLowestBitsSet<HuffmanCode>( fillerBitCount ) << length ) );
            assert( maximumPaddedCode < m_codeCache.size() );
            const auto increment = static_cast<HuffmanCode>( HuffmanCode( 1 ) << length );
            for ( auto paddedCode = reversedCode; paddedCode <= maximumPaddedCode; paddedCode += increment ) {
                m_codeCache[paddedCode] = symbol + 1;
            }
        }

        m_needsToBeZeroed = true;

        return Error::NONE;
    }

    [[nodiscard]] forceinline std::optional<Symbol>
    decode( gzip::BitReader& bitReader ) const
    {
        try {
            const auto value = bitReader.peek( this->m_maxCodeLength );

            assert( value < m_codeCache.size() );
            auto symbol = m_codeCache[(int)value];
            if ( symbol == 0 ) {
                /* This might happen for non-optimal Huffman trees out of which all except the case of a single
                 * symbol with bit length 1 are forbidden! */
                return std::nullopt;
            }

            const auto length = m_codeLengths[symbol];

            bitReader.seekAfterPeek( length );
            return symbol - 1;
        } catch ( const gzip::BitReader::EndOfFileReached& ) {
            /* Should only happen at the end of the file and probably not even there
             * because the gzip footer should be longer than the peek length. */
            return BaseType::decode( bitReader );
        }
    }

private:
    alignas( 8 ) std::array<uint8_t, MAX_SYMBOL_COUNT> m_codeLengths{};
    alignas( 8 ) std::array<Symbol, ( 1UL << MAX_CODE_LENGTH )> m_codeCache{};
    bool m_needsToBeZeroed{ false };
};
}  // namespace rapidgzip
