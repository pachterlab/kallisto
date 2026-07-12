#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include <huffman/HuffmanCodingSymbolsPerLength.hpp>
#include <rapidgzip/gzip/definitions.hpp>
#include <rapidgzip/gzip/RFCTables.hpp>


namespace rapidgzip::deflate
{
/**
 * This version uses a lookup table (LUT) to avoid repetitive loops over BitReader::read<1>() to speed things up a lot.
 * This started as a copy of @ref HuffmanCodingReversedBitsCached, however:
 * - It limits the LUT to a fixed size instead of caching everything up to MAX_CODE_LENGTH because that
 *   would be too large for bzip2, for which MAX_CODE_LENGTH is 20 instead of 16 for gzip.
 */
template<uint8_t LUT_BITS_COUNT>
class HuffmanCodingShortBitsMultiCached :
    public HuffmanCodingSymbolsPerLength<uint16_t, MAX_CODE_LENGTH, uint16_t, MAX_LITERAL_HUFFMAN_CODE_COUNT,
                                         /* CHECK_OPTIMALITY */ true>
{
public:
    using Symbol = uint16_t;
    using HuffmanCode = uint16_t;
    using BaseType = HuffmanCodingSymbolsPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_LITERAL_HUFFMAN_CODE_COUNT,
                                                   /* CHECK_OPTIMALITY */ true>;
    using BitCount = typename BaseType::BitCount;
    using CodeLengthFrequencies = typename BaseType::CodeLengthFrequencies;

    struct CacheEntry
    {
        bool needToReadDistanceBits : 1;
        uint8_t bitsToSkip : 6;  // ceil(log2 2*MAX_CODE_LENGTH(20)) = 6 bits would suffice
        uint8_t symbolCount : 2;
        /**
         * Contains one 8-9-bit symbol in the lower bits and one 8-9 bit symbol in the higher bits only if
         * the first symbol is an 8-bit literal.
         */
        uint32_t symbols : 18;
    };
    static_assert( sizeof( CacheEntry ) == 4 );

    using Symbols = std::pair</* packed symbols */ uint32_t, /* symbol count */ uint32_t>;

    static constexpr size_t DISTANCE_OFFSET = 254U;

public:
    [[nodiscard]] constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        if ( const auto errorCode = BaseType::initializeFromLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        m_lutBitsCount = std::min( LUT_BITS_COUNT, this->m_maxCodeLength );
        m_bitsToReadAtOnce = std::max( LUT_BITS_COUNT, this->m_minCodeLength );

        /* Initialize the cache. */
        if ( m_needsToBeZeroed ) {
            // Works constexpr
            for ( size_t symbol = 0; symbol < m_codeCache.size(); ++symbol ) {
                m_codeCache[symbol].bitsToSkip = 0;
            }
        }

        /* Compute Huffman values. */

        struct HuffmanEntry
        {
            HuffmanCode reversedCode{ 0 };
            Symbol symbol{ 0 };
            uint8_t length{ 0 };
        };

        std::array<HuffmanEntry, MAX_LITERAL_HUFFMAN_CODE_COUNT> huffmanTable{};
        size_t huffmanTableSize{ 0 };

        if ( codeLengths.size() > huffmanTable.size() ) {
            return Error::INVALID_CODE_LENGTHS;
        }

        auto codeValues = this->m_minimumCodeValuesPerLevel;
        for ( size_t symbol = 0; symbol < codeLengths.size(); ++symbol ) {
            const auto length = codeLengths[symbol];
            /* Note that the preemptive continue for large lengths is not a bug even if we do not increment
             * codeValues because all symbols of the same length will either be filtered or not, so that
             * the missing increment does not matter because that code is either always used or never. */
            if ( ( length == 0 ) || ( length > m_lutBitsCount ) ) {
                continue;
            }

            auto& entry = huffmanTable[huffmanTableSize++];
            entry.length = length;
            entry.symbol = static_cast<Symbol>( symbol );
            entry.reversedCode = reverseBits( codeValues[length - this->m_minCodeLength]++, length );
        }

        /* Fill up cache. */

        for ( size_t i = 0; i < huffmanTableSize; ++i ) {
            const auto& huffmanEntry = huffmanTable[i];

            CacheEntry cacheEntry{};
            cacheEntry.bitsToSkip = huffmanEntry.length;
            cacheEntry.symbols = huffmanEntry.symbol;
            cacheEntry.symbolCount = 1;
            cacheEntry.needToReadDistanceBits = huffmanEntry.symbol > END_OF_BLOCK_SYMBOL;

            if ( cacheEntry.needToReadDistanceBits ) {
                insertLengthSymbolIntoCache( huffmanEntry.reversedCode, cacheEntry );
                continue;
            }

            insertIntoCache( huffmanEntry.reversedCode, cacheEntry );
        }

        m_needsToBeZeroed = true;

        return Error::NONE;
    }

    template<typename BitReader>
    [[nodiscard]] forceinline Symbols
    decode( BitReader& bitReader ) const
    {
        try {
            const auto& cacheEntry = m_codeCache[bitReader.peek( m_lutBitsCount )];
            if ( cacheEntry.bitsToSkip == 0 ) {
                return decodeLong( bitReader );
            }
            bitReader.seekAfterPeek( cacheEntry.bitsToSkip );
            return {
                cacheEntry.needToReadDistanceBits
                ? readLength( cacheEntry.symbols, bitReader )
                : cacheEntry.symbols,
                cacheEntry.symbolCount
            };
        } catch ( const typename BitReader::EndOfFileReached& ) {
            /* Should only happen at the end of the file and probably not even there
             * because the bzip2 footer (EOS block) should be longer than the peek length. */
            const auto result = BaseType::decode( bitReader );
            if ( result ) {
                return { readLength( *result, bitReader ), 1 };
            }
            return { 0, 0 };
        }
    }

private:
    template<typename BitReader>
    [[nodiscard]] forceinline constexpr Symbols
    decodeLong( BitReader& bitReader ) const
    {
        HuffmanCode code = 0;

        for ( BitCount i = 0; i < m_bitsToReadAtOnce; ++i ) {
            code = ( code << 1U ) | ( bitReader.template read<1>() );
        }

        for ( BitCount k = m_bitsToReadAtOnce - this->m_minCodeLength;
              k <= this->m_maxCodeLength - this->m_minCodeLength; ++k )
        {
            const auto minCode = this->m_minimumCodeValuesPerLevel[k];
            if ( minCode <= code ) {
                const auto subIndex = this->m_offsets[k] + static_cast<size_t>( code - minCode );
                if ( subIndex < this->m_offsets[k + 1] ) {
                    return { readLength( this->m_symbolsPerLength[subIndex], bitReader ), 1 };
                }
            }

            code <<= 1;
            code |= bitReader.template read<1>();
        }

        return { 0, 0 };
    }

    template<typename BitReader>
    [[nodiscard]] forceinline constexpr Symbol
    readLength( const Symbol symbol,
                BitReader&   bitReader ) const
    {
        if ( symbol <= 256 ) {
            return symbol;
        }
        /* Basically the same as ( 1 << 8U ) | getLengthMinus3 */
        return getLength( symbol, bitReader ) + DISTANCE_OFFSET;
    }

    forceinline void
    insertIntoCache( HuffmanCode reversedCode,
                     CacheEntry  cacheEntry )
    {
        const auto length = cacheEntry.bitsToSkip;
        if ( length > m_lutBitsCount ) {
            return;
        }
        const auto fillerBitCount = m_lutBitsCount - length;

        const auto maximumPaddedCode = static_cast<HuffmanCode>(
            reversedCode | ( nLowestBitsSet<HuffmanCode>( fillerBitCount ) << length ) );
        assert( maximumPaddedCode < m_codeCache.size() );
        const auto increment = static_cast<HuffmanCode>( HuffmanCode( 1 ) << length );
        for ( auto paddedCode = reversedCode; paddedCode <= maximumPaddedCode; paddedCode += increment ) {
            m_codeCache[paddedCode] = cacheEntry;
        }
    }

    forceinline void
    insertLengthSymbolIntoCache( const HuffmanCode reversedCode,
                                 const CacheEntry& inputCacheEntry )
    {
        if ( !inputCacheEntry.needToReadDistanceBits ) {
            insertIntoCache( reversedCode, inputCacheEntry );
            return;
        }

        const auto previousBitCount = ( inputCacheEntry.symbolCount - 1U ) * BYTE_SIZE;
        const auto symbol = static_cast<uint32_t>( inputCacheEntry.symbols ) >> previousBitCount;
        const auto codeLength = inputCacheEntry.bitsToSkip;
        const auto previousSymbols = symbol & nLowestBitsSet<uint32_t>( previousBitCount );
        const auto prependLength = [=] ( uint32_t length ) { return previousSymbols | ( length << previousBitCount ); };

        auto cacheEntry = inputCacheEntry;
        if ( symbol <= 264U ) {
            cacheEntry.symbols = prependLength( symbol - 257U + 3U + DISTANCE_OFFSET );
            cacheEntry.needToReadDistanceBits = false;
            insertIntoCache( reversedCode, cacheEntry );
        } else if ( symbol < 285U ) {
            const auto lengthCode = static_cast<uint8_t>( symbol - 261U );
            const auto extraBitCount = lengthCode / 4;  /* <= 5 */
            if ( codeLength + extraBitCount <= m_lutBitsCount ) {
                cacheEntry.needToReadDistanceBits = false;
                cacheEntry.bitsToSkip = codeLength + extraBitCount;
                for ( uint8_t extraBits = 0; extraBits < ( 1U << extraBitCount ); ++extraBits ) {
                    cacheEntry.symbols = prependLength( calculateLength( lengthCode ) + extraBits + DISTANCE_OFFSET );
                    insertIntoCache( reversedCode | ( extraBits << codeLength ),
                                     cacheEntry );
                }
            } else {
                cacheEntry.symbols = prependLength( symbol - 254U + DISTANCE_OFFSET );
                insertIntoCache( reversedCode, cacheEntry );
            }
        } else if ( symbol == 285U ) {
            cacheEntry.needToReadDistanceBits = false;
            cacheEntry.symbols = prependLength( 258U + DISTANCE_OFFSET );
            insertIntoCache( reversedCode, cacheEntry );
        } else {
            throw std::logic_error( "Symbol count or symbols bit field is inconsistent!" );
        }
    }

private:
    alignas( 64 ) std::array<CacheEntry, ( 1UL << LUT_BITS_COUNT )> m_codeCache{};

    uint8_t m_lutBitsCount{ LUT_BITS_COUNT };
    uint8_t m_bitsToReadAtOnce{ LUT_BITS_COUNT };
    bool m_needsToBeZeroed{ false };
};
}  // namespace rapidgzip::deflate
