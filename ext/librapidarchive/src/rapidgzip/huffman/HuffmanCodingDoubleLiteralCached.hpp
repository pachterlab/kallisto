#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include <rapidgzip/gzip/definitions.hpp>
#include <rapidgzip/huffman/HuffmanCodingReversedCodesPerLength.hpp>


namespace rapidgzip
{
template<typename HuffmanCode,
         uint8_t  MAX_CODE_LENGTH,
         typename Symbol,
         size_t   MAX_SYMBOL_COUNT>
class HuffmanCodingDoubleLiteralCached :
    public HuffmanCodingReversedCodesPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>
{
public:
    using BaseType = HuffmanCodingReversedCodesPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>;
    using BitCount = typename BaseType::BitCount;

    /* Either ceil(log2(MAX_SYMBOL_COUNT)) or std::numeric_limits<Symbol>::digits - ceil(log2(MAX_CODE_LENGTH)),
     * but the ceil o log2 composition is hard to calculate at compile-time.
     * floor o log2 would be position of first non-zero bit. */
    static constexpr auto LENGTH_SHIFT = 10;

    static constexpr auto NONE_SYMBOL = std::numeric_limits<Symbol>::max();
    static_assert( MAX_SYMBOL_COUNT <= NONE_SYMBOL, "Not enough unused symbols for special none symbol!" );

public:
    /**
     * @note Reusing this struct by calling this method multiple times is allowed. All members will be reinitialized.
     */
    constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        //const auto t0 = now();

        if ( const auto errorCode = BaseType::initializeFromLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        /**
         * Forbid single-symbol Huffman codings for this implementation because:
         *  - this implementation is unable to detect invalid encoded symbols, which only are possible for the
         *    single symbol case.
         *  - this implementation is only used for the literal encoding right now for which a single symbol
         *    makes no sense because then that symbol would have to be the end of block symbol and why should there
         *    be empty dynamic blocks? Those are very space-wasting.
         */
        if ( ( this->m_minCodeLength == 1 ) && ( this->m_maxCodeLength == 1 ) && ( this->m_offsets[1] == 1 ) ) {
            return Error::INVALID_CODE_LENGTHS;
        }

        m_nextSymbol = NONE_SYMBOL;

#if 1
        /**
         * Size and decompressed base64 bandwidths:
         * @verbatim
         * 2 * minCodeLength     : 220.2 <= 220.6 +- 0.7 <= 221.4
         * 2 * minCodeLength + 1 : 252.5 <= 254.3 +- 1.7 <= 256
         * 2 * minCodeLength + 2 : 221.15 <= 221.18 +- 0.05 <= 221.24
         * @endverbatim
         * Urgh, it is difficult to find a stable formula for the optimal double cache size :/.
         * I might need the expected deflate block size as well as take into account the whole code length statistics
         * holistically.
         * E.g., 2 * m_minCodeLength allows for 2 min code length values to be cached.
         * 2 * m_minCodeLength + 1 allows the above **and** combinations of the next most frequent with the most-
         * frequent, which should still be pretty common.
         * 2 * m_minCodeLength + 2 probably has a bad performance because the accounted cases become increasingly
         * rare, e.g., it includes the case of two less common symbols, which are expected exponentially less than
         * a single one.
         * Furthermore, this optimum might change and maybe even become more stable, when implementing two-staged
         * lookup.
         * For 2 * minCodeLength + 1:
         *     readDynamicHuffmanCoding : 0.000170131 s (15.3217 %)
         *     readData                 : 0.00094026 s (84.6783 %)
         * minCodeLength is almost always 6 for base64 data in my tests but maxCodeLength can go up to 9.
         */
        m_cachedBitCount = std::min<uint32_t>( std::max<uint32_t>( this->m_maxCodeLength,
                                                                   2 * this->m_minCodeLength + 1 ),
                                               MAX_CODE_LENGTH );

        /* Measuring the time here instead of before the if above, leads to a 40% performance penalty!!??
         * I.e., measuring more code, yields faster times than measuring only a part of the whole...
         * Modern complex processors and compiler optimizations are fun. */
        /* const auto t0 = now() */

        // Works constexpr
        for ( auto& x : m_doubleCodeCache ) {
            x = NONE_SYMBOL;
        }
        // Fast
        //const auto nElements = m_doubleCodeCache.size() * sizeof( m_doubleCodeCache[0] ) / sizeof(uint64_t);
        //for ( size_t i = 0; i < nElements; ++i ) {
        //    reinterpret_cast<uint64_t*>( m_doubleCodeCache.data() )[i] = 0xFFFF'FFFF'FFFF'FFFFULL;
        //}
        // Slow
        //std::memset( m_doubleCodeCache.data(), 0xFF, m_doubleCodeCache.size() * sizeof( m_doubleCodeCache[0] ) );

        const auto size = this->m_offsets[this->m_maxCodeLength - this->m_minCodeLength + 1];
        auto length = this->m_minCodeLength;
        for ( size_t i = 0; i < size; ) {
            const auto reversedCode = this->m_codesPerLength[i];
            const auto symbol = this->m_symbolsPerLength[i];

            /* Do not greedily decode two symbols at once if the first symbol is a special deflate LZ77 symbol,
             * which will consume some of the next bits! */
            if ( ( length + this->m_minCodeLength > m_cachedBitCount ) || ( symbol >= 256 ) ) {
                const auto fillerBitCount = m_cachedBitCount - length;
                const auto symbolAndLength =
                    static_cast<Symbol>( symbol | ( static_cast<Symbol>( length ) << LENGTH_SHIFT ) );

                for ( uint32_t fillerBits = 0; fillerBits < ( uint32_t( 1 ) << fillerBitCount ); ++fillerBits ) {
                    const auto paddedCode = static_cast<HuffmanCode>( fillerBits << length ) | reversedCode;
                    m_doubleCodeCache[paddedCode * 2] = symbolAndLength;
                    //m_doubleCodeCache[paddedCode * 2 + 1] = NONE_SYMBOL;
                }
            } else {
                auto length2 = this->m_minCodeLength;
                for ( size_t i2 = 0; i2 < size; ) {
                    const auto reversedCode2 = this->m_codesPerLength[i2];
                    const auto symbol2 = this->m_symbolsPerLength[i2];

                    /* Store only one symbol if the Huffman code of the second would be truncated because of the
                     * limited bit count for the cache. */
                    const auto totalLength = static_cast<uint32_t>( length + length2 );
                    if ( totalLength > m_cachedBitCount ) {
                        assert( length <= m_cachedBitCount );
                        const auto paddedCode =
                            static_cast<HuffmanCode>(
                                static_cast<HuffmanCode>( reversedCode2 << length ) | reversedCode )
                            & nLowestBitsSet<HuffmanCode>( m_cachedBitCount );

                        m_doubleCodeCache[paddedCode * 2] =
                            static_cast<Symbol>( symbol | ( static_cast<Symbol>( length ) << LENGTH_SHIFT ) );
                        //m_doubleCodeCache[paddedCode * 2 + 1] = NONE_SYMBOL;
                    } else {
                        const auto fillerBitCount = m_cachedBitCount - totalLength;
                        const auto mergedCode = static_cast<HuffmanCode>( ( reversedCode2 << length ) | reversedCode );
                        const auto symbolAndLength =
                            static_cast<Symbol>( symbol | ( static_cast<Symbol>( totalLength ) << LENGTH_SHIFT ) );

                        /* Using SIMD for this loop actually worsens timings. Probably too short or because of the
                         * necessary code rearrangement for the while condition for the required canonical form. */
                        for ( uint32_t fillerBits = 0; fillerBits < ( 1U << fillerBitCount ); ++fillerBits ) {
                            const auto paddedCode = static_cast<HuffmanCode>( fillerBits << totalLength ) | mergedCode;
                            m_doubleCodeCache[paddedCode * 2] = symbolAndLength;
                            m_doubleCodeCache[paddedCode * 2 + 1] = symbol2;
                        }
                    }

                    ++i2;
                    if ( i2 >= size ) {
                        break;
                    }

                    /* Update the code length to reflect the next loop iteration. */
                    while ( this->m_offsets[length2 - this->m_minCodeLength + 1] == i2 ) {
                        length2++;
                    }
                }
            }

            ++i;
            if ( i >= size ) {
                break;
            }

            /* Update the code length to reflect the next loop iteration. */
            while ( this->m_offsets[length - this->m_minCodeLength + 1] == i ) {
                length++;
            }
        }

#else

        std::array<Symbol, 2> symbols = { NONE_SYMBOL, NONE_SYMBOL };
        size_t symbolSize = 0;

        const std::function<void( uint16_t, uint8_t )> fillCache =
            /**
             * @param mergedCode contains one or multiple Huffman codes in the @ref mergedCodeLength lowest bits.
             * @param symbolsOffset contains the index pointing to the symbol vector in @ref m_bitFieldCacheStorage
             *        represented by mergedCode. This is used to avoid reinserting the same symbol
             *        vector for all garbage bits for further Huffman codes not fitting inside CACHED_BIT_COUNT
             */
            [&] ( uint16_t mergedCode,
                  uint8_t  mergedCodeLength )
            {
                assert( mergedCodeLength <= CACHED_BIT_COUNT );

                const auto emptyBitCount = static_cast<uint8_t>( CACHED_BIT_COUNT - mergedCodeLength );

                if ( ( emptyBitCount < this->m_minCodeLength )
                     || ( symbolSize >= 2 )
                     || ( ( symbolSize > 0 ) && ( symbols[symbolSize - 1] >= 256 ) ) ) {
                    for ( HuffmanCode fillerBits = 1; fillerBits < ( 1UL << emptyBitCount ); ++fillerBits ) {
                        const auto paddedCode =
                            static_cast<HuffmanCode>( ( fillerBits << mergedCodeLength ) | mergedCode );

                        m_doubleCodeCache[paddedCode * 2] = ( mergedCodeLength << LENGTH_SHIFT ) | symbols[0];
                        m_doubleCodeCache[paddedCode * 2 + 1] = symbols[1];
                    }
                    return;
                }

                symbols[symbolSize++] = 0;
                for ( BitCount k = 0; k <= this->m_maxCodeLength - this->m_minCodeLength; ++k ) {
                    for ( HuffmanCode subIndex = 0;
                          subIndex < this->m_offsets[k + 1] - this->m_offsets[k];
                          ++subIndex )
                    {
                        const auto code = static_cast<HuffmanCode>( this->m_minimumCodeValuesPerLevel[k] + subIndex );
                        const auto symbol = this->m_symbolsPerLength[this->m_offsets[k] + subIndex];
                        const auto length = this->m_minCodeLength + k;

                        /* Reverse bits so that lookup is does not have to reverse. */
                        const auto reversedCode = reverseBits( code, length );

                        assert( ( reversedCode & nLowestBitsSet<decltype( reversedCode )>( length ) ) == reversedCode );
                        assert( ( mergedCode & nLowestBitsSet<decltype( mergedCode )>( mergedCodeLength ) )
                                == mergedCode );

                        /* Add to cache or append further Huffman codes recursively. */
                        const auto newMergedCode = static_cast<HuffmanCode>( ( reversedCode << mergedCodeLength )
                                                                             | mergedCode )
                                                   & nLowestBitsSet<HuffmanCode, CACHED_BIT_COUNT>();
                        if ( mergedCodeLength + length <= CACHED_BIT_COUNT ) {
                            symbols[symbolSize - 1] = symbol;
                            fillCache( newMergedCode, mergedCodeLength + length );
                        } else {  //if ( m_bitFieldCache[newMergedCode] == 0 ) {
                            /* The check that m_bitFieldCache[newMergedCode] does not already have data is necessary
                             * because all bits higher than CACHED_BIT_COUNT in newMergedCode are zeroed out, which
                             * can result in duplicates. However, preemptively stopping the iteration when the total
                             * code length has been exceeded would also be wrong because the Huffman tree is kinda
                             * traversed in a depth-first order while we are interested in the breadth-first values up
                             * to a threshold depth.
                             * Because the depth is now limited to 2, this should be relatively rare!
                             */
                            m_doubleCodeCache[2 * newMergedCode] = ( mergedCodeLength << LENGTH_SHIFT ) | symbols[0];
                            m_doubleCodeCache[2 * newMergedCode + 1] = NONE_SYMBOL;
                        }
                    }
                }
                symbolSize--;
            };

        fillCache( 0, 0 );

#endif

        //const auto t1 = now();
        //std::cerr << "Creating Huffman LUT took " << duration(t0,t1) << " s\n";

    #if 0
        std::array<uint16_t, 2 * MAX_CODE_LENGTH> mergedCodeLengthFrequencies = {};
        for ( size_t i = 0; i < m_doubleCodeCache.size(); i += 2 ) {
            mergedCodeLengthFrequencies[m_doubleCodeCache[i] >> LENGTH_SHIFT]++;
        }
        std::cerr << "Merged code length frequencies (out of " << (int)CACHED_BIT_COUNT << " cache key size):\n";
        for ( size_t i = 0; i < mergedCodeLengthFrequencies.size(); ++i ) {
            if ( mergedCodeLengthFrequencies[i] != 0 ) {
                std::cerr << " " << i << ":" << mergedCodeLengthFrequencies[i];
            }
        }
        std::cerr << "\n";
    #endif

        return Error::NONE;
    }

    [[nodiscard]] forceinline std::optional<Symbol>
    decode( gzip::BitReader& bitReader ) const
    {
        if ( m_nextSymbol != NONE_SYMBOL ) {
            const auto result = m_nextSymbol;
            m_nextSymbol = NONE_SYMBOL;
            return result;
        }

        try
        {
            const auto value = bitReader.peek( m_cachedBitCount );

            assert( value < m_doubleCodeCache.size() / 2 );
            /* Casting value to int improves speed ~2% probably because of the relaxed bound-checking because
             * signed integer overflow is undefined behavior (and therefore can be assumed to never happen by
             * the compiler). */
            auto symbol1 = m_doubleCodeCache[(int)value * 2];
            m_nextSymbol = m_doubleCodeCache[(int)value * 2 + 1];
            assert( static_cast<uint32_t>( symbol1 >> LENGTH_SHIFT ) <= m_cachedBitCount );
            bitReader.seekAfterPeek( symbol1 >> LENGTH_SHIFT );
            symbol1 &= nLowestBitsSet<Symbol, LENGTH_SHIFT>();

            return symbol1;
        } catch ( const gzip::BitReader::EndOfFileReached& ) {
            /* Should only happen at the end of the file and probably not even there
             * because the gzip footer should be longer than the peek length. */
            return BaseType::decode( bitReader );
        }
    }

private:
    uint32_t m_cachedBitCount{ 0 };
    mutable Symbol m_nextSymbol = NONE_SYMBOL;
    /* note that Symbol is uint16_t but MAX_SYMBOL_COUNT = 512 only requires 9 bits, i.e., we have 7 unused bits,
     * which can be used to store the code length, which only requires ceil(log2(15)) = 4 bits, or 5 bits
     * because we want to store the code length sum for both symbols in only the first symbol.
     * Using std::array<std::array<Symbol, 2>, ( 1UL << CACHED_BIT_COUNT )> instead of a one-dimensional array
     * with the same size reduces speed for base64.gz by 10%! */
    alignas( 8 ) std::array<Symbol, 2 * ( 1UL << MAX_CODE_LENGTH )> m_doubleCodeCache{};
};
}  // namespace rapidgzip
