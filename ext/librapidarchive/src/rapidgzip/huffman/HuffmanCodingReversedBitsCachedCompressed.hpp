#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <optional>
#include <stdexcept>

#include <huffman/HuffmanCodingSymbolsPerLength.hpp>
#include <rapidgzip/gzip/definitions.hpp>


namespace rapidgzip
{
/**
 * This version uses a large lookup table (LUT) to avoid loops over the BitReader to speed things up a lot.
 * The problem is that the LUT creation can take a while depending on the code lengths.
 * - During initialization, it creates a LUT. The index of that array are a fixed number of bits read from BitReader.
 *   To simplify things, the fixed bits must be larger or equal than the maximum code length.
 *   To fill the LUT, the higher bits the actual codes with shorter lengths are filled with all possible values
 *   and the LUT table result is duplicated for all those values. This process is slow.
 * - During decoding, it reads MAX_CODE_LENGTH bits from the BitReader and uses that value to access the LUT,
 *   which contains the symbol and the actual code length, which is <= MAX_CODE_LENGTH. The BitReader will be seeked
 *   by the actual code length.
 * The "compressed" part of the name references the fact that the symbol and code length are stored not as a pair
 * but in a bit-packed manner in the LUT. This reduces the LUT size by 50% for Symbol = uint16_t
 * (value is uint16_t instead of std::pair<uint8_t, uint16_t> and uint16_t aligns to 2 B, which effectively increases
 * the pair size to 4 B inside the array).
 */
template<typename HuffmanCode,
         uint8_t  MAX_CODE_LENGTH,
         typename Symbol,
         size_t   MAX_SYMBOL_COUNT>
class HuffmanCodingReversedBitsCachedCompressed :
    public HuffmanCodingSymbolsPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>
{
public:
    using BaseType = HuffmanCodingSymbolsPerLength<HuffmanCode, MAX_CODE_LENGTH, Symbol, MAX_SYMBOL_COUNT>;
    using BitCount = typename BaseType::BitCount;
    using CodeLengthFrequencies = typename BaseType::CodeLengthFrequencies;

    static constexpr auto LENGTH_SHIFT = requiredBits( MAX_SYMBOL_COUNT );
    static_assert( MAX_SYMBOL_COUNT <= ( 1UL << LENGTH_SHIFT ), "Not enough free bits to pack length into Symbol!" );

public:
    [[nodiscard]] constexpr Error
    initializeFromLengths( const VectorView<BitCount>& codeLengths )
    {
        if ( const auto errorCode = BaseType::initializeFromLengths( codeLengths );
             errorCode != Error::NONE )
        {
            return errorCode;
        }

        /* Initialize the cache.
         * In benchmarks, this takes 28µs out of ~ 30µs for total initialization.
         * And for decoding 13403 deflate blocks in 5.7s, this makes a total overhead of 0.38s (6.6%).
         * The actual block decoding as opposed to header reading, takes roughly 400µs (total over blocks: 5.3s)
         *  -> This adds up to the observed timings and shows that the header reading is still more than
         *     a magnitude faster and could still do some more setup if it reduces decoding more than that!.
         * So it isn't all that large but also doesn't improve speed by all that much either :(
         * Maybe try smaller lookup table to stay in L1 cache?
         * The test processor, a Ryzen 3900X has
         *   L1 Cache: 64K (per core)
         *   L2 Cache: 512K (per core)
         *   L3 Cache: 64MB (shared)
         * So, theoretically it shouldn't exceed the L1 cache size but who knows. */
        //const auto t0 = now();

        if ( m_needsToBeZeroed ) {
            // Works constexpr
            for ( size_t symbol = 0; symbol < ( 1ULL << this->m_maxCodeLength ); ++symbol ) {
                m_codeCache[symbol] = 0;
            }
        }

        auto codeValues = this->m_minimumCodeValuesPerLevel;
        for ( size_t symbol = 0; static_cast<size_t>( symbol ) < codeLengths.size(); ++symbol ) {
            const auto length = codeLengths[symbol];
            if ( length == 0 ) {
                continue;
            }

            const auto k = length - this->m_minCodeLength;
            const auto code = codeValues[k]++;
            const auto reversedCode = reverseBits( code, length );

            const auto fillerBitCount = this->m_maxCodeLength - length;
            const auto maximumPaddedCode = static_cast<HuffmanCode>(
                reversedCode | static_cast<HuffmanCode>( nLowestBitsSet<HuffmanCode>( fillerBitCount ) << length ) );
            assert( maximumPaddedCode < m_codeCache.size() );
            const auto increment = static_cast<HuffmanCode>( HuffmanCode( 1 ) << length );
            const auto value = static_cast<Symbol>( symbol | static_cast<Symbol>( length << LENGTH_SHIFT ) );
            assert( ( value >> LENGTH_SHIFT ) == length );
            assert( ( value & nLowestBitsSet<Symbol, LENGTH_SHIFT>() ) == symbol );
            for ( auto paddedCode = reversedCode; paddedCode <= maximumPaddedCode; paddedCode += increment ) {
                m_codeCache[paddedCode] = value;
            }
        }

        m_needsToBeZeroed = true;

        //const auto t1 = now();
        //std::cerr << "Creating Huffman LUT took " << duration(t0,t1) << " s\n";
        // Some timings for small.gz: 2.3676e-05 s, 2.186e-05 s, 2.3607e-05 s, 2.1791e-05 s, 2.3606e-05 s, 3.6038e-05 s

        return Error::NONE;
    }

    [[nodiscard]] forceinline std::optional<Symbol>
    decode( gzip::BitReader& bitReader ) const
    {
        try {
            const auto value = bitReader.peek( this->m_maxCodeLength );

            assert( value < m_codeCache.size() );
            auto symbol = m_codeCache[(int)value];
            const auto length = symbol >> LENGTH_SHIFT;
            symbol &= nLowestBitsSet<Symbol, LENGTH_SHIFT>();

            if ( length == 0 ) {
                /* This might happen for non-optimal Huffman trees out of which all except the case of a single
                 * symbol with bit length 1 are forbidden! */
                return std::nullopt;
            }

            bitReader.seekAfterPeek( length );
            return symbol;
        } catch ( const gzip::BitReader::EndOfFileReached& ) {
            /* Should only happen at the end of the file and probably not even there
             * because the gzip footer should be longer than the peek length. */
            return BaseType::decode( bitReader );
        }
    }

private:
    /* Note that Symbol is uint16_t but MAX_SYMBOL_COUNT = 512 only requires 9 bits, i.e., we have 7 unused bits,
     * which can be used to store the code length, which only requires ceil(log2(15)) = 4 bits.
     * This scheme is ~5% faster than storing the length and symbol as a pair probably because of multiple reasons:
     *  - any pair < 64-bit probably has to use some bit shifts anyway so not much more work
     *  - using 8-bit length and 16-bit symbol yields non-aligned access quite frequently
     *  - the space reduction by 33% might improve L1 cache hit rates or cache line utilization. */
    alignas( 8 ) std::array<Symbol, ( 1UL << MAX_CODE_LENGTH )> m_codeCache{};
    bool m_needsToBeZeroed{ false };
};
}  // namespace rapidgzip
