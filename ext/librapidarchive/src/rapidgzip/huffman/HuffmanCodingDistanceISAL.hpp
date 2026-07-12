#pragma once

#include <array>
#include <cstdint>
#include <utility>

#include <core/common.hpp>                  // forceinline
#include <core/Error.hpp>
#include <core/VectorView.hpp>
#include <rapidgzip/gzip/definitions.hpp>   // BitReader

#include <igzip_lib.h>


namespace rapidgzip
{
/**
 * A wrapper around the Huffman decoder for distance codes from ISA-l
 */
class HuffmanCodingDistanceISAL
{
public:
    static constexpr auto DIST_LEN = ISAL_DEF_DIST_SYMBOLS;
    static constexpr auto LIT_LEN = ISAL_DEF_LIT_LEN_SYMBOLS;
    static constexpr auto LIT_LEN_ELEMS = 514U;

public:
    [[nodiscard]] Error
    initializeFromLengths( const VectorView<uint8_t>& codeLengths )
    {
        std::array<huff_code, LIT_LEN_ELEMS /* 514 */> dist_huff{};
        std::array<uint16_t, 16> dist_count{};

        /* Decode the lit/len and dist huffman codes using the code huffman code */
        for ( size_t i = 0; i < codeLengths.size(); ++i ) {
            const auto symbol = codeLengths[i];

            dist_count[symbol]++;
            write_huff_code( &dist_huff[i], 0, symbol );
        }

        //std::cerr << "set_codes\n";
        if ( set_codes( dist_huff.data(), LIT_LEN, dist_count.data() ) != ISAL_DECOMP_OK ) {
            m_error = Error::INVALID_HUFFMAN_CODE;
            return m_error;
        }

        //std::cerr << "make_inflate_huff_code_dist\n";

        /* max_dist may also be derived from state->hist_bits for when the ISA-L API user configures
         * a smaller window size than 32 KiB. */
        make_inflate_huff_code_dist( &m_huffmanCode, dist_huff.data(), DIST_LEN, dist_count.data(),
                                     /* max_dist */ DIST_LEN );

        m_error = Error::NONE;
        return Error::NONE;
    }

    [[nodiscard]] bool
    isValid() const
    {
        return m_error == Error::NONE;
    }

private:
    static void
    write_huff_code( huff_code * const huff_code,
                     uint32_t    const code,
                     uint32_t    const length )
    {
        huff_code->code_and_length = code | ( length << 24U );
    }

public:
    /* Decodes the next symbol symbol in in_buffer using the huff code defined by
     * huff_code  and returns the value in next_lits and sym_count */
    [[nodiscard]] forceinline std::optional<uint16_t>
    decode( BitReader& bitReader ) const
    {
        static constexpr auto SMALL_SHORT_SYM_LEN = 9U;
        static constexpr auto SMALL_SHORT_SYM_MASK = ( ( 1U << SMALL_SHORT_SYM_LEN ) - 1U );
        static constexpr auto SMALL_SHORT_CODE_LEN_OFFSET = 11U;
        static constexpr auto SMALL_LONG_CODE_LEN_OFFSET = 10U;
        static constexpr auto SMALL_FLAG_BIT_OFFSET = 10U;
        static constexpr auto SMALL_FLAG_BIT = ( 1U << SMALL_FLAG_BIT_OFFSET );

        static constexpr auto DIST_SYM_LEN = 5U;
        static constexpr auto DIST_SYM_MASK = ( ( 1U << DIST_SYM_LEN ) - 1U );

        auto next_bits = bitReader.peek<ISAL_DECODE_SHORT_BITS /* 10 */>();

        /* next_sym is a possible symbol decoded from next_bits. If bit 15 is 0,
         * next_code is a symbol. Bits 9:0 represent the symbol, and bits 14:10
         * represent the length of that symbol's huffman code. If next_sym is not
         * a symbol, it provides a hint of where the large symbols containing
         * this code are located. Note the hint is at largest the location the
         * first actual symbol in the long code list.*/
        uint32_t next_sym = m_huffmanCode.short_code_lookup[next_bits];
        uint32_t bit_count{ 0 };
        if ( LIKELY( ( next_sym & SMALL_FLAG_BIT ) == 0 ) ) [[likely]] {
            /* Return symbol found if next_code is a complete huffman code
             * and shift in buffer over by the length of the next_code */
            bit_count = next_sym >> SMALL_SHORT_CODE_LEN_OFFSET;
            bitReader.seekAfterPeek( bit_count );
        } else {
            /* If a symbol is not found, do a lookup in the long code list. */
            next_bits = bitReader.peek( ( next_sym - SMALL_FLAG_BIT ) >> SMALL_SHORT_CODE_LEN_OFFSET );
            next_sym = m_huffmanCode.long_code_lookup[( next_sym & SMALL_SHORT_SYM_MASK )
                                                      + ( next_bits >> ISAL_DECODE_SHORT_BITS )];
            bit_count = next_sym >> SMALL_LONG_CODE_LEN_OFFSET;
            bitReader.seekAfterPeek( bit_count );
        }

        if ( bit_count == 0 ) {
            return std::nullopt;
        }
        return next_sym & DIST_SYM_MASK;
    }

private:
    Error m_error{ Error::INVALID_HUFFMAN_CODE };
    inflate_huff_code_small m_huffmanCode;
};
}  // namespace rapidgzip
