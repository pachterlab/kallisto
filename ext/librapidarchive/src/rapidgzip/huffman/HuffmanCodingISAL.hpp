#pragma once

#include <array>
#include <cstdint>
#include <utility>

#include <core/common.hpp>                  // forceinline
#include <core/Error.hpp>
#include <core/VectorView.hpp>
#include <huffman/HuffmanCodingBase.hpp>
#include <rapidgzip/gzip/definitions.hpp>   // BitReader

#include <igzip_lib.h>


namespace rapidgzip::deflate
{
/**
 * A wrapper around the Huffman decoder from ISA-l
 */
class HuffmanCodingISAL  // NOLINT(cppcoreguidelines-pro-type-member-init)
{
public:
    static constexpr auto LIT_LEN_ELEMS = 514U;

    static constexpr auto MAX_LIT_LEN_CODE_LEN = 21U;
    static constexpr auto MAX_LIT_LEN_COUNT = MAX_LIT_LEN_CODE_LEN + 2U;
    static constexpr auto LIT_LEN = ISAL_DEF_LIT_LEN_SYMBOLS;

    static constexpr std::array<uint8_t, 32> len_extra_bit_count = {
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x01, 0x01, 0x01, 0x01, 0x02, 0x02, 0x02, 0x02,
        0x03, 0x03, 0x03, 0x03, 0x04, 0x04, 0x04, 0x04,
        0x05, 0x05, 0x05, 0x05, 0x00, 0x00, 0x00, 0x00
    };

public:
    [[nodiscard]] Error
    initializeFromLengths( const VectorView<uint8_t>& codeLengths )
    {
        m_error = checkHuffmanCodeLengths<MAX_CODE_LENGTH>( codeLengths ) ? Error::NONE : Error::INVALID_CODE_LENGTHS;

        std::array<huff_code, LIT_LEN_ELEMS> lit_and_dist_huff{};
        std::array<uint16_t, MAX_LIT_LEN_COUNT> lit_count{};
        std::array<uint16_t, MAX_LIT_LEN_COUNT> lit_expand_count{};

        /* Decode the lit/len and dist huffman codes using the code huffman code */
        for ( size_t i = 0; i < codeLengths.size(); ++i ) {
            const auto symbol = codeLengths[i];

            lit_count[symbol]++;
            write_huff_code( &lit_and_dist_huff[i], 0, symbol );

            if ( ( symbol != 0 ) && ( i >= 264 /* Lit/Len with no extra bits */ ) ) {
                const auto extra_count = len_extra_bit_count[i - 257U];
                lit_expand_count[symbol]--;
                lit_expand_count[symbol + extra_count] += 1U << extra_count;
            }
        }

        std::array<uint32_t, LIT_LEN_ELEMS + 2> code_list{};  /* The +2 is for the extra codes in the static header */

        if ( set_and_expand_lit_len_huffcode( lit_and_dist_huff.data(), LIT_LEN,
                                              lit_count.data(), lit_expand_count.data(), code_list.data() )
             != ISAL_DECOMP_OK ) {
            m_error = Error::INVALID_HUFFMAN_CODE;
            return m_error;
        }

        make_inflate_huff_code_lit_len( &m_huffmanCode, lit_and_dist_huff.data(), LIT_LEN_ELEMS,
                                        lit_count.data(), code_list.data(), 0 );

        return m_error;
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
        huff_code->code_and_length = code | ( length << 24U );  // NOLINT(cppcoreguidelines-pro-type-union-access)
    }

public:
    /* Decodes the next symbol symbol in in_buffer using the huff code defined by
     * huff_code and returns the value in next_lits and sym_count */
    [[nodiscard]] forceinline std::pair<uint32_t, uint32_t>
    decode( gzip::BitReader& bitReader ) const
    {
        static constexpr auto LARGE_SHORT_SYM_LEN = 25U;
        static constexpr auto LARGE_SHORT_SYM_MASK = ( 1U << LARGE_SHORT_SYM_LEN ) - 1U;
        static constexpr auto LARGE_LONG_SYM_LEN = 10U;
        static constexpr auto LARGE_LONG_SYM_MASK = ( 1U << LARGE_LONG_SYM_LEN ) - 1U;
        static constexpr auto LARGE_FLAG_BIT = 1U << 25U;
        static constexpr auto LARGE_SHORT_CODE_LEN_OFFSET = 28U;
        static constexpr auto LARGE_SYM_COUNT_OFFSET = 26U;
        static constexpr auto LARGE_SYM_COUNT_MASK = ( 1U << 2U ) - 1U;
        static constexpr auto LARGE_SHORT_MAX_LEN_OFFSET = 26U;
        static constexpr auto LARGE_LONG_CODE_LEN_OFFSET = 10U;
        static constexpr auto INVALID_SYMBOL = 0x1FFFU;

        /**
         * I have also tested other fixed peek sizes, such as 24 and 48 but they all were significantly slower.
         * It seems that 32 is most amenable when it comes to refilling the bit buffer.
         * @verbatim
         * peek<32>                                      : 431.67 | 439.15 +- 0.09 | 442.87
         * peek<ISAL_DECODE_LONG_BITS(12)> and peek<...> : 419.05 | 424.46 +- 0.07 | 427.36
         * @endverbatim
         */
        uint64_t nextBits{ 0 };
        try
        {
            nextBits = bitReader.peek<32>();
        } catch ( const gzip::BitReader::EndOfFileReached& exception ) {
            /* This should only happen in the error case or for raw deflate streams because those don't have
             * any footer acting as a kind of buffer to ensure that peek always works. */
            const auto [availableBits, count] = bitReader.peekAvailable();
            if ( count == 0 ) {
                throw exception;
            }
            nextBits = availableBits;
        }

        /* next_sym is a possible symbol decoded from next_bits. If bit 15 is 0,
         * next_code is a symbol. Bits 9:0 represent the symbol, and bits 14:10
         * represent the length of that symbol's huffman code. If next_sym is not
         * a symbol, it provides a hint of where the large symbols containing
         * this code are located. Note the hint is at largest the location the
         * first actual symbol in the long code list.*/
        const auto next12Bits = nextBits & N_LOWEST_BITS_SET_LUT<uint64_t>[ISAL_DECODE_LONG_BITS /* 12 */];
        uint32_t next_sym = m_huffmanCode.short_code_lookup[next12Bits];
        if ( LIKELY( ( next_sym & LARGE_FLAG_BIT ) == 0 ) ) [[likely]] {
            /* Return symbol found if next_code is a complete huffman code
             * and shift in buffer over by the length of the next_code */
            const auto bit_count = next_sym >> LARGE_SHORT_CODE_LEN_OFFSET;
            bitReader.seekAfterPeek( bit_count );

            if ( bit_count == 0 ) {
                next_sym = INVALID_SYMBOL;
            }

            return { next_sym & LARGE_SHORT_SYM_MASK,
                     ( next_sym >> LARGE_SYM_COUNT_OFFSET ) & LARGE_SYM_COUNT_MASK };
        }

        /* If a symbol is not found, do a lookup in the long code list starting from the hint in next_sym.
         * > If bit 15 is set, the i  corresponds to the first DECODE_LOOKUP_SIZE bits of a Huffman code which has
         * > length longer than DECODE_LOOKUP_SIZE. In this case, bits 0 through 8
         * > represent an offset into long_code_lookup table and bits 9 through 12
         * > represent the maximum length of a Huffman code starting with the bits in the index i.
         * Ergo, the maximum length itself is stored in 4 bits and therefore cat at most be 16.
         * And probably, because it never is double-cached in the long table, it can at most be
         * deflate::MAX_CODE_LENGTH = 15.
         * In practice, I have encountered peek sizes of up to 20. So, maybe it also includes the distance count?
         * With a distance code, it should only need up to 13 further bits, so 32 bits in the peek above should
         * still be sufficient. But as I am not 100% sure, I cannot get rid of the fallback to bitReader.peek.
         * It probably will never be called and the if also does not seem to add any measurable overhead, probably
         * because branch prediction is never wrong.
         */
        const auto bitCount = next_sym >> LARGE_SHORT_MAX_LEN_OFFSET;
        if ( LIKELY( bitCount <= 32U ) ) {
            nextBits &= N_LOWEST_BITS_SET_LUT<uint64_t>[bitCount];
        } else {
            nextBits = bitReader.peek( bitCount );
        }
        next_sym = m_huffmanCode.long_code_lookup[( next_sym & LARGE_SHORT_SYM_MASK )
                                                  + ( nextBits >> ISAL_DECODE_LONG_BITS )];
        const auto bit_count = next_sym >> LARGE_LONG_CODE_LEN_OFFSET;
        bitReader.seekAfterPeek( bit_count );

        if ( bit_count == 0 ) {
            next_sym = INVALID_SYMBOL;
        }

        return { next_sym & LARGE_LONG_SYM_MASK, 1 };
    }

private:
    Error m_error{ Error::INVALID_HUFFMAN_CODE };
    inflate_huff_code_large m_huffmanCode;
};
}  // namespace rapidgzip::deflate
