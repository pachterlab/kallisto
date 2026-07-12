#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    #include <crc.h>
#endif


namespace rapidgzip
{
/* CRC32 according to RFC 1952 */

/* Size: 1 KiB */
using CRC32LookupTable = std::array<uint32_t, 256>;


static constexpr uint32_t CRC32_GENERATOR_POLYNOMIAL{ 0xEDB88320U };


[[nodiscard]] constexpr CRC32LookupTable
createCRC32LookupTable() noexcept
{
    CRC32LookupTable table{};
    for ( uint32_t n = 0; n < table.size(); ++n ) {
        auto c = static_cast<unsigned long int>( n );
        for ( int j = 0; j < 8; ++j ) {
            if ( ( c & 1UL ) != 0 ) {
                c = CRC32_GENERATOR_POLYNOMIAL ^ ( c >> 1U );
            } else {
                c >>= 1U;
            }
        }
        table[n] = c;
    }
    return table;
}


static constexpr int CRC32_LOOKUP_TABLE_SIZE = 256;

/* a small lookup table: raw data -> CRC32 value to speed up CRC calculation */
alignas( 8 ) constexpr static CRC32LookupTable CRC32_TABLE = createCRC32LookupTable();

[[nodiscard]] constexpr uint32_t
updateCRC32( uint32_t crc,
             uint8_t  data ) noexcept
{
    return ( crc >> 8U ) ^ CRC32_TABLE[( crc ^ data ) & 0xFFU];
}


static constexpr size_t MAX_CRC32_SLICE_SIZE = 64;

/**
 * @see https://ieeexplore.ieee.org/document/4531728
 * @see https://create.stephan-brumme.com/crc32/#slicing-by-16-overview
 * @note LUT[n + 1] contains the CRC32 of a byte steam consisting of n zero-bytes.
 * Size: 64 * 256 * 32 bit = 64 KiB
 */
alignas( 8 ) static constexpr std::array<std::array<uint32_t, 256>, MAX_CRC32_SLICE_SIZE> CRC32_SLICE_BY_N_LUT =
    [] ()
    {
        std::array<std::array<uint32_t, 256>, MAX_CRC32_SLICE_SIZE> lut{};
        lut[0] = CRC32_TABLE;
        for ( size_t i = 0; i < lut[0].size(); ++i ) {
            for ( size_t j = 1; j < lut.size(); ++j ) {
                lut[j][i] = updateCRC32( lut[j - 1][i], 0 );
            }
        }
        return lut;
    }();


template<unsigned int SLICE_SIZE>
[[nodiscard]] uint32_t
crc32SliceByN( uint32_t    crc,
               const char* data,
               size_t      size )
{
    static_assert( SLICE_SIZE % 4 == 0, "Chunk size must be divisible by 4 because of the loop unrolling." );
    static_assert( SLICE_SIZE > 0, "Chunk size must not be 0." );
    static_assert( SLICE_SIZE <= MAX_CRC32_SLICE_SIZE, "Chunk size must not exceed the lookup table size." );

    constexpr const auto& LUT = CRC32_SLICE_BY_N_LUT;
    /* Unrolling by 8 increases speed from 4 GB/s to 4.5 GB/s (+12.5%).
     * Might be CPU-dependent (instruction cache size, ...). */
    #pragma GCC unroll 8
    for ( size_t i = 0; i + SLICE_SIZE <= size; i += SLICE_SIZE ) {
        uint32_t firstDoubleWord{ 0 };
        std::memcpy( &firstDoubleWord, data + i, sizeof( uint32_t ) );
        crc ^= firstDoubleWord;

        alignas( 8 ) std::array<uint8_t, SLICE_SIZE> chunk{};
        std::memcpy( chunk.data(), &crc, sizeof( uint32_t ) );
        std::memcpy( chunk.data() + sizeof( uint32_t ),
                     data + i + sizeof( uint32_t ),
                     SLICE_SIZE - sizeof( uint32_t ) );

        uint32_t result = 0;
        /* Has no effect. I assume it is automatically unrolled with -O3 even without this. */
        #pragma GCC unroll 16
        for ( size_t j = 0; j < SLICE_SIZE; ++j ) {
            result ^= LUT[j][chunk[SLICE_SIZE - 1 - j]];
        }
        crc = result;
    }

    for ( size_t i = size - ( size % SLICE_SIZE ); i < size; ++i ) {
        crc = updateCRC32( crc, data[i] );
    }

    return crc;
}


template<unsigned int SLICE_SIZE = 16>
[[nodiscard]] uint32_t
updateCRC32( const uint32_t    crc,
             const char* const buffer,
             const size_t      size )
{
#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    return ~crc32_gzip_refl( ~crc, reinterpret_cast<const uint8_t*>( buffer ), size );
#else
    return crc32SliceByN<SLICE_SIZE>( crc, buffer, size );
#endif
}


/**
 * @return a(x) multiplied (polynomial multiplication) by b(x) modulo p(x)
 *
 * Example polynomial multiplication: 1101 * 1011 (calculate with x^n and modulo coefficients).
 * (Beware that this ad-hoc notation mixes the usual multiplication of real numbers with polynomial one.
 * @verbatim
 *  1101 * 1011 % p =     |     ( 1*x^3 + 0*x^2 + 1*x^1 + 1*x^0 ) * ( 1*x^3 + 1*x^2 + 0*x^1 + 1*x^0 ) =
 *   1 * ( 1101 << 0 )    |     1 * x^3 * ( 1*x^3 + 1*x^2 + 0*x^1 + 1*x^0 )
 * + 1 * ( 1101 << 1 )    |     1 * x^2 * ( 1*x^3 + 1*x^2 + 0*x^1 + 1*x^0 )
 * + 0 * ( 1101 << 2 )    |     0 * x^1 * ( 1*x^3 + 1*x^2 + 0*x^1 + 1*x^0 )
 * + 1 * ( 1101 << 3 )    |     1 * x^0 * ( 1*x^3 + 1*x^2 + 0*x^1 + 1*x^0 )
 * @endverbatim
 * This can be derived by expanding normal polynomial multiplication. The shift by n is congruent
 * to multiplying the left factor with one of the x^n parts of the polynomial.
 * Note that, in contrast to the example, this function works on the reflected polynomial representation.
 */
[[nodiscard]] constexpr uint32_t
polynomialMultiplyModulo( const uint32_t a,
                          uint32_t       b,
                          const uint32_t p )
{
    uint32_t result = 0;
    for ( auto coefficientPosition = uint32_t( 1 ) << 31U; coefficientPosition > 0; coefficientPosition >>= 1U ) {
        if ( ( a & coefficientPosition ) != 0 ) {
            result ^= b;
        }

        const auto overflows = ( b & 1U ) != 0U;
        b >>= 1U;
        if ( overflows ) {
            b ^= p;  // When it overflows, subtract the divisor / generator polynomial to get the remainder.
        }
    }
    return result;
}


/**
 * The n-th entry in this lookup table caches the result of q(x)^(2^n) % p where q(x) = x^1 is a polynomial.
 * This aids in the computation of x^m % p(x) for arbitrary m and fixed p as given for CRC32.
 * This formula arises as a constant when expanding the CRC32 computation of two "concatenated" numbers.
 * See the comment in @ref combineCRC32.
 */
static constexpr std::array<uint32_t, 32> X2N_LUT =
    [] ()
    {
        std::array<uint32_t, 32> result{};
        result[0] = uint32_t( 1 ) << 30U;  // x^1 (reflected notation)
        for ( size_t n = 1; n < result.size(); ++n ) {
            result[n] = polynomialMultiplyModulo( result[n - 1], result[n - 1], CRC32_GENERATOR_POLYNOMIAL );
        }
        return result;
    }();


/**
 * @return x^n % p(x).
 * In order to avoid n calls to @ref polynomialMultiplyModulo, batches of such multiplication have been
 * precomputed. To cover all possible lengths, the batch size grows exponentially, i.e.,
 * x^2, x^4, x^8, x^16, x^32, ..., x^2^32 have all been precomputed and can be multiplied together
 * in a maximum of 32 calls to get the answer for x^n, e.g., x^18 = ( ( x^16 % p ) * ( x^2 % p ) ) % p.
 */
[[nodiscard]] constexpr uint32_t
xPowerModulo( uint64_t exponent )
{
    auto p = uint32_t( 1 ) << 31U;  // x^0 (reflected notation)
    for ( size_t k = 0; exponent > 0; exponent >>= 1U, k++ ) {
        if ( ( exponent & 1U ) != 0U ) {
            p = polynomialMultiplyModulo( X2N_LUT[k % X2N_LUT.size()], p, CRC32_GENERATOR_POLYNOMIAL );
        }
    }
    return p;
}


/**
 * @return The combined CRC32 given two CRC32s for two subsequent parts of a larger stream.
 */
[[nodiscard]] constexpr uint32_t
combineCRC32( uint32_t crc1,
              uint32_t crc2,
              uint64_t crc32ByteStreamLength )
{
    /**
     * @see http://www.ross.net/crc/download/crc_v3.txt
     * CRC computation is getting the remainder of polynomial division without carry (mod 2 per coefficient)
     * using bitstream as an arbitrary length dividend and the generator polynomial as divisor.
     * For gzip the generator polynomial is 0xEDB88320.
     * Note that the highest bit of the generator polynomial corresponding to x^32 is usually omitted.
     * Addition and subtraction without carry are both equivalent to xor.
     * Multiplication without carry: Consider 1101 * 1011 = ( x^3 + x^2 + x^0 ) * ( x^3 + x^2 + x^0 ) =
     *   1 * ( 1101 << 0 )
     * + 1 * ( 1101 << 1 )
     * + 0 * ( 1101 << 2 )
     * + 1 * ( 1101 << 3 )
     * I.e., for each set bit in the second factor, shift the first factor by the bit position to the left and add/xor
     * all the summands.
     * Division is similar to long division / polynomial division without carry, which means that subtraction is xor,
     * and subtraction is possible iff the minuend is larger, which is the case iff the position of the highest set bit
     * is larger (equivalent to the degree of the polynomial).
     * From the mapping to existing operations like xor and modulo, follow some of the usual mathematical properties
     * like associativity, commutativity of xor and multiplication, affinity.
     * @see https://crypto.stackexchange.com/questions/34011/why-is-crc-said-to-be-linear
     *
     * Let there be two streams (arbitrary length numbers) a and b, then the concatenation of these streams is given
     * as c = ( b << len(a) ) xor a assuming that b contains the higher significant bits.
     * b << len(a) = b * ( 1 << len(a) ).
     * crc(c) := c % p, where p is the generator polynomial and % signifies the carryless polynomial division remainder.
     * Then crc(c) = ( b * ( 1 << len(a) ) xor a ) % p = [(b % p) * (( 1 << len(a) ) % p ) ) % p] xor a % p
     *             = [crc(b) * ((1 << len(a)) % p) % p] xor crc(a)
     * This makes use of the "compatibility with multiplication" property of modular arithmetic.
     *
     * We now need to speed up computation of (1 << len(a)) % p. For that, we can cache the results of powers of
     * two and then combine those results as required. E.g., if len(a) = 129, then we can calculate it as
     * [ ( 1 << 128 ) % p ] * [ ( 1 << 1 ) % p ] % p. In this manner, the complexity is only O(log_2(len(a))).
     * Note that we only need to cache results up to the degree of the polynomial (dop) because
     * ( 1 << ( n * dop + m ) ) % p = ( 1 << ( 1 << m ) ) % p.
     * This is because of the long division and because a xor p xor p = a. (Try to do one a % p division by hand,
     * maybe start with the normal polynomial division then go over to the carryless division).
     */
    return polynomialMultiplyModulo( xPowerModulo( crc32ByteStreamLength * 8 ), crc1, CRC32_GENERATOR_POLYNOMIAL )
                                     ^ ( crc2 & 0xFFFF'FFFFU );
}


class CRC32Calculator
{
public:
    void
    setEnabled( bool enabled ) noexcept
    {
        m_enabled = enabled;
    }

    [[nodiscard]] constexpr bool
    enabled() const noexcept
    {
        return m_enabled;
    }

    void
    reset()
    {
        m_crc32 = ~uint32_t( 0 );
        m_streamSizeInBytes = 0;
    }

    [[nodiscard]] uint32_t
    crc32() const noexcept
    {
        return ~m_crc32;
    }

    [[nodiscard]] uint64_t
    streamSize() const noexcept
    {
        return m_streamSizeInBytes;
    }

    void
    update( const void* data,
            size_t      size )
    {
        if ( enabled() ) {
            m_crc32 = updateCRC32( m_crc32, reinterpret_cast<const char*>( data ), size );
            m_streamSizeInBytes += size;
        }
    }

    /**
     * Throws on error.
     */
    bool  // NOLINT(modernize-use-nodiscard)
    verify( uint32_t crc32ToCompare ) const
    {
        if ( !enabled() || ( crc32() == crc32ToCompare ) ) {
            return true;
        }

        std::stringstream message;
        message << "Mismatching CRC32 (0x" << std::hex << crc32() << " <-> stored: 0x" << crc32ToCompare << ")!";
        throw std::domain_error( std::move( message ).str() );
        return false;
    }

    void
    append( const CRC32Calculator& toAppend )
    {
        if ( m_enabled != toAppend.m_enabled ) {
            return;
        }
        m_crc32 = ~combineCRC32( crc32(), toAppend.crc32(), toAppend.streamSize() );
        m_streamSizeInBytes += toAppend.streamSize();
    }

    void
    prepend( const CRC32Calculator& toPrepend )
    {
        if ( m_enabled != toPrepend.m_enabled ) {
            return;
        }
        m_crc32 = ~combineCRC32( toPrepend.crc32(), crc32(), m_streamSizeInBytes );
        m_streamSizeInBytes += toPrepend.streamSize();
    }

protected:
    uint64_t m_streamSizeInBytes{ 0 };
    uint32_t m_crc32{ ~uint32_t( 0 ) };
    bool m_enabled{ true };
};


[[nodiscard]] inline uint32_t
crc32( const void* const buffer,
       const size_t      size )
{
    return ~updateCRC32<>( ~uint32_t( 0 ), reinterpret_cast<const char*>( buffer ), size );
}
}  // namespace rapidgzip
