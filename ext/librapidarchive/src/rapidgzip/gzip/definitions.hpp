#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <string>

#include <filereader/BitReader.hpp>


namespace rapidgzip
{
namespace gzip
{
/**
 * Using 64-bit instead of 32-bit improved performance by 10% when it was introduced.
 * This might be because of rarer (but longer) refilling of the bit buffer, which might
 * improve pipelining and branch prediction a bit.
 */
using BitReader = BitReader<false, uint64_t>;
}  // namespace gzip


/** Note that this describes bytes in the data format not on the host system, which is CHAR_BIT and might differ. */
constexpr auto BYTE_SIZE = 8U;


/** For this namespace, refer to @see RFC 1951 "DEFLATE Compressed Data Format Specification version 1.3" */
namespace deflate
{
constexpr size_t MAX_WINDOW_SIZE = 32ULL * 1024ULL;
/** This is because the length of an uncompressed block is a 16-bit number. */
constexpr size_t MAX_UNCOMPRESSED_SIZE = std::numeric_limits<uint16_t>::max();
/** This is because the code length alphabet can't encode any higher value and because length 0 is ignored! */
constexpr uint8_t MAX_CODE_LENGTH = 15;

/* Precode constants. */
constexpr uint32_t PRECODE_COUNT_BITS = 4;  // The number of bits to encode the precode count
constexpr uint32_t MAX_PRECODE_COUNT = 19;  // The maximum precode count
constexpr uint32_t PRECODE_BITS = 3;        // The number of bits per precode (code length)
constexpr uint32_t MAX_PRECODE_LENGTH = ( 1U << PRECODE_BITS ) - 1U;
static_assert( MAX_PRECODE_LENGTH == 7 );
alignas( 8 ) static constexpr std::array<uint8_t, MAX_PRECODE_COUNT> PRECODE_ALPHABET = {
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
};

constexpr size_t MAX_LITERAL_OR_LENGTH_SYMBOLS = 286;
/**
 * Note that RFC1951 section 3.2.7 lists the range of HCDIST as 1-32, however section 3.2.6 states that:
 * > distance codes 30-31 will never actually occur in the compressed data.
 * This explains why we define MAX_DISTANCE_SYMBOL_COUNT as 30 instead of 32!
 */
constexpr uint8_t MAX_DISTANCE_SYMBOL_COUNT = 30;
/* next power of two (because binary tree) of MAX_LITERAL_OR_LENGTH_SYMBOLS. This is assuming that all symbols
 * are equally likely to appear, i.e., all codes would be encoded with the same number of bits (9). */
constexpr size_t MAX_LITERAL_HUFFMAN_CODE_COUNT = 512;
constexpr size_t MAX_RUN_LENGTH = 258;

constexpr uint16_t END_OF_BLOCK_SYMBOL = 256;

enum class CompressionType : uint8_t
{
    UNCOMPRESSED    = 0b00,
    FIXED_HUFFMAN   = 0b01,
    DYNAMIC_HUFFMAN = 0b10,
    RESERVED        = 0b11,
};


[[nodiscard]] inline std::string
toString( CompressionType compressionType ) noexcept
{
    switch ( compressionType )
    {
    case CompressionType::UNCOMPRESSED:
        return "Uncompressed";
    case CompressionType::FIXED_HUFFMAN:
        return "Fixed Huffman";
    case CompressionType::DYNAMIC_HUFFMAN:
        return "Dynamic Huffman";
    case CompressionType::RESERVED:
        return "Reserved";
    }
    return "Unknown";
}
}  // namespace deflate


/**
 * Used for GzipReader and IsalInflateWrapper to request preemptive stopping points from the decoder.
 */
enum StoppingPoint : uint32_t
{
    NONE                 = 0U,
    END_OF_STREAM_HEADER = 1U << 0U,
    END_OF_STREAM        = 1U << 1U,  // after gzip footer has been read
    END_OF_BLOCK_HEADER  = 1U << 2U,
    END_OF_BLOCK         = 1U << 3U,
    ALL                  = 0xFFFF'FFFFU,
};


[[nodiscard]] inline std::string
toString( StoppingPoint stoppingPoint )
{
    // *INDENT-OFF*
    switch ( stoppingPoint )
    {
    case StoppingPoint::NONE                 : return "None";
    case StoppingPoint::END_OF_STREAM_HEADER : return "End of Stream Header";
    case StoppingPoint::END_OF_STREAM        : return "End of Stream";
    case StoppingPoint::END_OF_BLOCK_HEADER  : return "End of Block Header";
    case StoppingPoint::END_OF_BLOCK         : return "End of Block";
    case StoppingPoint::ALL                  : return "All";
    };
    return "Unknown";
    // *INDENT-ON*
}


struct BlockBoundary
{
    size_t encodedOffset{ 0 };
    size_t decodedOffset{ 0 };

    [[nodiscard]] bool
    operator==( const BlockBoundary& other ) const
    {
        return ( encodedOffset == other.encodedOffset ) && ( decodedOffset == other.decodedOffset );
    }
};
}  // namespace rapidgzip
