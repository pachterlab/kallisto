#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <core/BitManipulation.hpp>
#include <core/Error.hpp>

#include "definitions.hpp"


namespace rapidgzip
{
enum class FileType
{
    NONE,
    BGZF,
    GZIP,
    ZLIB,
    DEFLATE,
    BZIP2,
};


[[nodiscard]] inline const char*
toString( FileType fileType )
{
    switch ( fileType )
    {
    case FileType::NONE:
        return "None";
    case FileType::BGZF:
        return "BGZF";
    case FileType::GZIP:
        return "GZIP";
    case FileType::ZLIB:
        return "ZLIB";
    case FileType::DEFLATE:
        return "DEFLATE";
    case FileType::BZIP2:
        return "BZIP2";
    }
    return "";
}


[[nodiscard]] inline bool
hasCRC32( FileType fileType )
{
    switch ( fileType )
    {
    case FileType::NONE:
    case FileType::BZIP2:
        /* > For example, the CRC32 used in Gzip and Bzip2 use the same polynomial,
         * > but Gzip employs reversed bit ordering, while Bzip2 does not. */
    case FileType::DEFLATE:
    case FileType::ZLIB:
        return false;
    case FileType::BGZF:
    case FileType::GZIP:
        return true;
    }

    std::stringstream message;
    message << "Invalid file type: " << static_cast<int>( fileType );
    throw std::invalid_argument( std::move( message ).str() );
}


namespace gzip
{
/** For this namespace, refer to @see RFC 1952 "GZIP File Format Specification" */

constexpr auto MAGIC_ID1 = 0x1FU;
constexpr auto MAGIC_ID2 = 0x8BU;
constexpr auto MAGIC_COMPRESSION = 0x08U;

/* Note that the byte order is reversed because of the LSB BitReader. */
constexpr auto MAGIC_BYTES_GZIP = 0x08'8B'1FU;

/* This is not a gzip-specific constant. It's such so that the decoder will not try to
 * read the whole file to memory for invalid data. */
constexpr auto MAX_ALLOWED_FIELD_SIZE = 1024 * 1024;


[[nodiscard]] inline std::string
getOperatingSystemName( uint8_t code ) noexcept
{
    switch ( code )
    {
    case   0: return "FAT filesystem (MS-DOS, OS/2, NT/Win32)";
    case   1: return "Amiga";
    case   2: return "VMS (or OpenVMS)";
    case   3: return "Unix";
    case   4: return "VM/CMS";
    case   5: return "Atari TOS";
    case   6: return "HPFS filesystem (OS/2, NT)";
    case   7: return "Macintosh";
    case   8: return "Z-System";
    case   9: return "CP/M";
    case  10: return "TOPS-20";
    case  11: return "NTFS filesystem (NT)";
    case  12: return "QDOS";
    case  13: return "Acorn RISCOS";
    default:
    case 255: return "unknown";
    }
    return std::string( "undefined (" ) + std::to_string( code ) + ")";
}


[[nodiscard]] inline std::string
getExtraFlagsDescription( uint8_t code ) noexcept
{
    switch ( code )
    {
    case   0: return "none";
    case   2: return "compressor used maximum compression, slowest algorithm";
    case   4: return "compressor used fastest algorithm";
    default: break;
    }
    return std::string( "undefined (" ) + std::to_string( code ) + ")";
}


struct Header
{
    uint32_t modificationTime{ 0 };
    uint8_t operatingSystem{ 0 };
    /**
     * 2: compressor used maximum compression, slowest algorithm
     * 4: compressor used fastest algorithm
     */
    uint8_t extraFlags{ 0 };

    bool isLikelyASCII{ false };
    std::optional<std::vector<uint8_t> > extra;
    std::optional<std::string> fileName;
    std::optional<std::string> comment;
    std::optional<uint16_t> crc16;
};


struct Footer
{
    uint32_t crc32{ 0 };
    uint32_t uncompressedSize{ 0 };  // If larger than UINT32_MAX, then contains the modulo.
};


inline std::pair<Header, Error>
readHeader( gzip::BitReader& bitReader )
{
    Header header;

    try {
        bitReader.peek<1>();
    } catch ( const gzip::BitReader::EndOfFileReached& ) {
        return { header, Error::END_OF_FILE };
    }

    try {
        const auto readBytes = bitReader.read<3 * BYTE_SIZE>();
        if ( readBytes != MAGIC_BYTES_GZIP ) {
            return { header, Error::INVALID_GZIP_HEADER };
        }

        const auto flags = bitReader.read<BYTE_SIZE>();
        header.modificationTime = static_cast<uint32_t>( bitReader.read<4 * BYTE_SIZE>() );
        header.extraFlags = static_cast<uint8_t>( bitReader.read<BYTE_SIZE>() );
        header.operatingSystem = static_cast<uint8_t>( bitReader.read<BYTE_SIZE>() );

        header.isLikelyASCII = ( flags & 1U ) != 0U;

        const auto readZeroTerminatedString =
            [&bitReader] () -> std::pair<std::string, Error>
            {
                std::string result;
                for ( size_t i = 0; i < MAX_ALLOWED_FIELD_SIZE; ++i ) {
                    if ( bitReader.eof() ) {
                        return { result, Error::EOF_ZERO_STRING };
                    }

                    const auto toAppend = static_cast<char>( bitReader.read<BYTE_SIZE>() );
                    if ( toAppend == 0 ) {
                        break;
                    }

                    result.push_back( toAppend );
                }
                return { result, Error::NONE };
            };

        if ( ( flags & ( 1U << 2U ) ) != 0 ) {
            const auto length = bitReader.read<16>();
            std::vector<uint8_t> extraData( static_cast<size_t>( length ) );
            for ( auto& extraByte : extraData ) {
                extraByte = static_cast<uint8_t>( bitReader.read<BYTE_SIZE>() );
            }
            header.extra = std::move( extraData );
        }

        Error error = Error::NONE;

        if ( ( flags & ( 1U << 3U ) ) != 0 ) {
            std::tie( header.fileName, error ) = readZeroTerminatedString();
            if ( error != Error::NONE ) {
                return { header, error };
            }
        }

        if ( ( flags & ( 1U << 4U ) ) != 0 ) {
            std::tie( header.comment, error ) = readZeroTerminatedString();
            if ( error != Error::NONE ) {
                return { header, error };
            }
        }

        if ( ( flags & ( 1U << 1U ) ) != 0 ) {
            header.crc16 = bitReader.read<16>();
        }
    } catch ( const gzip::BitReader::EndOfFileReached& ) {
        return { header, Error::INCOMPLETE_GZIP_HEADER };
    }

    return { header, Error::NONE };
}


inline Error
checkHeader( gzip::BitReader& bitReader )
{
    const auto readBytes = bitReader.read<3 * BYTE_SIZE>();
    if ( readBytes != MAGIC_BYTES_GZIP ) {
        return Error::INVALID_GZIP_HEADER;
    }

    const auto flags = bitReader.read<BYTE_SIZE>();
    bitReader.read<4 * BYTE_SIZE>();  // modification time
    bitReader.read<BYTE_SIZE>();  // extra flags
    bitReader.read<BYTE_SIZE>();  // OS identifier

    const auto skipZeroTerminatedString =
        [&bitReader] () -> Error
        {
            for ( size_t i = 0; i < MAX_ALLOWED_FIELD_SIZE; ++i ) {
                if ( bitReader.eof() ) {
                    return Error::EOF_ZERO_STRING;
                }

                const auto toAppend = bitReader.read<BYTE_SIZE>();
                if ( toAppend == 0 ) {
                    break;
                }
            }
            return Error::NONE;
        };

    if ( ( flags & ( 1U << 2U ) ) != 0 ) {
        const auto length = bitReader.read<16>();
        bitReader.seek( static_cast<long long int>( length ) * BYTE_SIZE, SEEK_CUR );
    }

    Error error = Error::NONE;

    if ( ( flags & ( 1U << 3U ) ) != 0 ) {
        error = skipZeroTerminatedString();
        if ( error != Error::NONE ) {
            return error;
        }
    }

    if ( ( flags & ( 1U << 4U ) ) != 0 ) {
        error = skipZeroTerminatedString();
        if ( error != Error::NONE ) {
            return error;
        }
    }

    if ( ( flags & ( 1U << 1U ) ) != 0 ) {
        bitReader.read<16>();  // CRC16
    }

    return Error::NONE;
}


inline Footer
readFooter( gzip::BitReader& bitReader )
{
    if ( bitReader.tell() % BYTE_SIZE != 0 ) {
        bitReader.read( BYTE_SIZE - ( bitReader.tell() % BYTE_SIZE ) );
    }

    Footer footer;
    footer.crc32 = static_cast<uint32_t>( bitReader.read<32>() );
    footer.uncompressedSize = static_cast<uint32_t>( bitReader.read<32>() );
    return footer;
}
}  // namespace gzip
}  // namespace rapidgzip


namespace rapidgzip::zlib
{
enum class CompressionLevel
{
    FASTEST = 0,
    FAST = 1,
    DEFAULT = 2,
    SLOWEST = 3,  /* maximum compression */
};


[[nodiscard]] constexpr const char*
toString( CompressionLevel compressionLevel )
{
    switch ( compressionLevel )
    {
    case CompressionLevel::FASTEST: return "Fastest";
    case CompressionLevel::FAST:    return "Fast";
    case CompressionLevel::DEFAULT: return "Default";
    case CompressionLevel::SLOWEST: return "Slowest";
    }
    return "";
}


struct Header
{
    uint16_t windowSize{ 0 };
    CompressionLevel compressionLevel{ CompressionLevel::DEFAULT };
    uint32_t dictionaryID{ 1 /* ADLER32 of empty data stream */ };
};


struct Footer
{
    uint32_t adler32{ 1 };
};


inline std::pair<Header, Error>
readHeader( const std::function<uint64_t()>& readByte )
{
    Header header;
    bool readPartialHeader{ false };

    try {
        const auto cmf = readByte();
        readPartialHeader = true;
        const auto compressionMethod = cmf & nLowestBitsSet<uint64_t, 4>();
        if ( compressionMethod != /* deflate */ 8 ) {
            return { header, Error::INVALID_GZIP_HEADER };
        }

        /* > For CM = 8, CINFO is the base-2 logarithm of the LZ77 window size, minus eight (CINFO=7 indicates
         * > a 32K window size). Values of CINFO above 7 are not allowed in this version of the specification. */
        const auto compressionInfo = cmf >> 4U;
        if ( compressionInfo > 7 ) {
            return { header, Error::INVALID_GZIP_HEADER };
        }
        header.windowSize = static_cast<uint16_t>( 2U << ( 8U + compressionInfo ) );

        const auto flags = readByte();
        if ( ( ( cmf << 8U ) + flags ) % 31 != 0 ) {
            return { header, Error::INVALID_GZIP_HEADER };
        }

        const auto usesDictionary = ( ( flags >> 5U ) & 1U ) != 0;
        if ( usesDictionary ) {
            header.dictionaryID = 0;
            for ( size_t i = 0; i < 4U; ++i ) {
                header.dictionaryID = ( header.dictionaryID << BYTE_SIZE ) | readByte();
            }
            /* For now, dictionaries are not supported because there is no centralized database for dictionary IDs
             * and no API to set dictionary-ID-to-dictionary-contents mappings. */
            return { header, Error::INVALID_GZIP_HEADER };
        }

        header.compressionLevel = static_cast<CompressionLevel>( ( flags >> 6U ) & 0b11U );
    } catch ( const gzip::BitReader::EndOfFileReached& ) {
        return { header, readPartialHeader ? Error::INCOMPLETE_GZIP_HEADER : Error::END_OF_FILE };
    }

    return { header, Error::NONE };
}


inline std::pair<Header, Error>
readHeader( gzip::BitReader& bitReader )
{
    return readHeader( [&] () { return bitReader.read<BYTE_SIZE>(); } );
}


inline Footer
readFooter( gzip::BitReader& bitReader )
{
    if ( bitReader.tell() % BYTE_SIZE != 0 ) {
        bitReader.read( BYTE_SIZE - ( bitReader.tell() % BYTE_SIZE ) );
    }

    Footer footer;
    footer.adler32 = static_cast<uint32_t>( bitReader.read<32>() );
    return footer;
}
}  // namespace rapidgzip::zlib


namespace rapidgzip
{
struct Footer
{
    /**
     * The blockBoundary used to aid block splitting in order to split after a gzip footer because then the window
     * is known to be empty, which would save space and time.
     * The uncompressed block boundary offset is unambiguous and may even be set to 0, e.g., by the InflateWrappers.
     * The compressed block boundary is more ambiguous. There are three possibilities:
     *  - The end of the preceding deflate block. The footer start is then the next byte-aligned boundary.
     *  - The byte-aligned footer start.
     *  - The byte-aligned footer end, which is the file end or the next gzip stream start.
     *    For gzip, it is exactly FOOTER_SIZE bytes after the footer start.
     * Thoughts about the choice:
     *  - The offset after the footer is more relevant to the intended block splitting improvement.
     *  - The previous deflate block end contains the most information because the other two possible
     *    choices can be derived from it by rounding up and adding FOOTER_SIZE. The inverse is not true.
     *  - The previous block end might be the most stable choice because stopping at that boundary is
     *    already a requirement for using ISA-l without an exact untilOffset. Stopping at the footer end
     *    might not work perfectly and might already have read some of the next block.
     * Currently, the unit tests, test that all possibilities to derive the footer offsets: GzipReader, decodeBlock,
     * decodeBlockWithInflateWrapper with ISA-L or zlib, return the same value.
     * That value is currently the footer end because it seemed easier to implement. This might be subject to
     * change until it is actually used for something (e.g. smarter block splitting).
     * The most complicated to implement but least ambiguous solution would be to add all three boundaries to
     * this struct.
     */
    BlockBoundary blockBoundary;
    gzip::Footer gzipFooter;
    zlib::Footer zlibFooter;
};
}
