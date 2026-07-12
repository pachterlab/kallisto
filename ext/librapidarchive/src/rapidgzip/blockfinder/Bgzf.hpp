#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include <filereader/FileReader.hpp>

#include "Interface.hpp"


namespace rapidgzip::blockfinder
{
/**
 * @see https://www.ietf.org/rfc/rfc1952.txt
 * Each member has the following structure:
 *
 *    +---+---+---+---+---+---+---+---+---+---+
 *    |ID1|ID2|CM |FLG|     MTIME     |XFL|OS | (more-->)
 *    +---+---+---+---+---+---+---+---+---+---+
 *
 * (if FLG.FEXTRA set)
 *
 *    +---+---+=================================+
 *    | XLEN  |...XLEN bytes of "extra field"...| (more-->)
 *    +---+---+=================================+
 *
 * ID1 (IDentification 1)
 * ID2 (IDentification 2)
 *    These have the fixed values ID1 = 31 (0x1f, \037), ID2 = 139
 *    (0x8b, \213), to identify the file as being in gzip format.
 *
 * CM (Compression Method)
 *    This identifies the compression method used in the file.  CM
 *    = 0-7 are reserved.  CM = 8 denotes the "deflate"
 *    compression method, which is the one customarily used by
 *    gzip and which is documented elsewhere.
 *
 * FLG (FLaGs)
 *    This flag byte is divided into individual bits as follows:
 *
 *       bit 0   FTEXT
 *       bit 1   FHCRC
 *       bit 2   FEXTRA
 *       bit 3   FNAME
 *       bit 4   FCOMMENT
 *       bit 5   reserved
 *       bit 6   reserved
 *       bit 7   reserved
 *
 * If the FLG.FEXTRA bit is set, an "extra field" is present in
 * the header, with total length XLEN bytes.  It consists of a
 * series of subfields, each of the form:
 *
 *    +---+---+---+---+==================================+
 *    |SI1|SI2|  LEN  |... LEN bytes of subfield data ...|
 *    +---+---+---+---+==================================+
 *
 * @see http://samtools.github.io/hts-specs/SAMv1.pdf
 *
 * Each BGZF block contains a standard gzip file header with the following standard-compliant extensions:
 *
 *  - The F.EXTRA bit in the header is set to indicate that extra fields are present.
 *  - The extra field used by BGZF uses the two subfield ID values 66 and 67 (ASCII ‘BC’).
 *  - The length of the BGZF extra field payload (field LEN in the gzip specification) is 2 (two bytes of payload).
 *  - The payload of the BGZF extra field is a 16-bit unsigned integer in little endian format.
 *    This integer gives the size of the containing BGZF block minus one.
 *
 * => 10 byte gzip header + 8 bytes FEXTRA field.
 *
 * An end-of-file (EOF) trailer or marker block should be written at the end of BGZF files, so that unintended
 * file truncation can be easily detected. The EOF marker block is a particular empty BGZF block encoded
 * with the default zlib compression level settings, and consists of the following 28 hexadecimal bytes:
 * 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00
 * The presence of this EOF marker at the end of a BGZF file indicates that the immediately following physical
 * EOF is the end of the file as intended by the program that wrote it. Empty BGZF blocks are not otherwise
 * special; in particular, the presence of an EOF marker block does not by itself signal end of file.
 */
class Bgzf final :
    public Interface
{
public:
    using HeaderBytes = std::array<uint8_t, 18>;
    using FooterBytes = std::array<uint8_t, 28>;

    static constexpr FooterBytes BGZF_FOOTER = {
        0x1F, 0x8B, 0x08,                   /* gzip magic bytes */
        0x04,                               /* Flags with FEXTRA set */
        0x00, 0x00, 0x00, 0x00,             /* Modification time (dummy) */
        0x00,                               /* Extra flags */
        0xFF,                               /* Unknown OS */
        0x06, 0x00,                         /* Length of extra field */
        0x42, 0x43, 0x02, 0x00, 0x1B, 0x00, /* Extra field with subfield ID "BC" = 0x42 0x43 */
        0x03,                               /* Fixed Huffman compressed deflate block with final bit set
                                             * and a single EOB character, i.e., no contents. */
        0x00,                               /* Part of EOB (257 == 0b000'0000 (7 bits)) plus byte padding */
        0x00, 0x00, 0x00, 0x00,             /* gzip footer CRC32 */
        0x00, 0x00, 0x00, 0x00              /* gzip footer uncompressed size */
    };

public:
    explicit
    Bgzf( UniqueFileReader fileReader ) :
        m_fileReader( std::move( fileReader ) ),
        m_currentBlockOffset( m_fileReader->tell() )
    {
        HeaderBytes header;
        const auto nBytesRead = m_fileReader->read( reinterpret_cast<char*>( header.data() ), header.size() );
        if ( nBytesRead != header.size() ) {
            throw std::invalid_argument( "Could not read enough data from given file!" );
        }

        if ( !isBgzfHeader( header ) ) {
            throw std::invalid_argument( "Given file does not start with a BGZF header!" );
        }

        /* Check the footer, but only if it does not result in buffering the whole file as in SinglePassReader. */
        if ( m_fileReader->seekable() && m_fileReader->size().has_value() ) {
            FooterBytes footer;
            m_fileReader->seek( -static_cast<long long int>( footer.size() ), SEEK_END );
            const auto nBytesReadFooter = m_fileReader->read( reinterpret_cast<char*>( footer.data() ), footer.size() );
            if ( nBytesReadFooter != footer.size() ) {
                throw std::invalid_argument( "Could not read enough data from given file for BGZF footer!" );
            }

            if ( footer != BGZF_FOOTER ) {
                throw std::invalid_argument( "Given file does not end with a BGZF footer!" );
            }

            m_fileReader->seekTo( m_currentBlockOffset );
        }
    }

    [[nodiscard]] static bool
    isBgzfFile( const UniqueFileReader& file )
    {
        const auto oldPos = file->tell();

        HeaderBytes header;
        const auto nBytesRead = file->read( reinterpret_cast<char*>( header.data() ), header.size() );
        if ( ( nBytesRead != header.size() ) || !isBgzfHeader( header ) ) {
            file->seekTo( oldPos );
            return false;
        }

        /* Check the footer, but only if it does not result in buffering the whole file as in SinglePassReader. */
        if ( file->seekable() && file->size().has_value() ) {
            FooterBytes footer;
            file->seek( -static_cast<long long int>( footer.size() ), SEEK_END );
            const auto nBytesReadFooter = file->read( reinterpret_cast<char*>( footer.data() ), footer.size() );
            if ( ( nBytesReadFooter != footer.size() ) || ( footer != BGZF_FOOTER ) ) {
                file->seekTo( oldPos );
                return false;
            }
        }

        file->seekTo( oldPos );
        return true;
    }

    [[nodiscard]] static bool
    isBgzfHeader( const HeaderBytes& header )
    {
        return    ( header[ 0] == 0x1F )
               && ( header[ 1] == 0x8B )
               && ( header[ 2] == 0x08 )
               && ( ( header[3] & ( 1U << 2U ) ) != 0 )
               && ( header[10] == 0x06 )   // length of extra field is 6B
               && ( header[11] == 0x00 )
               && ( header[12] == 'B'  )   // subfield ID "BC"
               && ( header[13] == 'C'  )
               && ( header[14] == 0x02 )   // subfield length is 2B
               && ( header[15] == 0x00 );
    }

    /**
     * @return Size of the deflate block size - 1!
     *         This does not include the size of the gzip stream header and footer!
     */
    [[nodiscard]] static std::optional<uint16_t>
    getBgzfCompressedSize( const HeaderBytes& header )
    {
        if ( isBgzfHeader( header ) ) {
            /* The two-stepped cast is necessary for correct overflows! Number is in little-endian. */
            return ( static_cast<uint16_t>( header[17] ) << 8U ) + header[16];
        }
        return std::nullopt;
    }

    /**
     * @return offset of deflate block in bits (not the gzip stream offset!).
     */
    [[nodiscard]] size_t
    find() override
    {
        if ( m_currentBlockOffset == std::numeric_limits<size_t>::max() ) {
            return m_currentBlockOffset;
        }

        auto result = ( m_currentBlockOffset + HeaderBytes().size() ) * 8;

        m_fileReader->seekTo( m_currentBlockOffset );
        HeaderBytes header;
        const auto nBytesRead = m_fileReader->read( reinterpret_cast<char*>( header.data() ), header.size() );
        if ( nBytesRead == header.size() ) {
            const auto blockSize = getBgzfCompressedSize( header );
            if ( blockSize ) {
                m_currentBlockOffset += *blockSize + 1;
                const auto fileSize = m_fileReader->size();
                if ( fileSize && ( m_currentBlockOffset >= *fileSize ) ) {
                    m_currentBlockOffset = std::numeric_limits<size_t>::max();
                }
            } else {
                if ( !m_fileReader->eof() ) {
                    std::cerr << "Ignoring all junk data after invalid block offset "
                              << m_currentBlockOffset << " B!\n";
                }
                std::cerr << "Failed to get Bgzf metadata!\n";
                m_currentBlockOffset = std::numeric_limits<size_t>::max();
            }
        } else {
            if ( nBytesRead > 0 ) {
                std::cerr << "Got partial header!\n";
            }
            m_currentBlockOffset = std::numeric_limits<size_t>::max();
        }

        return result;
    }

private:
    const UniqueFileReader m_fileReader;
    size_t m_currentBlockOffset = 0;  /**< in bytes because these are gzip stream offsets */
};
}  // rapidgzip::blockfinder
