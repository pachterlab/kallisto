#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <zlib.h>

#include <core/VectorView.hpp>

#include "definitions.hpp"
#include "gzip.hpp"


namespace rapidgzip
{
enum class CompressionStrategy : int
{
    DEFAULT             = Z_DEFAULT_STRATEGY,
    FILTERED            = Z_FILTERED,
    RUN_LENGTH_ENCODING = Z_RLE,
    HUFFMAN_ONLY        = Z_HUFFMAN_ONLY,
    FIXED_HUFFMAN       = Z_FIXED,
};


[[nodiscard]] constexpr std::string_view
toString( const CompressionStrategy compressionStrategy )
{
    using namespace std::literals;

    switch ( compressionStrategy )
    {
    case CompressionStrategy::DEFAULT: return "Default"sv;
    case CompressionStrategy::FILTERED: return "Filtered"sv;
    case CompressionStrategy::RUN_LENGTH_ENCODING: return "Run-Length Encoding"sv;
    case CompressionStrategy::HUFFMAN_ONLY: return "Huffman Only"sv;
    case CompressionStrategy::FIXED_HUFFMAN: return "Fixed Huffman"sv;
    }
    return {};
}


enum class ContainerFormat : int
{
    DEFLATE,
    ZLIB,
    GZIP,
};


template<typename ResultContainer = std::vector<uint8_t> >
[[nodiscard]] ResultContainer
compressWithZlib( const VectorView<uint8_t> toCompress,
                  const CompressionStrategy compressionStrategy = CompressionStrategy::DEFAULT,
                  const VectorView<uint8_t> dictionary = {},
                  const ContainerFormat     containerFormat = ContainerFormat::GZIP )
{
    ResultContainer output;
    output.reserve( toCompress.size() );

    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    stream.avail_in = toCompress.size();
    stream.next_in = const_cast<Bytef*>( reinterpret_cast<const Bytef*>( toCompress.data() ) );
    stream.avail_out = 0;
    stream.next_out = nullptr;

    /* > Add 16 to windowBits to write a simple gzip header and trailer around the
     * > compressed data instead of a zlib wrapper.
     * > windowBits can also be -8..-15 for raw deflate. In this case, -windowBits determines the window size. */
    auto windowBits = MAX_WBITS;
    switch ( containerFormat )
    {
    case ContainerFormat::DEFLATE:
        windowBits = -windowBits;
        break;
    case ContainerFormat::ZLIB:
        break;
    case ContainerFormat::GZIP:
        windowBits += 16;
        break;
    }

    deflateInit2( &stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                  windowBits, /* memLevel */ 8, static_cast<int>( compressionStrategy ) );

    if ( !dictionary.empty() ) {
        deflateSetDictionary( &stream, const_cast<Bytef*>( reinterpret_cast<const Bytef*>( dictionary.data() ) ),
                              dictionary.size() );
    }

    auto status = Z_OK;
    constexpr auto CHUNK_SIZE = 1_Mi;
    while ( status == Z_OK ) {
        output.resize( output.size() + CHUNK_SIZE );
        stream.next_out = reinterpret_cast<Bytef*>( output.data() + output.size() - CHUNK_SIZE );
        stream.avail_out = CHUNK_SIZE;
        status = ::deflate( &stream, Z_FINISH );
    }

    deflateEnd( &stream );

    output.resize( stream.total_out );
    output.shrink_to_fit();

    return output;
}


/**
 * This is a small wrapper around zlib. It is able to:
 *  - work on BitReader as input
 *  - start at deflate block offset as opposed to gzip start
 */
class ZlibInflateWrapper
{
public:
    explicit
    ZlibInflateWrapper( gzip::BitReader&& bitReader,
                        const size_t      untilOffset = std::numeric_limits<size_t>::max() ) :
        m_bitReader( std::move( bitReader ) ),
        m_encodedStartOffset( m_bitReader.tell() ),
        m_encodedUntilOffset(
            [untilOffset] ( const auto& size ) {
                return size ? std::min( *size, untilOffset ) : untilOffset;
            } ( m_bitReader.size() )
        )
    {
        initStream();
        if ( inflateInit2( &m_stream, m_windowFlags ) != Z_OK ) {
            throw std::runtime_error( "Probably encountered invalid deflate data!" );
        }
    }

    ~ZlibInflateWrapper()
    {
        inflateEnd( &m_stream );
    }

    void
    initStream();

    void
    refillBuffer();

    void
    setWindow( VectorView<uint8_t> const& window )
    {
        m_setWindowSize = window.size();
        if ( inflateSetDictionary( &m_stream, window.data(), window.size() ) != Z_OK ) {
            throw std::runtime_error( "Failed to set back-reference window in zlib!" );
        }
    }

    /**
     * May return fewer bytes than requested. Only reads one deflate stream per call so that it can return
     * the gzip footer appearing after each deflate stream.
     */
    [[nodiscard]] std::pair<size_t, std::optional<Footer> >
    readStream( uint8_t* output,
                size_t   outputSize );

    [[nodiscard]] size_t
    tellCompressed() const
    {
        return m_bitReader.tell() - getUnusedBits();
    }

    void
    setFileType( FileType fileType )
    {
        m_fileType = fileType;
    }

    /**
     * For legacy reasons, this class is always intended to start decompression at deflate boundaries.
     * The file type will only be used when the end of the deflate stream is reached and there is still data
     * to decode. If there is a header at the beginning, you can call this method with argument "true".
     */
    void
    setStartWithHeader( bool enable )
    {
        m_needToReadHeader = enable;
    }

private:
    [[nodiscard]] size_t
    getUnusedBits() const
    {
        /* > on return inflate() always sets strm->data_type to the
         * > number of unused bits in the last byte taken from strm->next_in, plus 64 if
         * > inflate() is currently decoding the last block in the deflate stream [...]
         * > The number of unused bits may in general be greater than seven, except when bit 7 of
         * > data_type is set, in which case the number of unused bits will be less than eight.
         * @see https://github.com/madler/zlib/blob/09155eaa2f9270dc4ed1fa13e2b4b2613e6e4851/zlib.h#L443C22-L444C66 */
        return m_stream.avail_in * BYTE_SIZE + ( static_cast<uint32_t>( m_stream.data_type ) & 0b11'1111U );
    }

    template<size_t SIZE>
    std::array<std::byte, SIZE>
    readBytes();

    [[nodiscard]] Footer
    readGzipFooter();

    [[nodiscard]] Footer
    readZlibFooter();

    [[nodiscard]] Footer
    readDeflateFooter()
    {
        /* Effectively skip over some bits to align to the next byte. */
        readBytes<0>();
        return Footer{};
    }

    Footer
    readFooter()
    {
        switch ( m_fileType )
        {
        case FileType::NONE:
        case FileType::DEFLATE:
            return readDeflateFooter();
        case FileType::GZIP:
        case FileType::BGZF:
            return readGzipFooter();
        case FileType::ZLIB:
            return readZlibFooter();
        case FileType::BZIP2:
            break;
        }
        throw std::logic_error( "[ZlibInflateWrapper::readFooter] Invalid file type!" );
    }

    void
    readHeader();

private:
    gzip::BitReader m_bitReader;
    const size_t m_encodedStartOffset;
    const size_t m_encodedUntilOffset;
    std::optional<size_t> m_setWindowSize;
    bool m_needToReadHeader{ false };

    /* 2^15 = 32 KiB window buffer and minus signaling raw deflate stream to decode.
     * n in [8,15]
     * -n for raw inflate, not looking for zlib/gzip header and not generating a check value!
     * n + 16 for gzip decoding but not zlib
     * n + 32 for gzip or zlib decoding with automatic detection
     * 0 for automatic window size detection based on the zlib header.
     * We set it to -15 to always force raw deflate decoding so that we can decode the header and footer ourselves. */
    const int m_windowFlags{ -15 };
    z_stream m_stream{};
    /* Loading the whole encoded data (multiple MiB) into memory first and then
     * decoding it in one go is 4x slower than processing it in chunks of 128 KiB! */
    std::array<char, 128_Ki> m_buffer{};

    FileType m_fileType{ FileType::GZIP };
};


inline void
ZlibInflateWrapper::initStream()
{
    m_stream = {};

    m_stream.zalloc = Z_NULL;     /* used to allocate the internal state */
    m_stream.zfree = Z_NULL;      /* used to free the internal state */
    m_stream.opaque = Z_NULL;     /* private data object passed to zalloc and zfree */

    m_stream.avail_in = 0;        /* number of bytes available at next_in */
    m_stream.next_in = Z_NULL;    /* next input byte */

    m_stream.avail_out = 0;       /* remaining free space at next_out */
    m_stream.next_out = Z_NULL;   /* next output byte will go here */
    m_stream.total_out = 0;       /* total amount of bytes read */

    m_stream.msg = nullptr;
}


inline void
ZlibInflateWrapper::refillBuffer()
{
    if ( ( m_stream.avail_in > 0 ) || ( m_bitReader.tell() >= m_encodedUntilOffset ) ) {
        return;
    }

    if ( m_bitReader.tell() % BYTE_SIZE != 0 ) {
        /* This might be called at the very first refillBuffer call when it does not start on a byte-boundary. */
        const auto nBitsToPrime = BYTE_SIZE - ( m_bitReader.tell() % BYTE_SIZE );
        if ( inflatePrime( &m_stream, int( nBitsToPrime ), int( m_bitReader.read( nBitsToPrime ) ) ) != Z_OK ) {
            throw std::runtime_error( "InflatePrime failed!" );
        }
        assert( m_bitReader.tell() % BYTE_SIZE == 0 );
    } else if ( const auto remainingBits = m_encodedUntilOffset - m_bitReader.tell(); remainingBits < BYTE_SIZE ) {
        /* This might be called at the very last refillBuffer call, when it does not end on a byte-boundary. */
        if ( inflatePrime( &m_stream, int( remainingBits ), int( m_bitReader.read( remainingBits ) ) ) != Z_OK ) {
            throw std::runtime_error( "InflatePrime failed!" );
        }
        return;
    }

    /* This reads byte-wise from BitReader. */
    m_stream.avail_in = m_bitReader.read(
        m_buffer.data(), std::min( ( m_encodedUntilOffset - m_bitReader.tell() ) / BYTE_SIZE, m_buffer.size() ) );
    m_stream.next_in = reinterpret_cast<unsigned char*>( m_buffer.data() );
}


[[nodiscard]] inline std::pair<size_t, std::optional<Footer> >
ZlibInflateWrapper::readStream( uint8_t* const output,
                                size_t   const outputSize )
{
    m_stream.next_out = output;
    m_stream.avail_out = outputSize;
    m_stream.total_out = 0;

    if ( m_needToReadHeader ) {
        readHeader();
        m_needToReadHeader = false;
    }

    size_t decodedSize{ 0 };
    /* Do not check for avail_out == 0 here so that progress can still be made on empty blocks as might
     * appear in pigz files or at the end of BGZF files. Note that zlib's inflate should return Z_BUF_ERROR
     * anyway if the output buffer is full. It might be that the check was only here to idiomatically avoid
     * a "while ( true )". */
    while ( true ) {
        refillBuffer();

        const auto oldUnusedBits = getUnusedBits();
        const auto oldTotalOut = m_stream.total_out;

        /* == actual zlib inflate call == */
        const auto errorCode = ::inflate( &m_stream, Z_BLOCK );

        /**
         * > Z_BUF_ERROR if no progress was possible or if there was not enough room in the output
         * > buffer when Z_FINISH is used
         * @see https://github.com/madler/zlib/blob/09155eaa2f9270dc4ed1fa13e2b4b2613e6e4851/zlib.h#L511C68-L513C31
         */
        if ( errorCode == Z_BUF_ERROR ) {
            break;
        }

        if ( ( errorCode != Z_OK ) && ( errorCode != Z_STREAM_END ) ) {
            std::stringstream message;
            message << "[ZlibInflateWrapper][Thread " << std::this_thread::get_id() << "] "
                    << "Decoding failed with error code " << errorCode << " "
                    << ( m_stream.msg == nullptr ? "" : m_stream.msg ) << "! "
                    << "Already decoded " << m_stream.total_out << " B. "
                    << "Read " << formatBits( oldUnusedBits - getUnusedBits() ) << " during the failing isal_inflate "
                    << "from offset " << formatBits( m_bitReader.tell() - oldUnusedBits ) << ". "
                    << "Bit range to decode: [" << m_encodedStartOffset << ", " << m_encodedUntilOffset << "]. "
                    << "BitReader::size: " << m_bitReader.size().value_or( 0 ) << ".";

            if ( m_setWindowSize ) {
                message << " Set window size: " << *m_setWindowSize << " B.";
            } else {
                message << " No window was set.";
            }

        #ifndef NDEBUG
            message << " First bytes: 0x\n";
            const auto oldOffset = m_bitReader.tell();
            m_bitReader.seek( m_bitReader.tell() - oldUnusedBits );
            size_t nPrintedBytes{ 0 };
            for ( size_t offset = m_bitReader.tell();
                  ( !m_bitReader.size() || ( offset < *m_bitReader.size() ) ) && ( nPrintedBytes < 128 );
                  offset += BYTE_SIZE, ++nPrintedBytes )
            {
                if ( ( offset / BYTE_SIZE ) % 16 == 0 ) {
                    message << '\n';
                } else if ( ( offset / BYTE_SIZE ) % 8 == 0 ) {
                    message << ' ';
                }
                message << ' ' << std::setw( 2 ) << std::setfill( '0' ) << std::hex << m_bitReader.read<BYTE_SIZE>();
            }
            m_bitReader.seek( oldOffset );
        #endif

            throw std::runtime_error( std::move( message ).str() );
        }

        if ( decodedSize + m_stream.total_out > outputSize ) {
            throw std::logic_error( "Decoded more than fits into the output buffer!" );
        }

        const auto progressedBits = oldUnusedBits != getUnusedBits();
        const auto progressedOutput = m_stream.total_out != oldTotalOut;

        if ( errorCode == Z_STREAM_END ) {
            if ( ( m_stream.total_out == 0 ) && !progressedBits ) {
                break;
            }
            decodedSize += m_stream.total_out;

            /* If we started with raw deflate, then we also have to skip over the gzip footer.
             * Assuming we are decoding gzip and not zlib or multiple raw deflate streams. */
            std::optional<Footer> footer;
            if ( m_windowFlags < 0 ) {
                footer = readFooter();
                readHeader();
            }

            m_stream.next_out = output + decodedSize;
            m_stream.avail_out = outputSize - decodedSize;

            return { decodedSize, footer };
        }

        if ( !progressedBits && !progressedOutput ) {
            break;
        }
    }

    return { decodedSize + m_stream.total_out, std::nullopt };
}


template<size_t SIZE>
std::array<std::byte, SIZE>
ZlibInflateWrapper::readBytes()
{
    std::array<std::byte, SIZE> buffer{};
    size_t footerSize{ 0 };
    for ( auto stillToRemove = buffer.size(); stillToRemove > 0; ) {
        if ( m_stream.avail_in >= stillToRemove ) {
            std::memcpy( buffer.data() + footerSize, m_stream.next_in, stillToRemove );
            footerSize += stillToRemove;

            m_stream.avail_in -= stillToRemove;
            m_stream.next_in += stillToRemove;
            stillToRemove = 0;
        } else {
            std::memcpy( buffer.data() + footerSize, m_stream.next_in, m_stream.avail_in );
            footerSize += m_stream.avail_in;

            stillToRemove -= m_stream.avail_in;
            m_stream.avail_in = 0;
            refillBuffer();
            if ( m_stream.avail_in == 0 ) {
                throw gzip::BitReader::EndOfFileReached();
            }
        }
    }

    return buffer;
}


inline Footer
ZlibInflateWrapper::readGzipFooter()
{
    const auto footerBuffer = readBytes<8U>();
    gzip::Footer footer{ 0, 0 };

    /* Get CRC32 and size machine-endian-agnostically. */
    for ( auto i = 0U; i < 4U; ++i ) {
        const auto subbyte = static_cast<uint8_t>( footerBuffer[i] );
        footer.crc32 += static_cast<uint32_t>( subbyte ) << ( i * BYTE_SIZE );
    }
    for ( auto i = 0U; i < 4U; ++i ) {
        const auto subbyte = static_cast<uint8_t>( footerBuffer[4U + i] );
        footer.uncompressedSize += static_cast<uint32_t>( subbyte ) << ( i * BYTE_SIZE );
    }

    Footer result;
    result.gzipFooter = footer;
    result.blockBoundary.encodedOffset = tellCompressed();
    result.blockBoundary.decodedOffset = 0;  // Should not be used by the caller.
    return result;
}


inline Footer
ZlibInflateWrapper::readZlibFooter()
{
    const auto footerBuffer = readBytes<4U>();
    zlib::Footer footer;

    /* Get Adler32 machine-endian-agnostically. */
    for ( auto i = 0U; i < 4U; ++i ) {
        const auto subbyte = static_cast<uint8_t>( footerBuffer[i] );
        footer.adler32 += static_cast<uint32_t>( subbyte ) << ( i * BYTE_SIZE );
    }

    Footer result;
    result.zlibFooter = footer;
    result.blockBoundary.encodedOffset = tellCompressed();
    result.blockBoundary.decodedOffset = 0;  // Should not be used by the caller.
    return result;
}


inline int
getZlibWindowBits( FileType fileType,
                   int      windowSize )
{
    switch ( fileType )
    {
    case FileType::NONE:
    case FileType::BZIP2:
        break;
    case FileType::BGZF:
    case FileType::GZIP:
        return 16 + windowSize;
    case FileType::DEFLATE:
        return -windowSize;
    case FileType::ZLIB:
        /* > windowBits can also be zero to request that inflate use the window size in
         * > the zlib header of the compressed stream. */
        return 0;
    }
    throw std::logic_error( "[getZlibWindowBits] Invalid file type!" );
}


/**
 * It really only reads the header and then proceeds to reinitialize the stream for raw deflate decoding so
 * that we can decode the footer ourselves.
 */
inline void
ZlibInflateWrapper::readHeader()
{
    auto* const oldNextOut = m_stream.next_out;

    switch ( m_fileType )
    {
    case FileType::NONE:
    case FileType::BZIP2:
        throw std::logic_error( "[ZlibInflateWrapper::readHeader] Invalid file type!" );
    case FileType::DEFLATE:
        break;
    case FileType::BGZF:
    case FileType::GZIP:
    {
        /* Note that inflateInit and inflateReset set total_out to 0 among other things. */
        if ( inflateReset2( &m_stream, getZlibWindowBits( m_fileType, /* 2^15 buffer */ 15 ) ) != Z_OK ) {
            throw std::logic_error( "Probably encountered invalid gzip header!" );
        }

        gz_header gzipHeader{};
        const auto getHeaderError = inflateGetHeader( &m_stream, &gzipHeader );
        if ( getHeaderError != Z_OK ) {
            std::stringstream message;
            message << "Failed to initialize gzip header structure (error: " << getHeaderError
                    << "). Inconsistent zlib stream object?";
            throw std::logic_error( std::move( message ).str() );
        }

        refillBuffer();
        while ( ( m_stream.avail_in > 0 ) && ( gzipHeader.done == 0 ) ) {
            const auto errorCode = ::inflate( &m_stream, Z_BLOCK );
            if ( errorCode != Z_OK ) {
                /* Even Z_STREAM_END would be unexpected here because we test for avail_in > 0. */
                throw std::runtime_error( "Failed to parse gzip header!" );
            }

            /* > As inflate() processes the gzip stream, head->done is zero until the header
             * > is completed, at which time head->done is set to one.
             * > If a zlib stream is being decoded, then head->done is set to -1. */
            if ( gzipHeader.done != 0 ) {
                break;
            }

            refillBuffer();
        }

        if ( m_stream.next_out != oldNextOut ) {
            throw std::logic_error( "Zlib wrote some output even though we only wanted to read the gzip header!" );
        }
        break;
    }
    case FileType::ZLIB:
    {
        const auto& [header, error] = zlib::readHeader( [this] () {
            return static_cast<uint64_t>( readBytes<1U>().front() );
        } );
        if ( error == Error::END_OF_FILE ) {
            return;
        }
        if ( error != Error::NONE ) {
            std::stringstream message;
            message << "Error reading zlib header: " << toString( error );
            throw std::logic_error( std::move( message ).str() );
        }
        break;
    }
    }

    if ( inflateReset2( &m_stream, m_windowFlags ) != Z_OK ) {
        throw std::logic_error( "Probably encountered invalid gzip header!" );
    }
}
}  // namespace rapidgzip
