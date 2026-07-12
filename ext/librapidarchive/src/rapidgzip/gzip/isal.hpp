#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <igzip_lib.h>

#include <core/VectorView.hpp>

#include "definitions.hpp"
#include "gzip.hpp"


namespace rapidgzip
{
/**
 * This is a small wrapper around ISA-l. It is able to:
 *  - work on BitReader as input
 *  - start at deflate block offset as opposed to gzip start
 */
class IsalInflateWrapper
{
public:
    using CompressionType = deflate::CompressionType;

public:
    explicit
    IsalInflateWrapper( gzip::BitReader&& bitReader,
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
    }

    void
    initStream();

    void
    refillBuffer();

    void
    setWindow( VectorView<uint8_t> const& window )
    {
        m_setWindowSize = window.size();
        if ( isal_inflate_set_dict( &m_stream, window.data(), window.size() ) != COMP_OK ) {
            throw std::runtime_error( "Failed to set back-reference window in ISA-l!" );
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
    setStoppingPoints( StoppingPoint stoppingPoints )
    {
        m_stream.points_to_stop_at = static_cast<isal_stopping_point>( stoppingPoints );
    }

    [[nodiscard]] StoppingPoint
    stoppedAt() const
    {
        return static_cast<StoppingPoint>( m_stream.stopped_at );
    }

    [[nodiscard]] bool
    isFinalBlock() const
    {
        return m_stream.bfinal != 0;
    }

    [[nodiscard]] std::optional<CompressionType>
    compressionType() const
    {
        if ( stoppedAt() != StoppingPoint::END_OF_BLOCK_HEADER ) {
            return std::nullopt;
        }

        switch ( m_stream.btype )
        {
        case 0: return CompressionType::UNCOMPRESSED;
        case 1: return CompressionType::FIXED_HUFFMAN;
        case 2: return CompressionType::DYNAMIC_HUFFMAN;
        default: break;
        }

        return std::nullopt;
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

    [[nodiscard]] static std::string_view
    getErrorString( int errorCode ) noexcept;

private:
    [[nodiscard]] size_t
    getUnusedBits() const
    {
        return m_stream.avail_in * BYTE_SIZE + m_stream.read_in_length;
    }

    [[nodiscard]] bool
    hasInput() const
    {
        return ( m_stream.avail_in > 0 ) || ( m_stream.read_in_length > 0 );
    }

    void
    inflatePrime( size_t   nBitsToPrime,
                  uint64_t bits )
    {
        m_stream.read_in |= bits << static_cast<uint8_t>( m_stream.read_in_length );
        m_stream.read_in_length += static_cast<int32_t>( nBitsToPrime );
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
        throw std::logic_error( "[IsalInflateWrapper::readFooter] Invalid file type!" );
    }

    [[nodiscard]] bool
    readHeader();

    template<typename Header,
             typename GetHeader>
    bool
    readIsalHeader( Header*          header,
                    const GetHeader& getHeader );

private:
    gzip::BitReader m_bitReader;
    const size_t m_encodedStartOffset;
    const size_t m_encodedUntilOffset;
    std::optional<size_t> m_setWindowSize;

    inflate_state m_stream{};
    /* Loading the whole encoded data (multiple MiB) into memory first and then
     * decoding it in one go is 4x slower than processing it in chunks of 128 KiB! */
    std::array<char, 128_Ki> m_buffer{};

    std::optional<StoppingPoint> m_currentPoint;
    bool m_needToReadHeader{ false };
    FileType m_fileType{ FileType::GZIP };
};


inline void
IsalInflateWrapper::initStream()
{
    isal_inflate_init( &m_stream );
    m_stream.crc_flag = ISAL_DEFLATE;  // This way no gzip header or footer is read
    /* The next_in, avail_in, next_out, avail_out "interface" is the same as zlib. */
    m_stream.next_in = nullptr;
    m_stream.avail_in = 0;
    m_stream.read_in = 0;
    m_stream.read_in_length = 0;
}


inline void
IsalInflateWrapper::refillBuffer()
{
    if ( ( m_stream.avail_in > 0 ) || ( m_bitReader.tell() >= m_encodedUntilOffset ) ) {
        return;
    }

    if ( m_bitReader.tell() % BYTE_SIZE != 0 ) {
        /* This might be called at the very first refillBuffer call when it does not start on a byte-boundary. */
        const auto nBitsToPrime = BYTE_SIZE - ( m_bitReader.tell() % BYTE_SIZE );
        inflatePrime( nBitsToPrime, m_bitReader.read( nBitsToPrime ) );
        assert( m_bitReader.tell() % BYTE_SIZE == 0 );
    } else if ( const auto remainingBits = m_encodedUntilOffset - m_bitReader.tell(); remainingBits < BYTE_SIZE ) {
        /* This might be called at the very last refillBuffer call, when it does not end on a byte-boundary. */
        inflatePrime( remainingBits, m_bitReader.read( remainingBits ) );
        return;
    }

    /* This reads byte-wise from BitReader. */
    m_stream.avail_in = m_bitReader.read(
        m_buffer.data(), std::min( ( m_encodedUntilOffset - m_bitReader.tell() ) / BYTE_SIZE, m_buffer.size() ) );
    m_stream.next_in = reinterpret_cast<unsigned char*>( m_buffer.data() );
}


[[nodiscard]] inline std::pair<size_t, std::optional<Footer> >
IsalInflateWrapper::readStream( uint8_t* const output,
                                size_t   const outputSize )
{
    m_stream.next_out = output;
    m_stream.avail_out = outputSize;
    m_stream.total_out = 0;
    /* isal_inflate also clears it like this but isal_inflate might not be called in this function! */
    m_stream.stopped_at = ISAL_STOPPING_POINT_NONE;

    if ( m_needToReadHeader ) {
        const auto headerSuccess = readHeader();
        if ( !headerSuccess ) {
            return { 0, std::nullopt };
        }
        m_needToReadHeader = false;
        if ( ( m_stream.points_to_stop_at & ISAL_STOPPING_POINT_END_OF_STREAM_HEADER ) != 0 ) {
            m_stream.stopped_at = ISAL_STOPPING_POINT_END_OF_STREAM_HEADER;
            return { 0, std::nullopt };
        }
    }

    size_t decodedSize{ 0 };
    while ( ( decodedSize + m_stream.total_out < outputSize ) && ( m_stream.avail_out > 0 ) ) {
        refillBuffer();
        /* Note that even with avail_in == 0 and read_in_length == 0, there still might be new output from
         * an inflate call because of data in m_stream.tmp_out_buffer. Instead of checking almost non-public
         * members whether progress is possible, simply call inflate and check whether progress was actually done.
         * Of course, the latter assumes will always be done */
        /*if ( ( m_stream.avail_in == 0 ) && ( m_stream.read_in_length == 0 ) && ( m_stream.tmp_out_valid == 0 ) ) {
            break;
        }*/

        /**
         * > If the crc_flag is set to ISAL_GZIP or ISAL_ZLIB, the
         * > gzip/zlib header is parsed, state->crc is set to the appropriate checksum,
         * > and the checksum is verified. If the crc_flag is set to ISAL_DEFLATE
         * > (default), then the data is treated as a raw deflate block.
         * Not that m_stream has a tmp_out_buffer member, which might store some arbitrary amount of decompressed
         * data that will be returned written on the next call possibly without making any progress on the input!
         * Note that in some very rare cases, the call to isal_inflate only moves bytes from avail_in to read_in!
         * This seems to happen right before the EOB symbol. Therefore, do not check the getUnusedBits but instead
         * check the unused bytes and bits separately! For this pathological case, see the test that uses:
         * src/tests/data/wikidata-20220103-all.json.gz-379508635534b--379510732698b.deflate
         */
        const auto oldPosition = std::make_tuple( m_stream.avail_in, m_stream.read_in_length, m_stream.total_out );
        const auto oldUnusedBits = getUnusedBits();  // Only used for debug output.

        /* == actual ISA-L inflate call == */
        const auto errorCode = isal_inflate( &m_stream );

        /* isal_inflate maps all other valid (>= 0) return codes, such as ISAL_OUT_OVERFLOW to
         * ISAL_DECOMP_OK (0). See the code comment in igzip_lib.h. */
        if ( errorCode < 0 ) {
            std::stringstream message;
            message << "[IsalInflateWrapper][Thread " << std::this_thread::get_id() << "] "
                    << "Decoding failed with error code " << errorCode << ": " << getErrorString( errorCode )
                    << "! Already decoded " << m_stream.total_out << " B. "
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
            m_bitReader.seek( m_encodedStartOffset );
            size_t nPrintedBytes{ 0 };
            for ( size_t offset = m_encodedStartOffset;
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

        if ( m_stream.stopped_at != ISAL_STOPPING_POINT_NONE ) {
            break;
        }

        const auto newPosition = std::make_tuple( m_stream.avail_in, m_stream.read_in_length, m_stream.total_out );
        const auto progressed = oldPosition != newPosition;

        if ( m_stream.block_state == ISAL_BLOCK_FINISH ) {
            decodedSize += m_stream.total_out;

            /* If we started with raw deflate, then we also have to skip over the gzip footer.
             * Assuming we are decoding gzip and not zlib or multiple raw deflate streams. */
            const auto footer = readFooter();  // This resets m_stream.total_out
            if ( ( m_stream.points_to_stop_at & ISAL_STOPPING_POINT_END_OF_STREAM ) != 0 ) {
                m_needToReadHeader = true;
                m_stream.stopped_at = ISAL_STOPPING_POINT_END_OF_STREAM;
            } else {
                const auto headerSuccess = readHeader();
                if ( headerSuccess
                     && ( ( m_stream.points_to_stop_at & ISAL_STOPPING_POINT_END_OF_STREAM_HEADER ) != 0 ) )
                {
                    m_stream.stopped_at = ISAL_STOPPING_POINT_END_OF_STREAM_HEADER;
                }
            }

            m_stream.next_out = output + decodedSize;
            m_stream.avail_out = outputSize - decodedSize;

            return { decodedSize, footer };
        }

        if ( !progressed ) {
            break;
        }
    }

    return { decodedSize + m_stream.total_out, std::nullopt };
}


template<size_t SIZE>
std::array<std::byte, SIZE>
IsalInflateWrapper::readBytes()
{
    const auto remainingBits = static_cast<uint8_t>( m_stream.read_in_length % BYTE_SIZE );
    m_stream.read_in >>= remainingBits;
    m_stream.read_in_length -= remainingBits;

    std::array<std::byte, SIZE> buffer{};
    for ( auto stillToRemove = buffer.size(); stillToRemove > 0; ) {
        const auto footerSize = buffer.size() - stillToRemove;
        if ( m_stream.read_in_length > 0 ) {
            /* This should be ensured by making read_in_length % BYTE_SIZE == 0 prior. */
            assert( m_stream.read_in_length >= static_cast<int>( BYTE_SIZE ) );

            buffer[footerSize] = static_cast<std::byte>( m_stream.read_in & 0xFFU );
            m_stream.read_in >>= BYTE_SIZE;
            m_stream.read_in_length -= BYTE_SIZE;
            --stillToRemove;
        } else if ( m_stream.avail_in >= stillToRemove ) {
            std::memcpy( buffer.data() + footerSize, m_stream.next_in, stillToRemove );
            m_stream.avail_in -= stillToRemove;
            m_stream.next_in += stillToRemove;
            stillToRemove = 0;
        } else {
            if ( m_stream.avail_in > 0 ) {
                std::memcpy( buffer.data() + footerSize, m_stream.next_in, m_stream.avail_in );
            }
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
IsalInflateWrapper::readGzipFooter()
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
IsalInflateWrapper::readZlibFooter()
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


inline bool
IsalInflateWrapper::readHeader()
{
    /* Note that isal_inflate_init and isal_inflate_reset set total_out to 0.
     * Unfortunately, isal_inflate_reset also resets read_in and read_in_length to 0 thereby effectively
     * skipping bits! */
    const auto oldConfiguration = m_stream;
    isal_inflate_reset( &m_stream );
    m_stream.crc_flag = ISAL_DEFLATE;  // This way no gzip header or footer is read
    m_stream.points_to_stop_at = oldConfiguration.points_to_stop_at;
    m_stream.read_in           = oldConfiguration.read_in & nLowestBitsSet<uint64_t>( oldConfiguration.read_in_length );
    m_stream.read_in_length    = oldConfiguration.read_in_length;
    m_stream.avail_in          = oldConfiguration.avail_in;
    m_stream.next_in           = oldConfiguration.next_in;

    switch ( m_fileType )
    {
    case FileType::NONE:
    case FileType::BZIP2:
        break;

    case FileType::DEFLATE:
        /* There is no outer header to read. We need to directly read the deflate stream next. */
        return true;

    case FileType::BGZF:
    case FileType::GZIP:
    {
        isal_gzip_header gzipHeader{};
        isal_gzip_header_init( &gzipHeader );
        return readIsalHeader( &gzipHeader, isal_read_gzip_header );
    }

    case FileType::ZLIB:
    {
        const auto& [header, error] = zlib::readHeader( [this] () {
            return static_cast<uint64_t>( readBytes<1U>().front() );
        } );
        if ( error == Error::END_OF_FILE ) {
            return false;
        }
        if ( error != Error::NONE ) {
            std::stringstream message;
            message << "Error reading zlib header: " << toString( error );
            throw std::logic_error( std::move( message ).str() );
        }
        return true;
    }
    }

    throw std::logic_error( "[IsalInflateWrapper::readHeader] Invalid file type!" );
}


template<typename Header,
         typename GetHeader>
inline bool
IsalInflateWrapper::readIsalHeader( Header* const    header,
                                    const GetHeader& getHeader )
{
    auto* const oldNextOut = m_stream.next_out;

    refillBuffer();
    if ( !hasInput() ) {
        return false;
    }

    while ( hasInput() ) {
        const auto errorCode = getHeader( &m_stream, header );
        if ( errorCode == ISAL_DECOMP_OK ) {
            break;
        }

        if ( errorCode != ISAL_END_INPUT ) {
            std::stringstream message;
            message << "Failed to parse gzip/zlib header (" << errorCode << ": " << getErrorString( errorCode ) << ")!";
            throw std::runtime_error( std::move( message ).str() );
        }

        refillBuffer();
    }

    if ( m_stream.next_out != oldNextOut ) {
        throw std::logic_error( "ISA-l wrote some output even though we only wanted to read the gzip header!" );
    }

    return true;
}


[[nodiscard]] inline std::string_view
IsalInflateWrapper::getErrorString( int errorCode ) noexcept
{
    switch ( errorCode )
    {
    case ISAL_DECOMP_OK          /*  0 */ : return "No errors encountered while decompressing";
    case ISAL_END_INPUT          /*  1 */ : return "End of input reached";
    case ISAL_OUT_OVERFLOW       /*  2 */ : return "End of output reached";
    case ISAL_NAME_OVERFLOW      /*  3 */ : return "End of gzip name buffer reached";
    case ISAL_COMMENT_OVERFLOW   /*  4 */ : return "End of gzip comment buffer reached";
    case ISAL_EXTRA_OVERFLOW     /*  5 */ : return "End of extra buffer reached";
    case ISAL_NEED_DICT          /*  6 */ : return "Stream needs a dictionary to continue";
    case ISAL_INVALID_BLOCK      /* -1 */ : return "Invalid deflate block found";
    case ISAL_INVALID_SYMBOL     /* -2 */ : return "Invalid deflate symbol found";
    case ISAL_INVALID_LOOKBACK   /* -3 */ : return "Invalid lookback distance found";
    case ISAL_INVALID_WRAPPER    /* -4 */ : return "Invalid gzip/zlib wrapper found";
    case ISAL_UNSUPPORTED_METHOD /* -5 */ : return "Gzip/zlib wrapper specifies unsupported compress method";
    case ISAL_INCORRECT_CHECKSUM /* -6 */ : return "Incorrect checksum found";
    default: break;
    }
    return "Unknown Error";
}


template<typename ResultContainer = std::vector<uint8_t> >
[[nodiscard]] ResultContainer
compressWithIsal( const VectorView<uint8_t> toCompress,
                  const VectorView<uint8_t> dictionary = {} )
{
    /* > Parameter avail_out must be large enough to fit the entire compressed output.
     * > Max expansion is limited to the input size plus the header size of a stored/raw block.
     * The computation in igzip/igzip_rand_test.c is not much more sophisticated, it basically is
     * 2 * in_size + 300 + extra/trailer bytes (<20).
     * igzip_cli.c uses this for the multi-threaded code! It initializes the output buffer size to 2 * BLOCK_SIZE,
     * while the input is limited to BLOCK_SIZE. */
    ResultContainer compressed( toCompress.size() + 1000 );

    isal_zstream stream{};
    isal_deflate_stateless_init( &stream );

    if ( !dictionary.empty() ) {
        isal_deflate_set_dict( &stream, const_cast<uint8_t*>( dictionary.data() ), dictionary.size() );
    }
    stream.level = 1;
    std::array<uint8_t, ISAL_DEF_LVL1_DEFAULT /* 282624 */> compressionBuffer{};
    stream.level_buf = compressionBuffer.data();
    stream.level_buf_size = compressionBuffer.size();

    stream.next_in = const_cast<uint8_t*>( reinterpret_cast<const uint8_t*>( toCompress.data() ) );
    stream.avail_in = toCompress.size();
    stream.next_out = const_cast<uint8_t*>( reinterpret_cast<const uint8_t*>( compressed.data() ) );
    stream.avail_out = compressed.size();
    stream.gzip_flag = IGZIP_GZIP;
    const auto result = isal_deflate_stateless( &stream );
    if ( result != COMP_OK ) {
        throw std::runtime_error( "Compression failed with error code: " + std::to_string( result ) );
    }
    if ( stream.avail_out >= compressed.size() ) {
        std::stringstream message;
        message << "Something went wrong. Avail_out should be smaller or equal than it was before, but it gew from "
                << formatBytes( compressed.size() ) << " to " << formatBytes( stream.avail_out );
        throw std::logic_error( std::move( message ).str() );
    }
    compressed.resize( compressed.size() - stream.avail_out );
    compressed.shrink_to_fit();

    return compressed;
}


/**
 * If @p decompressedSize is not known, use IsalInflateWrapper and gzip::readHeader (if necessary) instead.
 */
template<typename Container>
[[nodiscard]] Container
inflateWithIsal( const Container& toDecompress,
                 const size_t     decompressedSize,
                 const FileType   fileType = FileType::GZIP )
{
    Container decompressed( decompressedSize );

    inflate_state stream{};
    isal_inflate_init( &stream );

    stream.next_in = const_cast<uint8_t*>( reinterpret_cast<const uint8_t*>( toDecompress.data() ) );
    stream.avail_in = toDecompress.size();
    stream.next_out = const_cast<uint8_t*>( reinterpret_cast<const uint8_t*>( decompressed.data() ) );
    stream.avail_out = decompressed.size();

    switch ( fileType )
    {
    case FileType::BGZF:
    case FileType::GZIP:
    {
        isal_gzip_header header{};
        isal_read_gzip_header( &stream, &header );
        break;
    }
    case FileType::ZLIB:
    {
        isal_zlib_header header{};
        isal_read_zlib_header( &stream, &header );
        break;
    }
    case FileType::DEFLATE:
        break;
    default:
        throw std::invalid_argument( std::string( "Unsupported file type for inflating with ISA-L: " )
                                     + toString( fileType ) );
    }

    const auto result = isal_inflate_stateless( &stream );
    if ( result != ISAL_DECOMP_OK ) {
        std::stringstream message;
        message << "Decompression of " << toDecompress.size() << "B sized vector failed with error code: "
                << IsalInflateWrapper::getErrorString( result ) << " (" << std::to_string( result ) << ")";
        throw std::runtime_error( std::move( message ).str() );
    }
    if ( stream.avail_out > 0 ) {
        std::stringstream message;
        message << "Something went wrong. Decompressed only " << formatBytes( decompressedSize - stream.avail_out )
                << " out of " << formatBytes( decompressedSize ) << " requested!";
        throw std::logic_error( std::move( message ).str() );
    }

    return decompressed;
}
}  // namespace rapidgzip
