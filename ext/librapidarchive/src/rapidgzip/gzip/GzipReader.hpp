#pragma once

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>

#include <core/AffinityHelpers.hpp>
#include <core/BlockMap.hpp>
#include <core/FileUtils.hpp>
#include <filereader/FileReader.hpp>
#include <filereader/Shared.hpp>
#include <rapidgzip/DecodedDataView.hpp>
#include <rapidgzip/gzip/crc32.hpp>
#include <rapidgzip/gzip/definitions.hpp>
#include <rapidgzip/gzip/deflate.hpp>
#include <rapidgzip/gzip/format.hpp>
#include <rapidgzip/gzip/gzip.hpp>
#include <rapidgzip/IndexFileFormat.hpp>
#include <rapidgzip/WindowMap.hpp>

#ifdef WITH_PYTHON_SUPPORT
    #include <filereader/Python.hpp>
    #include <filereader/Standard.hpp>
#endif


namespace rapidgzip
{
/**
 * A strictly sequential gzip interface that can iterate over multiple gzip streams and of course deflate blocks.
 * It cannot seek back nor is it parallelized but it can be used to implement a parallelization scheme.
 */
class GzipReader :
    public FileReader
{
public:
    using DeflateBlock = typename deflate::Block<>;
    using WriteFunctor = std::function<void ( const void*, uint64_t )>;

public:
    explicit
    GzipReader( UniqueFileReader fileReader ) :
        m_fileReader( ensureSharedFileReader( std::move( fileReader ) ) ),
        m_fileType( [this] () {
            const auto result = determineFileTypeAndOffset( m_fileReader->clone() );
            /* Simply assume GZIP if it cannot be determined to show a more useful error message
             * when trying to read the header. */
            return result.has_value() ? result->first : FileType::GZIP;
        } () ),
        m_bitReader( m_fileReader->clone() )
    {}

    explicit
    GzipReader( const GzipReader& other ) :
        m_fileReader( ensureSharedFileReader( other.m_fileReader->clone() ) ),
        m_fileType( other.m_fileType ),
        m_bitReader( other.m_bitReader ),
        m_currentPosition( other.m_currentPosition ),
        m_atEndOfFile( other.m_atEndOfFile ),
        m_currentDeflateBlock( other.m_currentDeflateBlock ),
        m_lastBlockData( other.m_lastBlockData ),
        m_currentPoint( other.m_currentPoint ),
        m_streamBytesCount( other.m_streamBytesCount ),
        m_offsetInLastBuffers( other.m_offsetInLastBuffers ),
        m_crc32Calculator( other.m_crc32Calculator ),
        m_blockMap( other.m_blockMap ),
        m_windowMap( other.m_windowMap ),
        m_didReadHeader( other.m_didReadHeader )
    {}

#ifdef WITH_PYTHON_SUPPORT
    explicit
    GzipReader( const std::string& filePath ) :
        GzipReader( std::make_unique<StandardFileReader>( filePath ) )
    {}

    explicit
    GzipReader( int fileDescriptor ) :
        GzipReader( std::make_unique<StandardFileReader>( fileDescriptor ) )
    {}

    explicit
    GzipReader( PyObject* pythonObject ) :
        GzipReader( std::make_unique<PythonFileReader>( pythonObject ) )
    {}
#endif

    /* FileReader finals */

    [[nodiscard]] int
    fileno() const final
    {
        throw std::logic_error( "This is a virtual file object, which has no corresponding file descriptor!" );
    }

    [[nodiscard]] bool
    seekable() const final
    {
        return m_blockMap && m_blockMap->finalized() && m_bitReader.seekable();
    }

    void
    close() final
    {
        m_bitReader.close();
    }

    [[nodiscard]] bool
    closed() const final
    {
        return m_bitReader.closed();
    }

    [[nodiscard]] bool
    eof() const final
    {
        return m_atEndOfFile;
    }

    [[nodiscard]] bool
    fail() const final
    {
        throw std::logic_error( "Not implemented!" );
    }

    [[nodiscard]] size_t
    tell() const final
    {
        return m_currentPosition;
    }

    [[nodiscard]] std::optional<size_t>
    size() const final
    {
        if ( m_atEndOfFile ) {
            return m_currentPosition;
        }

        if ( m_blockMap && m_blockMap->finalized() ) {
            return m_blockMap->back().second;
        }

        throw std::invalid_argument( "Can't get stream size when not finished reading at least once!" );
    }

    size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) final
    {
        if ( closed() ) {
            throw std::invalid_argument( "You may not call seek on closed GzipReader!" );
        }

        if ( !m_blockMap ) {
            throw std::invalid_argument( "Block map is empty while seeking in GzipRader!" );
        }

        const auto positiveOffset = effectiveOffset( offset, origin );

        if ( positiveOffset == tell() ) {
            /* This extra check is necessary for empty files! */
            m_atEndOfFile = m_currentPosition >= m_blockMap->back().second;
            return positiveOffset;
        }

        if ( !seekable() ) {
            throw std::invalid_argument( "Cannot seek with non-seekable input or without an index!" );
        }

        if ( !m_windowMap ) {
            throw std::invalid_argument( "Window map is empty while seeking in GzipRader!" );
        }

        if ( !m_currentDeflateBlock.has_value() ) {
            readHeader();
        }

        const auto blockInfo = m_blockMap->findDataOffset( positiveOffset );
        if ( !blockInfo.contains( positiveOffset ) ) {
            throw std::logic_error( "BlockMap returned unwanted block!" );
        }

        const auto window = m_windowMap->get( blockInfo.encodedOffsetInBits );
        if ( window ) {
            auto decompressed = window->decompress();
            m_currentDeflateBlock->reset( *decompressed );
        } else {
            m_currentDeflateBlock->reset();
        }

        m_currentPosition = blockInfo.decodedOffsetInBytes;
        m_atEndOfFile = m_currentPosition >= m_blockMap->back().second;
        m_bitReader.seekTo( blockInfo.encodedOffsetInBits );
        readBlockHeader();
        m_didReadHeader = false;
        read( -1, nullptr, positiveOffset - m_currentPosition );
        return m_currentPosition;
    }

    void
    clearerr() final
    {
        m_bitReader.clearerr();
        m_atEndOfFile = false;
        throw std::invalid_argument( "Not fully tested!" );
    }

    size_t
    read( char*  outputBuffer,
          size_t nBytesToRead ) final
    {
        return read( -1, outputBuffer, nBytesToRead );
    }

    /* Gzip specific methods */

    /**
     * @return number of processed bits of compressed input file stream.
     * @note It's only useful for a rough estimate because of buffering and because deflate is block based.
     */
    [[nodiscard]] size_t
    tellCompressed() const
    {
        return m_bitReader.tell();
    }

    [[nodiscard]] std::optional<StoppingPoint>
    currentPoint() const
    {
        return m_currentPoint;
    }

    [[nodiscard]] const auto&
    currentDeflateBlock() const
    {
        return m_currentDeflateBlock;
    }

    /**
     * @param[out] outputBuffer should at least be large enough to hold @p nBytesToRead bytes
     * @return number of bytes written
     */
    size_t
    read( const int     outputFileDescriptor = -1,
          char* const   outputBuffer         = nullptr,
          const size_t  nBytesToRead         = std::numeric_limits<size_t>::max(),
          StoppingPoint stoppingPoints       = StoppingPoint::NONE )
    {
        const auto writeFunctor =
            [nBytesDecoded = uint64_t( 0 ), outputFileDescriptor, outputBuffer]
            ( const void* const buffer,
              uint64_t    const size ) mutable
            {
                auto* const currentBufferPosition = outputBuffer == nullptr ? nullptr : outputBuffer + nBytesDecoded;
                /**
                 * @note We cannot splice easily here because we don't use std::shared_ptr for the data and therefore
                 *       cannot easily extend the lifetime of the spliced data as necessary. It also isn't as
                 *       important as for the multi-threaded version because decoding is the bottleneck for the
                 *       sequential version.
                 */
                const auto errorCode = writeAll( outputFileDescriptor, currentBufferPosition, buffer, size );
                if ( errorCode != 0 ) {
                    std::stringstream message;
                    message << "Failed to write all bytes because of: " << strerror( errorCode )
                            << " (" << errorCode << ")";
                    throw std::runtime_error( std::move( message ).str() );
                }

                nBytesDecoded += size;
            };

        return read( writeFunctor, nBytesToRead, stoppingPoints );
    }

    size_t
    read( const WriteFunctor& writeFunctor,
          const size_t        nBytesToRead = std::numeric_limits<size_t>::max(),
          const StoppingPoint stoppingPoints = StoppingPoint::NONE )
    {
        size_t nBytesDecoded = 0;

        /* This loop is basically a state machine over m_currentPoint and will process different things
         * depending on m_currentPoint and after each processing step it needs to recheck for EOF!
         * First read metadata so that even with nBytesToRead == 0, the position can be advanced over those. */
        while ( ( hasDataToFlush() || !m_bitReader.eof() ) && !eof() ) {
            if ( !m_currentPoint.has_value() || ( *m_currentPoint == StoppingPoint::END_OF_BLOCK_HEADER ) ) {
                const auto nBytesDecodedInStep = readBlock( writeFunctor, nBytesToRead - nBytesDecoded );

                nBytesDecoded += nBytesDecodedInStep;
                m_streamBytesCount += nBytesDecodedInStep;

                /* After this call to readBlock, m_currentPoint is either unchanged END_OF_BLOCK_HEADER,
                 * std::nullopt (block not fully read) or END_OF_BLOCK. In the last case, we should try to read
                 * possible gzip footers and headers even if we already have the requested amount of bytes. */

                if ( !m_currentPoint.has_value() || ( *m_currentPoint == StoppingPoint::END_OF_BLOCK_HEADER ) ) {
                    if ( nBytesDecoded >= nBytesToRead ) {
                        break;
                    }

                    if ( nBytesDecodedInStep == 0 ) {
                        /* We did not advance after the readBlock call and did not even read any amount of bytes.
                         * Something went wrong with flushing. Break to avoid infinite loop. */
                        break;
                    }
                }
            } else {
                /* This else branch only handles headers and footers and will always advance
                 * the current point while not actually decoding any bytes. */
                switch ( *m_currentPoint )
                {
                case StoppingPoint::NONE:
                case StoppingPoint::END_OF_STREAM:
                    readHeader();
                    break;

                case StoppingPoint::END_OF_STREAM_HEADER:
                case StoppingPoint::END_OF_BLOCK:
                    if ( m_currentDeflateBlock.has_value() && m_currentDeflateBlock->eos() ) {
                        readFooter();
                    } else {
                        readBlockHeader();
                    }
                    break;

                // NOLINTBEGIN(bugprone-branch-clone)
                case StoppingPoint::END_OF_BLOCK_HEADER:
                    assert( false && "Should have been handled before the switch!" );
                    break;

                case StoppingPoint::ALL:
                    assert( false && "Should only be specified by the user not appear internally!" );
                    break;
                // NOLINTEND(bugprone-branch-clone)
                }
            }

        #ifdef WITH_PYTHON_SUPPORT
            checkPythonSignalHandlers();
        #endif

            if ( m_currentPoint.has_value() && testFlags( *m_currentPoint, stoppingPoints ) ) {
                break;
            }
        }

        if ( !hasDataToFlush() && m_bitReader.eof() ) {
            m_atEndOfFile = true;
        }

        m_currentPosition += nBytesDecoded;
        return nBytesDecoded;
    }

    void
    setCRC32Enabled( bool enabled )
    {
        m_crc32Calculator.setEnabled( enabled );
    }

    void
    importIndex( UniqueFileReader&& indexFile, size_t parallelization = 0 )
    {
        setBlockOffsets( readGzipIndex( std::move( indexFile ), m_fileReader->clone(),
                                        parallelization == 0 ? availableCores() : parallelization ) );
    }

private:
    [[nodiscard]] UniqueFileReader
    cloneRaw() const override
    {
        return std::make_unique<GzipReader>( *this );
    }

    /**
     * @note Only to be used by readBlock!
     * @return The number of actually flushed bytes, which might be hindered,
     *         e.g., if the output file descriptor can't be written to!
     */
    size_t
    flushOutputBuffer( const WriteFunctor& writeFunctor,
                       size_t              maxBytesToFlush );

    [[nodiscard]] bool
    hasDataToFlush() const
    {
        return m_offsetInLastBuffers && m_currentDeflateBlock
               && m_currentDeflateBlock->isValid()
               && ( *m_offsetInLastBuffers < m_lastBlockData.size() );
    }

    void
    readBlockHeader()
    {
        if ( !m_currentDeflateBlock.has_value() ) {
            throw std::logic_error( "Call readHeader before calling readBlockHeader!" );
        }
        const auto error = m_currentDeflateBlock->readHeader( m_bitReader );
        if ( error != Error::NONE ) {
            std::stringstream message;
            message << "Encountered error: " << toString( error ) << " while trying to read deflate header!";
            throw std::domain_error( std::move( message ).str() );
        }
        m_currentPoint = StoppingPoint::END_OF_BLOCK_HEADER;
    }

    /**
     * Decodes data from @ref m_currentDeflateBlock and writes it to the file descriptor and/or the output buffer.
     * It will either return when the whole block has been read or when the requested amount of bytes has been read.
     */
    [[nodiscard]] size_t
    readBlock( const WriteFunctor& writeFunctor,
               size_t              nMaxBytesToDecode );

    void
    readHeader()
    {
        switch ( m_fileType )
        {
        case FileType::NONE:
        case FileType::BGZF:
        case FileType::GZIP:
        {
            auto [header, error] = gzip::readHeader( m_bitReader );
            if ( error != Error::NONE ) {
                std::stringstream message;
                message << "Encountered error: " << toString( error ) << " while trying to read gzip header!";
                throw std::domain_error( std::move( message ).str() );
            }
            break;
        }
        case FileType::ZLIB:
        {
            auto [header, error] = zlib::readHeader( m_bitReader );
            if ( error != Error::NONE ) {
                std::stringstream message;
                message << "Encountered error: " << toString( error ) << " while trying to read zlib header!";
                throw std::domain_error( std::move( message ).str() );
            }
            break;
        }
        case FileType::DEFLATE:
            break;
        case FileType::BZIP2:
            throw std::domain_error( "Bzip2 not supported by this class!" );
        }

        m_currentDeflateBlock.emplace();
        m_currentDeflateBlock->setInitialWindow();
        m_streamBytesCount = 0;
        m_currentPoint = StoppingPoint::END_OF_STREAM_HEADER;
        m_crc32Calculator.reset();
        m_didReadHeader = true;
    }

    void
    readFooter();

    [[nodiscard]] bool
    bufferHasBeenFlushed() const
    {
        return !m_offsetInLastBuffers.has_value();
    }

    [[nodiscard]] bool
    endOfStream() const
    {
        return !m_currentDeflateBlock.has_value() || !m_currentDeflateBlock->isValid()
               || ( bufferHasBeenFlushed() && m_currentDeflateBlock->eos() );
    }

    void
    setBlockOffsets( const GzipIndex& index )
    {
        if ( index.checkpoints.empty() || !index.windows ) {
            return;
        }

        const auto lockedWindows = index.windows->data();
        if ( lockedWindows.second == nullptr ) {
            throw std::invalid_argument( "Index window map must be a valid pointer!" );
        }

        const auto lessOffset =
            [] ( const auto& a, const auto& b ) {
                return a.uncompressedOffsetInBytes < b.uncompressedOffsetInBytes;
            };
        if ( !std::is_sorted( index.checkpoints.begin(), index.checkpoints.end(), lessOffset ) ) {
            throw std::invalid_argument( "Index checkpoints must be sorted by uncompressed offsets!" );
        }

        if ( index.hasLineOffsets ) {
            throw std::invalid_argument( "Index with line offsets is not supported!" );
        }

        /* Generate simple compressed to uncompressed offset map from index. */
        std::map<size_t, size_t> newBlockOffsets;
        auto windowMap = std::make_shared<WindowMap>();
        for ( const auto& checkpoint: index.checkpoints ) {
            newBlockOffsets.emplace( checkpoint.compressedOffsetInBits, checkpoint.uncompressedOffsetInBytes );
            if ( const auto window = lockedWindows.second->find( checkpoint.compressedOffsetInBits );
                 window != lockedWindows.second->end() )
            {
                windowMap->emplaceShared( checkpoint.compressedOffsetInBits, window->second );
            }
        }

        /* Input file-end offset if not included in checkpoints. */
        if ( const auto fileEndOffset = newBlockOffsets.find( index.compressedSizeInBytes * 8 );
             fileEndOffset == newBlockOffsets.end() )
        {
            newBlockOffsets.emplace( index.compressedSizeInBytes * 8, index.uncompressedSizeInBytes );
            windowMap->emplace( index.compressedSizeInBytes * 8, {}, CompressionType::NONE );
        } else if ( fileEndOffset->second != index.uncompressedSizeInBytes ) {
            throw std::invalid_argument( "Index has contradicting information for the file end information!" );
        }
        m_windowMap = std::move( windowMap );

        setBlockOffsets( newBlockOffsets );
    }

    void
    setBlockOffsets( const std::map<size_t, size_t>& offsets )
    {
        if ( offsets.empty() ) {
            if ( !m_blockMap || m_blockMap->dataBlockCount() == 0 ) {
                return;
            }
            throw std::invalid_argument( "May not clear offsets. Construct a new GzipReader instead!" );
        }

        if ( offsets.size() < 2 ) {
            throw std::invalid_argument( "Block offset map must contain at least one valid block and one EOS block!" );
        }
        auto blockMap = std::make_shared<BlockMap>();
        blockMap->setBlockOffsets( offsets );
        m_blockMap = std::move( blockMap );
    }

private:
    const std::unique_ptr<SharedFileReader> m_fileReader;
    const FileType m_fileType;

    gzip::BitReader m_bitReader;

    size_t m_currentPosition{ 0 };  /** the current position as can only be modified with read or seek calls. */
    bool m_atEndOfFile{ false };

    /**
     * The deflate block will be reused during a gzip stream because each block depends on the last output
     * of the previous block. But after the gzip stream end, this optional will be cleared and in case of
     * another concatenated gzip stream, it will be created anew.
     */
    std::optional<DeflateBlock> m_currentDeflateBlock;
    /** Holds non-owning views to the data decoded in the last call to m_currentDeflateBlock.read. */
    deflate::DecodedDataView m_lastBlockData;

    /**
     * If m_currentPoint has no value, then it means it is currently inside a deflate block.
     * Because a gzip file can contain multiple streams, the file beginning can generically be treated
     * as being at the end of a previous (empty) stream.
     * m_currentPoint may only every have exactly one StoppingPoint set, it may not contain or'ed values!
     */
    std::optional<StoppingPoint> m_currentPoint{ END_OF_STREAM };

    size_t m_streamBytesCount{ 0 };

    /* These are necessary states to return partial results and resume returning further ones.
     * I.e., things which would not be necessary with coroutines supports. This optional has no value
     * iff there is no current deflate block or if we have read all data from it already. */
    std::optional<size_t> m_offsetInLastBuffers;

    CRC32Calculator m_crc32Calculator;

    std::shared_ptr<const BlockMap> m_blockMap;
    std::shared_ptr<const WindowMap> m_windowMap;
    bool m_didReadHeader{ false };
};


inline size_t
GzipReader::flushOutputBuffer( const WriteFunctor& writeFunctor,
                               const size_t        maxBytesToFlush )
{
    if ( !m_offsetInLastBuffers.has_value()
         || !m_currentDeflateBlock.has_value()
         || !m_currentDeflateBlock->isValid() ) {
        return 0;
    }

    size_t totalBytesFlushed = 0;
    size_t bufferOffset = 0;
    for ( const auto& buffer : m_lastBlockData.data ) {
        if ( ( *m_offsetInLastBuffers >= bufferOffset ) && ( *m_offsetInLastBuffers < bufferOffset + buffer.size() ) ) {
            const auto offsetInBuffer = *m_offsetInLastBuffers - bufferOffset;
            const auto nBytesToWrite = std::min( buffer.size() - offsetInBuffer, maxBytesToFlush - totalBytesFlushed );

            m_crc32Calculator.update( reinterpret_cast<const char*>( buffer.data() + offsetInBuffer ), nBytesToWrite );

            if ( writeFunctor ) {
                writeFunctor( buffer.data() + offsetInBuffer, nBytesToWrite );
            }

            *m_offsetInLastBuffers += nBytesToWrite;
            totalBytesFlushed += nBytesToWrite;
        }

        bufferOffset += buffer.size();
    }

    /* Reset optional offset if end has been reached. */
    size_t totalBufferSize = 0;
    for ( const auto& buffer : m_lastBlockData.data ) {
        totalBufferSize += buffer.size();
    }
    if ( m_offsetInLastBuffers >= totalBufferSize ) {
        m_offsetInLastBuffers = std::nullopt;
    }

    return totalBytesFlushed;
}


inline size_t
GzipReader::readBlock( const WriteFunctor& writeFunctor,
                       const size_t        nMaxBytesToDecode )
{
    if ( eof() || ( nMaxBytesToDecode == 0 ) ) {
        return 0;
    }

    /* Try to flush remnants in output buffer from interrupted last call. */
    size_t nBytesDecoded = flushOutputBuffer( writeFunctor, nMaxBytesToDecode );
    if ( !bufferHasBeenFlushed() ) {
        return nBytesDecoded;
    }

    while ( true ) {
        if ( bufferHasBeenFlushed() ) {
            if ( !m_currentDeflateBlock.has_value() || !m_currentDeflateBlock->isValid() ) {
                throw std::logic_error( "Call readHeader and readBlockHeader before calling readBlock!" );
            }

            if ( m_currentDeflateBlock->eob() ) {
                m_currentPoint = StoppingPoint::END_OF_BLOCK;
                return nBytesDecoded;
            }

            /* Decode more data from current block. It can then be accessed via Block::lastBuffers. */
            auto error = Error::NONE;
            std::tie( m_lastBlockData, error ) =
                m_currentDeflateBlock->read( m_bitReader, std::numeric_limits<size_t>::max() );
            if ( error != Error::NONE ) {
                std::stringstream message;
                message << "Encountered error: " << toString( error ) << " while decompressing deflate block.";
                throw std::domain_error( std::move( message ).str() );
            }

            if ( ( m_lastBlockData.size() == 0 ) && !m_currentDeflateBlock->eob() ) {
                throw std::logic_error( "Could not read anything so it should be the end of the block!" );
            }
            m_offsetInLastBuffers = 0;
        }

        if ( nBytesDecoded >= nMaxBytesToDecode ) {
            break;
        }

        m_currentPoint = {};

        const auto flushedCount = flushOutputBuffer( writeFunctor, nMaxBytesToDecode - nBytesDecoded );

        if ( ( flushedCount == 0 ) && !bufferHasBeenFlushed() ) {
            /* Something went wrong with flushing and this would lead to an infinite loop. */
            break;
        }
        nBytesDecoded += flushedCount;
    }

    return nBytesDecoded;
}


inline void
GzipReader::readFooter()
{
    switch ( m_fileType )
    {
    case FileType::NONE:
    case FileType::BGZF:
    case FileType::GZIP:
    {
        const auto footer = gzip::readFooter( m_bitReader );

        if ( m_didReadHeader && static_cast<uint32_t>( m_streamBytesCount ) != footer.uncompressedSize ) {
            std::stringstream message;
            message << "Mismatching size (" << static_cast<uint32_t>( m_streamBytesCount ) << " <-> footer: "
                    << footer.uncompressedSize << ") for gzip stream!";
            throw std::domain_error( std::move( message ).str() );
        }

        if ( !m_currentDeflateBlock.has_value() || !m_currentDeflateBlock->isValid() ) {
            /* A gzip stream should at least contain an end-of-stream block! */
            throw std::logic_error( "Call readHeader and readBlockHeader before readFooter!" );
        }

        if ( m_didReadHeader ) {
            m_crc32Calculator.verify( footer.crc32 );
        }
        break;
    }
    case FileType::ZLIB:
    {
        zlib::readFooter( m_bitReader );

        if ( !m_currentDeflateBlock.has_value() || !m_currentDeflateBlock->isValid() ) {
            /* A gzip stream should at least contain an end-of-stream block! */
            throw std::logic_error( "Call readHeader and readBlockHeader before readFooter!" );
        }

        break;
    }
    case FileType::DEFLATE:
    {
        const auto remainder = m_bitReader.tell() % BYTE_SIZE;
        if ( remainder != 0 ) {
            /* Try to read the remaining bits and the start of the next byte to determine whether we should
             * regard the remaining bits as end-of-file padding for the end of the raw deflate stream. */
            try {
                m_bitReader.peek<7>();
            } catch ( const gzip::BitReader::EndOfFileReached& ) {
                /* Skip the padding bits to get into correct EOF state. */
                m_bitReader.read( BYTE_SIZE - remainder );
            }
        }
        break;
    }
    case FileType::BZIP2:
        throw std::domain_error( "Bzip2 not supported by this class!" );
    }

    if ( m_bitReader.eof() ) {
        m_atEndOfFile = true;
    }

    m_currentPoint = StoppingPoint::END_OF_STREAM;
    m_didReadHeader = false;
}
}  // namespace rapidgzip
