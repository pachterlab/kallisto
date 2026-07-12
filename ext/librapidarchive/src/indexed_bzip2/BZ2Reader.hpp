#pragma once

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <core/common.hpp>
#include <filereader/BitReader.hpp>
#include <filereader/FileReader.hpp>

#include "bzip2.hpp"
#include "BZ2ReaderInterface.hpp"

#ifdef WITH_PYTHON_SUPPORT
    #include <filereader/Python.hpp>
#endif


namespace indexed_bzip2
{
class BZ2Reader final :
    public BZ2ReaderInterface
{
public:
    using BlockHeader = bzip2::Block;
    using BitReader = bzip2::BitReader;

public:
    static constexpr size_t IOBUF_SIZE = 4096;

public:
    explicit
    BZ2Reader( rapidgzip::UniqueFileReader fileReader ) :
        /* Without this, the test_joining_archive test inside ratarmountcore's test_factory.py fails.
         * This is because tell() throws an exception because m_file.tell() is in an unexpected position.
         * I'm still not sure why but maybe some external caller resets the file position to 0.
         * As we are returning fileno to the underlying file and because the given PyObject input file object
         * can have shared ownership, we need to ensure that this class works even when the file object is
         * seeked from outside. */
        m_bitReader( ensureSharedFileReader( std::move( fileReader ) ) )
    {}

#ifdef WITH_PYTHON_SUPPORT
    explicit
    BZ2Reader( const std::string& filePath ) :
        BZ2Reader( std::make_unique<rapidgzip::StandardFileReader>( filePath ) )
    {}

    explicit
    BZ2Reader( int fileDescriptor ) :
        BZ2Reader( std::make_unique<rapidgzip::StandardFileReader>( fileDescriptor ) )
    {}

    explicit
    BZ2Reader( PyObject* pythonObject ) :
        BZ2Reader( std::make_unique<rapidgzip::PythonFileReader>( pythonObject ) )
    {}
#endif

    /* Forbid copying because it is hard to get right and there is not much use for it right now. */
    BZ2Reader( const BZ2Reader& ) = delete;

    BZ2Reader&
    operator=( const BZ2Reader& ) = delete;

    /* Forbid moving because it is not used right now but it could probably be defaulted as long
     * as BitReader has a working move constructor. */
    BZ2Reader( BZ2Reader&& ) = delete;

    BZ2Reader&
    operator=( BZ2Reader&& ) = delete;

    ~BZ2Reader() override
    {
        if ( m_showProfileOnDestruction ) {
            const auto& durations = m_statistics.durations;
            std::cerr << "[BZ2Reader] Time spent:\n";
            std::cerr << "    decodeBlock                   : " << durations.decodeBlock               << "s\n";
            std::cerr << "    readBlockHeader               : " << durations.readBlockHeader           << "s\n";
            std::cerr << "        readSymbolMaps            : " << durations.readSymbolMaps            << "s\n";
            std::cerr << "        readSelectors             : " << durations.readSelectors             << "s\n";
            std::cerr << "        readTrees                 : " << durations.readTrees                 << "s\n";
            std::cerr << "        createHuffmanTable        : " << durations.createHuffmanTable        << "s\n";
            std::cerr << "        burrowsWheelerPreparation : " << durations.burrowsWheelerPreparation << "s\n";
            std::cerr << std::endl;
        }
    }

    /* FileReader overrides */

    [[nodiscard]] int
    fileno() const override
    {
        throw std::logic_error( "This is a virtual file object, which has no corresponding file descriptor!" );
    }

    [[nodiscard]] bool
    seekable() const override
    {
        return m_bitReader.seekable();
    }

    void
    close() override
    {
        m_bitReader.close();
    }

    [[nodiscard]] bool
    closed() const override
    {
        return m_bitReader.closed();
    }

    [[nodiscard]] bool
    eof() const override
    {
        return m_atEndOfFile;
    }

    [[nodiscard]] bool
    fail() const override
    {
        throw std::logic_error( "Not implemented!" );
    }

    [[nodiscard]] size_t
    tell() const override
    {
        if ( m_atEndOfFile ) {
            const auto fileSize = size();
            if ( !fileSize ) {
                throw std::logic_error( "When the file end has been reached, the block map should have been finalized "
                                        "and the file size should be available!" );
            }
            return *fileSize;
        }
        return m_currentPosition;
    }

    [[nodiscard]] std::optional<size_t>
    size() const override
    {
        if ( !m_blockToDataOffsetsComplete ) {
            return std::nullopt;
        }
        return m_blockToDataOffsets.rbegin()->second;
    }

    size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) override;

    void
    clearerr() override
    {
        m_bitReader.clearerr();
        m_atEndOfFile = false;
        throw std::invalid_argument( "Not fully tested!" );
    }

    /* BZip2 specific methods */

    [[nodiscard]] uint32_t
    crc() const
    {
        return m_calculatedStreamCRC;
    }

    [[nodiscard]] bool
    blockOffsetsComplete() const override
    {
        return m_blockToDataOffsetsComplete;
    }

    /**
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    [[nodiscard]] std::map<size_t, size_t>
    blockOffsets() override
    {
        if ( !m_blockToDataOffsetsComplete ) {
            read();
        }

        return m_blockToDataOffsets;
    }

    /**
     * Same as @ref blockOffsets but it won't force calculation of all blocks and simply returns
     * what is available at call time.
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    [[nodiscard]] std::map<size_t, size_t>
    availableBlockOffsets() const override
    {
        return m_blockToDataOffsets;
    }

    void
    setBlockOffsets( std::map<size_t, size_t> offsets ) override
    {
        if ( offsets.size() < 2 ) {
            throw std::invalid_argument( "Block offset map must contain at least one valid block and one EOS block!" );
        }
        m_blockToDataOffsetsComplete = true;
        m_blockToDataOffsets = std::move( offsets );
    }

    /**
     * @return number of processed bits of compressed bzip2 input file stream
     * @note Bzip2 is block based and blocks are currently read fully, meaning that the granularity
     *       of the returned position is ~100-900kB. It's only useful for a rough estimate.
     */
    [[nodiscard]] size_t
    tellCompressed() const override
    {
        return m_bitReader.tell();
    }

    using BZ2ReaderInterface::read;

    size_t
    read( const WriteFunctor& writeFunctor,
          const size_t        nBytesToRead = std::numeric_limits<size_t>::max() ) override
    {
        size_t nBytesDecoded = 0;
        while ( ( nBytesDecoded < nBytesToRead ) && !m_bitReader.eof() && !eof() ) {
            /* The input may be a concatenation of multiple BZip2 files (like produced by pbzip2).
             * Therefore, iterate over those multiple files and decode them to the specified output. */
            if ( m_bitReader.tell() == 0 ) {
                readBzip2Header();
            } else if ( m_lastHeader.eos() ) {
                try {
                    readBzip2Header();
                } catch ( const std::domain_error& ) {
                    /* @TODO MIGHT THIS lead to a bug if the bzip2 file ends perfectly on a byte boundary
                     * such that m_bitReader.eof() will be true before this code part has been reached?! */
                    std::cerr << "[Warning] Trailing garbage after EOF ignored!\n";
                    m_atEndOfFile = true;
                    m_blockToDataOffsetsComplete = true;
                    break;
                }
            }

            nBytesDecoded += decodeStream( writeFunctor, nBytesToRead - nBytesDecoded );

        #ifdef WITH_PYTHON_SUPPORT
            rapidgzip::checkPythonSignalHandlers();
        #endif
        }
        m_currentPosition += nBytesDecoded;
        return nBytesDecoded;
    }

private:
    /**
     * @param  outputBuffer A char* to which the data is written.
     *                      You should ensure that at least @p maxBytesToFlush bytes can fit there!
     * @return The number of actually flushed bytes, which might be hindered,
     *         e.g., if the output file descriptor can't be written to!
     */
    size_t
    flushOutputBuffer( const WriteFunctor& writeFunctor = {},
                       size_t              maxBytesToFlush = std::numeric_limits<size_t>::max() );

    /**
     * Undo burrows-wheeler transform on intermediate buffer @ref dbuf to @ref outBuf
     *
     * Burrows-wheeler transform is described at:
     * @see http://dogma.net/markn/articles/bwt/bwt.htm
     * @see http://marknelson.us/1996/09/01/bwt/
     *
     * @return number of actually decoded bytes
     */
    [[nodiscard]] size_t
    decodeStream( const WriteFunctor& writeFunctor,
                  size_t              nMaxBytesToDecode = std::numeric_limits<size_t>::max() );

    BlockHeader
    readBlockHeader( size_t offsetBits );

protected:
    void
    readBzip2Header()
    {
        m_blockSize100k = bzip2::readBzip2Header( m_bitReader );
        m_calculatedStreamCRC = 0;
    }

protected:
    BitReader m_bitReader;

    uint8_t m_blockSize100k = 0;
    uint32_t m_streamCRC = 0;  /** CRC of stream as last block says */
    uint32_t m_calculatedStreamCRC = 0;
    bool m_blockToDataOffsetsComplete = false;
    size_t m_currentPosition = 0;  /** the current position as can only be modified with read or seek calls. */
    bool m_atEndOfFile = false;

    std::map<size_t, size_t> m_blockToDataOffsets;

private:
    BlockHeader m_lastHeader;

    /* This buffer is needed for decoding because decoding of runtime length encoded strings might lead to more
     * output data for the copies of a character than we can write out. In order to not save the outstanding copies
     * and the character to copy, this buffer acts as a kind of generalized "current decoder state". */
    std::vector<char> m_decodedBuffer = std::vector<char>( IOBUF_SIZE );
    /* it's strictly increasing during decoding and no previous data in m_decodedBuffer is acced,
     * so we can almost at any position clear m_decodedBuffer and set m_decodedBufferPos to 0, which is done for flushing! */
    size_t m_decodedBufferPos = 0;

    /** The sum over all decodeBuffer calls. This is used to create the block offset map */
    size_t m_decodedBytesCount = 0;

    BlockHeader::Statistics m_statistics;
};


inline size_t
BZ2Reader::seek( long long int offset,
                 int           origin )
{
    if ( origin == SEEK_END ) {
        /* size() requires the block offsets to be available! */
        if ( !m_blockToDataOffsetsComplete ) {
            read();
        }
    }

    const auto positiveOffset = effectiveOffset( offset, origin );

    if ( positiveOffset == tell() ) {
        return positiveOffset;
    }

    /* When block offsets are not complete yet, emulate forward seeking with a read. */
    if ( !m_blockToDataOffsetsComplete && ( positiveOffset > tell() ) ) {
        read( -1, nullptr, positiveOffset - tell() );
        return tell();
    }

    /* size() and then seeking requires the block offsets to be available! */
    if ( !m_blockToDataOffsetsComplete ) {
        read();
    }

    m_currentPosition = positiveOffset;

    flushOutputBuffer();  // ensure that no old data is left over

    m_atEndOfFile = positiveOffset >= size();
    if ( m_atEndOfFile ) {
        return tell();
    }

    /* find offset from map (key and values are sorted, so we can bisect!) */
    const auto blockOffset = std::lower_bound(
        m_blockToDataOffsets.rbegin(), m_blockToDataOffsets.rend(), std::make_pair( 0, positiveOffset ),
        [] ( std::pair<size_t, size_t> a, std::pair<size_t, size_t> b ) { return a.second > b.second; } );

    if ( ( blockOffset == m_blockToDataOffsets.rend() ) || ( positiveOffset < blockOffset->second ) ) {
        throw std::runtime_error( "Could not find block to seek to for given offset" );
    }
    const auto nBytesSeekInBlock = positiveOffset - blockOffset->second;

    m_statistics.merge( m_lastHeader.statistics );
    m_lastHeader = readBlockHeader( blockOffset->first );
    m_lastHeader.readBlockData();
    /* no decodeBzip2 necessary because we only seek inside one block! */
    const auto nBytesDecoded = decodeStream( {}, nBytesSeekInBlock );

    if ( nBytesDecoded != nBytesSeekInBlock ) {
        std::stringstream msg;
        msg << "Could not read the required " << nBytesSeekInBlock
            << " to seek in block but only " << nBytesDecoded << "\n";
        throw std::runtime_error( std::move( msg ).str() );
    }

    return m_currentPosition;
}


inline BZ2Reader::BlockHeader
BZ2Reader::readBlockHeader( size_t offsetBits )
{
    /* note that blocks are NOT byte-aligned! Only the end of the stream has a necessary padding. */
    if ( !m_blockToDataOffsetsComplete ) {
        m_blockToDataOffsets.insert( { offsetBits, m_decodedBytesCount } );
    }

    m_bitReader.seekTo( offsetBits );
    BlockHeader header( m_bitReader );

    if ( header.eos() ) {
        /* EOS block contains CRC for whole stream */
        m_streamCRC = header.headerCRC();

        if ( !m_blockToDataOffsetsComplete && ( m_streamCRC != m_calculatedStreamCRC ) ) {
            std::stringstream msg;
            msg << "[BZip2 block header] Stream CRC 0x" << std::hex << m_streamCRC
                << " does not match calculated CRC 0x" << m_calculatedStreamCRC;
            throw std::runtime_error( std::move( msg ).str() );
        }
    }

    m_atEndOfFile = header.eof();
    if ( header.eof() ) {
        m_blockToDataOffsetsComplete = true;
    }

    return header;
}


inline size_t
BZ2Reader::flushOutputBuffer( WriteFunctor const& writeFunctor,
                              size_t       const  maxBytesToFlush )
{
    const auto nBytesToFlush = std::min( m_decodedBufferPos, maxBytesToFlush );

    if ( writeFunctor ) {
        writeFunctor( m_decodedBuffer.data(), nBytesToFlush );
    }

    if ( nBytesToFlush > 0 ) {
        m_decodedBytesCount += nBytesToFlush;
        m_decodedBufferPos  -= nBytesToFlush;
        std::memmove( m_decodedBuffer.data(), m_decodedBuffer.data() + nBytesToFlush, m_decodedBufferPos );
    }

    return nBytesToFlush;
}


inline size_t
BZ2Reader::decodeStream( WriteFunctor const& writeFunctor,
                         size_t       const  nMaxBytesToDecode )
{
    if ( eof() || ( nMaxBytesToDecode == 0 ) ) {
        return 0;
    }

    /* try to flush remnants in output buffer from interrupted last call */
    size_t nBytesDecoded = flushOutputBuffer( writeFunctor, nMaxBytesToDecode );

    while ( nBytesDecoded < nMaxBytesToDecode ) {
        /* If we need to refill dbuf, do it. Only won't be required for resuming interrupted decodations. */
        if ( m_lastHeader.eob() ) {
            m_statistics.merge( m_lastHeader.statistics );
            m_lastHeader = readBlockHeader( m_bitReader.tell() );
            if ( m_lastHeader.eos() ) {
                return nBytesDecoded;
            }
            m_lastHeader.readBlockData();
        }

        /* m_decodedBufferPos should either be cleared by flush after or by flush before while!
         * It might happen that this is not the case when, e.g., the output file descriptor can't be written to.
         * However, if this happens, nBytesDecoded is very likely to not grow anymore and thereby we have to
         * throw to exit the infinite loop. */
        if ( m_decodedBufferPos > 0 ) {
            throw std::runtime_error( "[BZ2Reader::decodeStream] Could not write any of the decoded bytes to the "
                                      "file descriptor or buffer!" );
        }

        const auto nBytesToDecode = std::min( m_decodedBuffer.size(), nMaxBytesToDecode - nBytesDecoded );
        m_decodedBufferPos = m_lastHeader.read( nBytesToDecode, m_decodedBuffer.data() );

        if ( m_lastHeader.eob() && !m_blockToDataOffsetsComplete ) {
            m_calculatedStreamCRC = ( ( m_calculatedStreamCRC << 1U ) | ( m_calculatedStreamCRC >> 31U ) )
                                    ^ m_lastHeader.dataCRC();
        }

        /* required for correct data offsets in readBlockHeader and for while condition of course */
        nBytesDecoded += flushOutputBuffer( writeFunctor, nMaxBytesToDecode - nBytesDecoded );
    }

    return nBytesDecoded;
}
}  // namespace indexed_bzip2
