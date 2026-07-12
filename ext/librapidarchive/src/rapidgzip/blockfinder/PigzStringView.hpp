#pragma once

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

#include <core/common.hpp>
#include <filereader/Buffered.hpp>
#include <filereader/FileReader.hpp>
#include <rapidgzip/gzip/definitions.hpp>
#include <rapidgzip/gzip/gzip.hpp>

#include "Interface.hpp"


/**
 * A simplistic pigz block finder reaches 1.3 GB/s
 * A parallel implementation of the naive pigz block finder reached 2.3 GB/s.
 * This pigz blockfinder makes use of std::string_view::find to reach 8 GB/s.
 */
namespace rapidgzip::blockfinder
{
class PigzStringView final :
    public Interface
{
public:
    /**
     * Should probably be larger than the I/O block size of 4096 B and smaller than most L1 cache sizes.
     * Not fitting into L1 cache isn't as bad as thought but increasing the size past 16 kiB also does not improve
     * the timings anymore on my Ryzen 3900X.
     */
    static constexpr size_t BUFFER_SIZE = 16_Ki;
    static constexpr uint8_t MAGIC_BIT_STRING_SIZE = 35;

public:
    explicit
    PigzStringView( UniqueFileReader fileReader ) :
        m_fileReader( std::move( fileReader ) ),
        m_fileSize( m_fileReader->size() )
    {}

    /**
     * @return offset of deflate block in bits (not the gzip stream offset!).
     */
    [[nodiscard]] size_t
    find() override
    {
        while ( m_blockOffsets.empty() && !m_fileReader->eof() && !m_fileReader->fail() && !m_fileReader->closed() ) {
            if ( foundFirstBlock ) {
                analyzeNextChunk();
            } else {
                findFirstBlock();
            }
        }

        if ( m_blockOffsets.empty() ) {
            return std::numeric_limits<std::size_t>::max();
        }

        m_lastReturnedBlockOffset = m_blockOffsets.back() * CHAR_BIT;
        m_blockOffsets.pop_back();
        return m_lastReturnedBlockOffset;
    }

private:
    void
    findBlockOffsets( const std::string_view& stringView,
                      std::size_t             offset )
    {
        for ( auto position = stringView.find( EMPTY_DEFLATE_BLOCK.data(), 0, EMPTY_DEFLATE_BLOCK.size() );
              position != std::string_view::npos;
              position = stringView.find( EMPTY_DEFLATE_BLOCK.data(), position + 1, EMPTY_DEFLATE_BLOCK.size() ) )
        {
            if ( ( position >= 1 )
                 /* Note that the additional check of the three-bits only works if the padding is filled with zeros. */
                 && ( ( static_cast<uint8_t>( stringView[position - 1] ) & 0b1110'0000 ) == 0 ) )
            {
                const auto totalOffset = offset + position + EMPTY_DEFLATE_BLOCK.size();
                auto fileSize = m_fileSize ? m_fileSize : m_fileReader->size();
                if ( !fileSize || ( totalOffset < *fileSize ) ) {
                    m_blockOffsets.push_back( totalOffset );
                }
            }
        }
    }

    void
    analyzeNextChunk()
    {
        constexpr std::size_t nBytesToRetain = ceilDiv( MAGIC_BIT_STRING_SIZE, CHAR_BIT ) - 1;
        static_assert( nBytesToRetain == 4, "Assuming bit string size of 35 for empty deflate block." );
        const auto checkBoundary = m_bufferSize > 0;
        /* We want to be able to find strings even if only their first byte is in the last buffer
         * or even if only their last byte is in the next buffer and of course all cases between. */
        std::array<char, 2 * nBytesToRetain> boundaryBuffer{};
        std::size_t boundaryBufferSize = 0;

        if ( checkBoundary ) {
            boundaryBufferSize = std::min( m_bufferSize, nBytesToRetain );
            for ( std::size_t i = 0; i < boundaryBufferSize; ++i ) {
                boundaryBuffer[i] = m_buffer[i + ( m_bufferSize - nBytesToRetain )];
            }
        }

        /* Always read chunks of BUFFER_SIZE in order to have aligned I/O and memory accesses. */
        const auto bufferOffset = m_fileReader->tell();
        const auto boundaryBufferOffset = bufferOffset - boundaryBufferSize;
        m_bufferSize = m_fileReader->read( m_buffer.data(), BUFFER_SIZE );

        if ( checkBoundary ) {
            boundaryBufferSize += std::min( nBytesToRetain, m_bufferSize );
            for ( std::size_t i = 0; i < std::min( nBytesToRetain, m_bufferSize ); ++i ) {
                boundaryBuffer[nBytesToRetain + i] = m_buffer[i];
            }

            findBlockOffsets( { boundaryBuffer.data(), boundaryBufferSize }, boundaryBufferOffset );
        }

        findBlockOffsets( { m_buffer.data(), m_bufferSize }, bufferOffset );
    }

    void
    findFirstBlock()
    {
        #if 0
        /**
         * @todo This requires the buffer to be larger than the first gzip header may be.
         * Theoretically, the user could store arbitrary amount of data in the zero-terminated file name
         * and file comment ...
         * @todo Make clone work here
         */
        gzip::BitReader bitReader( m_fileReader->clone() );

        #else

        BufferedFileReader::AlignedBuffer buffer( BUFFER_SIZE );
        buffer.resize( m_fileReader->read( buffer.data(), buffer.size() ) );
        gzip::BitReader bitReader( std::make_unique<BufferedFileReader>( std::move( buffer ) ) );

        #endif

        if ( ( rapidgzip::gzip::checkHeader( bitReader ) == rapidgzip::Error::NONE )
             && ( bitReader.tell() % CHAR_BIT == 0 ) ) {
            m_blockOffsets.push_back( bitReader.tell() / CHAR_BIT );
            /* Do not seek directly to one byte after the found offset in order to keep I/O aligned. */
            m_fileReader->seekTo( 0 );
            m_bufferSize = 0;
            foundFirstBlock = true;
            return;
        }

        /* If the first block couldn't be found, then don't even try to search for the others because the
         * first block would be missing, i.e., it would lead to incomplete data! */
        m_fileReader->seek( 0, SEEK_END );
    }

private:
    const UniqueFileReader m_fileReader;
    const std::optional<std::size_t> m_fileSize;

    alignas( 64 ) std::array<char, BUFFER_SIZE> m_buffer{};
    size_t m_bufferSize{ 0 };

    bool foundFirstBlock{ false };
    std::vector<std::size_t> m_blockOffsets;
    std::size_t m_lastReturnedBlockOffset{ 0 };

    static constexpr std::string_view EMPTY_DEFLATE_BLOCK{
        "\0\0\xFF\xFF", 4  /* required or else strlen is used resulting in zero */
    };
};
}  // rapidgzip::blockfinder
