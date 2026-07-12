#pragma once

#include <algorithm>
#include <cstring>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <core/AlignedAllocator.hpp>
#include <core/common.hpp>
#include <core/VectorView.hpp>

#include "FileReader.hpp"


namespace rapidgzip
{
class BufferedFileReader :
    public FileReader
{
public:
    using AlignedBuffer = AlignedVector<char>;

public:
    explicit
    BufferedFileReader( UniqueFileReader fileReader,
                        size_t           bufferSize = 128_Ki ) :
        m_maxBufferSize( bufferSize ),
        m_file( std::move( fileReader ) )
    {}

    explicit
    BufferedFileReader( const std::vector<char>& inMemoryFileContents,
                        size_t                   bufferSize = 128_Ki ) :
        m_maxBufferSize( bufferSize ),
        m_buffer( inMemoryFileContents.begin(), inMemoryFileContents.end() )
    {}

    explicit
    BufferedFileReader( const VectorView<char>& inMemoryFileContents,
                        size_t                  bufferSize = 128_Ki ) :
        m_maxBufferSize( bufferSize ),
        m_buffer( inMemoryFileContents.begin(), inMemoryFileContents.end() )
    {}

    explicit
    BufferedFileReader( AlignedBuffer inMemoryFileContents,
                        size_t        bufferSize = 128_Ki ) :
        m_maxBufferSize( bufferSize ),
        m_buffer( std::move( inMemoryFileContents ) )
    {}

    [[nodiscard]] UniqueFileReader
    cloneRaw() const override
    {
        if ( m_file ) {
            return std::make_unique<BufferedFileReader>( m_file->clone() );
        }
        return std::make_unique<BufferedFileReader>( m_buffer, m_maxBufferSize );
    }

    /* Copying is simply not allowed because that might interfere with the file position state, use SharedFileReader! */

    void
    close() override
    {
        if ( m_file ) {
            m_file->close();
        }
        m_buffer.clear();
    }

    [[nodiscard]] bool
    closed() const override
    {
        return ( !m_file || m_file->closed() ) && m_buffer.empty();
    }

    [[nodiscard]] bool
    eof() const override
    {
        return ( !m_file || m_file->eof() ) && ( m_bufferPosition >= m_buffer.size() );
    }

    [[nodiscard]] bool
    fail() const override
    {
        return m_file && m_file->fail();
    }

    [[nodiscard]] int
    fileno() const override
    {
        if ( m_file ) {
            return m_file->fileno();
        }
        throw std::invalid_argument( "Trying to get fileno of an in-memory or closed file!" );
    }

    [[nodiscard]] bool
    seekable() const override
    {
        return !m_file || m_file->seekable();
    }

    [[nodiscard]] size_t
    read( char*  buffer,
          size_t nMaxBytesToRead ) override
    {
        if ( closed() ) {
            throw std::invalid_argument( "Cannot read from closed file!" );
        }

        if ( nMaxBytesToRead == 0 ) {
            return 0;
        }

        /* Read from buffer as much as possible. */
        auto nBytesRead = readFromBuffer( buffer, nMaxBytesToRead );
        if ( nBytesRead >= nMaxBytesToRead ) {
            return nBytesRead;
        }

        /* If we cannot refill the buffer, then that was it. */
        if ( !m_file ) {
            return nBytesRead;
        }

        /* Skip buffering step for very large reads. */
        if ( nMaxBytesToRead - nBytesRead >= m_maxBufferSize ) {
            const auto nBytesReadFromFile = m_file->read( buffer == nullptr ? buffer : buffer + nBytesRead,
                                                          nMaxBytesToRead - nBytesRead );
            m_originalBufferOffset += nBytesReadFromFile;
            return nBytesRead + nBytesReadFromFile;
        }

        refillBuffer();
        return nBytesRead + readFromBuffer( buffer + nBytesRead, nMaxBytesToRead - nBytesRead );
    }

    size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) override
    {
        if ( closed() ) {
            throw std::invalid_argument( "Cannot seek closed file!" );
        }

        const auto newBufferPosition = effectiveOffset( offset, origin );

        /* Check if we can simply seek inside the buffer. */
        if ( newBufferPosition <= m_buffer.size() ) {
            m_bufferPosition = newBufferPosition;
            return tell();
        }

        if ( !m_file ) {
            throw std::invalid_argument( "Cannot seek outside of in-memory file range!" );
        }

        m_originalBufferOffset = m_file->seek( static_cast<long long int>( tell() ) + offset, SEEK_SET );
        /* Clear buffer AFTER calling tell(), or else the position will be different! */
        m_bufferPosition = 0;
        m_buffer.clear();

        return tell();
    }

    [[nodiscard]] std::optional<size_t>
    size() const override
    {
        return m_file ? m_file->size() : m_buffer.size();
    }

    [[nodiscard]] size_t
    tell() const override
    {
        return m_originalBufferOffset + m_bufferPosition;
    }

    void
    clearerr() override
    {
        if ( m_file ) {
            m_file->clearerr();
        }
    }

private:
    void
    refillBuffer()
    {
        if ( !m_file ) {
            throw std::invalid_argument( "Cannot refill buffer for buffer-only file!" );
        }

        m_bufferPosition = 0;
        m_originalBufferOffset = m_file->seek( m_originalBufferOffset + m_buffer.size() );
        m_buffer.resize( m_maxBufferSize );
        const auto nBytesRead = m_file->read( m_buffer.data(), m_buffer.size() );
        m_buffer.resize( nBytesRead );
    }

    size_t
    readFromBuffer( char*  buffer,
                    size_t nMaxBytesToRead )
    {
        if ( m_bufferPosition >= m_buffer.size() ) {
            return 0;
        }

        const auto nBytesToReadFromBuffer = std::min( m_buffer.size() - m_bufferPosition, nMaxBytesToRead );
        if ( buffer != nullptr ) {
            std::memcpy( buffer, m_buffer.data() + m_bufferPosition, nBytesToReadFromBuffer );
        }
        m_bufferPosition += nBytesToReadFromBuffer;
        return nBytesToReadFromBuffer;
    }

protected:
    const size_t m_maxBufferSize;
    UniqueFileReader m_file;

    size_t m_originalBufferOffset{ 0 };
    AlignedBuffer m_buffer;
    size_t m_bufferPosition{ 0 };
};
}  // namespace rapidgzip
