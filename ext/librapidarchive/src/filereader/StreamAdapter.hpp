#pragma once

#include <ostream>
#include <vector>

#include "FileReader.hpp"


namespace rapidgzip
{
[[nodiscard]] int
toOrigin( std::ios_base::seekdir anchor )
{
    switch ( anchor )
    {
    case std::ios_base::beg: return SEEK_SET;
    case std::ios_base::cur: return SEEK_CUR;
    case std::ios_base::end: return SEEK_END;
    default: break;
    }
    return SEEK_SET;
}


/**
 * This class implements the std::streambuf interface by forwarding calls to the appropriate FileReader methods.
 *
 * It manages a buffer called "associated character sequence" buffering the access to the underlying FileReader,
 * called "controlled character sequence".
 * The interface is range-like and consists of a begin, end, and current pointer for the get and put area.
 * Because the FileReader is read-only, the put area as well as the overflow method is not supported.
 *
 * @todo Implement xsgetn, putbackfail for more efficient reading.
 * @todo Specialize for ParallelGzipReader after adding an API that returns the ChunkData pointer.
 *       This would avoid another unnecessary copy and the ChunkData would simply be the buffer.
 */
class FileReaderStreamBuffer :
    public std::streambuf
{
public:
    static constexpr size_t BUFFER_SIZE = 8ULL * 1024ULL;

    explicit
    FileReaderStreamBuffer( UniqueFileReader file ) :
        m_file( std::move( file ) )
    {
        if ( !m_file ) {
            throw std::invalid_argument( "May only be opened with a valid FileReader!" );
        }

        /* Already point to the buffer but signal that it does not yet contain data. */
        clearGetArea();
    }

    virtual
    ~FileReaderStreamBuffer()
    {
        sync();
    }

protected:
    /**
     * @return Number of characters that can be read from the internal buffer instead of the underlying file.
     * @note This implementation seems trivial to me. I don't understand why the default implementation only
     *       returns -1 instead of this generic implementation, probably because there is no interface for
     *       something like "is open".
     * @note "If showmanyc returns -1, underflow() and uflow() will definitely return Traits::eof or throw."
     *       So, it seems that -1 is a stronger 0, i.e., if not only the associated but also the controlled
     *       character sequence is at the end.
     * @note showmanyc stands for "Stream: HOW MANY Characters?".
     */
    std::streamsize
    showmanyc() override
    {
        if ( m_file->closed() ) {
            return -1;
        }
        return ( ( gptr() != nullptr ) && ( gptr() < egptr() ) ) ? std::streamsize( egptr() - gptr() ) : 0;
    }

    /**
     * Ensures that at least one character is available in the input area by updating the pointers
     * to the input area (if needed) and reading more data in from the input sequence (if applicable).
     *
     * @return The value of the character pointed to by the get pointer after the call on success,
     *         or Traits::eof() otherwise.
     */
    int_type
    underflow() override
    {
        if ( m_file->closed() ) {
            return traits_type::eof();
        }

        /* Still check even if:
         * > The public functions of std::streambuf call this function only
         * > if gptr() == nullptr or gptr() >= egptr(). */
        if ( ( gptr() != nullptr ) && ( gptr() < egptr() ) ) {
            return traits_type::to_int_type( *gptr() );
        }

        const auto nBytesRead = m_file->read( m_buffer.data(), m_buffer.size() );
        if ( nBytesRead == 0 ) {
            clearGetArea();
            return traits_type::eof();
        }

        setg( /* begin */ m_buffer.data(), /* current */ m_buffer.data(), /* end */ m_buffer.data() + nBytesRead );
        return traits_type::to_int_type( *gptr() );
    }

    /* Other overrides that are not required because this class is read_only.
     * putbackfail, sync: sufficient default implementations.
     * uflow: no override required because this class sets the get area pointers. */

    int_type
    overflow( int_type = traits_type::eof() ) override
    {
        throw std::runtime_error( "Writing is not supported!" );
    }

    pos_type
    seekoff( off_type                offset,
             std::ios_base::seekdir  anchor,
             std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out ) override
    {
        if ( ( mode & std::ios_base::out ) != 0 ) {
            throw std::runtime_error( "Writing is not supported!" );
        }

        /* We cannot simply forward to FileReader::tell or seek because the FileReader is positioned
         * at the end of the buffer! While this class' position somewhere inside the buffer. This means
         * that SEEK_CUR would be relative to some different offset and return wrong absolute offsets!
         * std::istream::tellg calls seekoff( 0, std::ios_base::cur, std::ios_base::in ) and should not
         * clear the buffer in general. */
        if ( anchor == std::ios_base::cur ) {
            /* Check if we can simply seek inside the buffer. */
            const auto bufferSize = static_cast<long long int>( egptr() - eback() );
            if ( bufferSize == 0 ) {
                return m_file->tell();
            }

            const auto bufferPosition = static_cast<long long int>( gptr() - eback() );
            const auto newPosition = bufferPosition + offset;
            /* It is allowed to seek at the end of the buffer. In that case, the next read call should simply
             * trigger underflow. However, it is not allowed to seek even further past the end of the buffer!
             * This <= bufferSize comparison instead of < also allows to handle tellg with an empty buffer,
             * e.g., before the first read, more gracefully. */
            if ( ( newPosition >= 0 ) && ( newPosition <= bufferSize ) ) {
                setg( /* begin */ eback(), /* current */ eback() + newPosition, /* end */ egptr() );
                const auto tellOffset = /* file offset corresponds to past-end of buffer offset */ m_file->tell()
                                        - ( egptr() - gptr() );
                return static_cast<pos_type>( tellOffset );
            }

            return seekpos( m_file->tell() - ( egptr() - gptr() ) + offset, mode );
        }

        /* Signal that an underflow is necessary. */
        clearGetArea();
        return static_cast<pos_type>( m_file->seek( offset, toOrigin( anchor ) ) );
    }

    /**
     * Same as @ref seekoff called with std::ios_base::beg.
     */
    pos_type
    seekpos( pos_type                offset,
             std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out ) override
    {
        if ( ( mode & std::ios_base::out ) != 0 ) {
            throw std::runtime_error( "Writing is not supported!" );
        }

        /* Signal that an underflow is necessary. */
        clearGetArea();
        return static_cast<pos_type>( m_file->seek( offset, SEEK_SET ) );
    }

private:
    void
    clearGetArea()
    {
        setg( /* begin */ m_buffer.data(), /* current */ m_buffer.data(), /* end */ m_buffer.data() );
    };

protected:
    const UniqueFileReader m_file;
    std::vector<char_type> m_buffer = std::vector<char_type>( BUFFER_SIZE );
};


class FileReaderStream :
    public FileReaderStreamBuffer,
    public std::istream
{
public:
    explicit
    FileReaderStream( UniqueFileReader file ) :
        FileReaderStreamBuffer( std::move( file ) ),
        std::istream( static_cast<std::streambuf*>( this ) )
    {}

    [[nodiscard]] bool
    is_open() const
    {
        return !m_file->closed();
    }

    /**
     * @note The destructor also closes, so calling this method should be rarely necessary.
     */
    void
    close() const
    {
        m_file->close();
    }
};
}  // namespace rapidgzip
