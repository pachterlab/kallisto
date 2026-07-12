#pragma once

#include <cstddef>
#include <cstdio>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include <core/common.hpp>


namespace rapidgzip
{
class FileReader;

using UniqueFileReader = std::unique_ptr<FileReader>;


/**
 * This file interface is heavily inspired by IOBase from Python because it is to be used from Python.
 * However, it strips anything related to file modification resulting a read-only file object.
 * @see https://docs.python.org/3/library/io.html#class-hierarchy
 */
class FileReader
{
public:
    FileReader() = default;

    virtual
    ~FileReader() = default;

    /* Delete copy constructors and assignments for performance and to avoid slicing. */

    FileReader( const FileReader& ) = delete;

    FileReader&
    operator=( const FileReader& ) = delete;

    /* Move constructors also are affected by slicing, but for performance, do not delete them to avoid
     * defaulted move constructores in derived classes not being created! We also do not want to manually
     * implement move constructors because it is error-prone. The only choice here is between bad and worse :(. */

    FileReader( FileReader&& ) = default;

    /* I almost never use move assignments, but move constructors occur during function calls. */
    FileReader&
    operator=( FileReader&& ) = delete;

    /**
     * @return A std::unique_ptr to a copy of this FileReader. The copied file reader should have the same
     *         file position if possible, i.e., if not closed.
     */
    [[nodiscard]] UniqueFileReader
    clone() const
    {
        auto fileReader = cloneRaw();
        if ( !fileReader->closed() && ( fileReader->tell() != tell() ) ) {
            fileReader->seekTo( tell() );
        }
        return fileReader;
    }

    virtual void
    close() = 0;

    [[nodiscard]] virtual bool
    closed() const = 0;

    [[nodiscard]] virtual bool
    eof() const = 0;

    [[nodiscard]] virtual bool
    fail() const = 0;

    [[nodiscard]] virtual int
    fileno() const = 0;

    [[nodiscard]] virtual bool
    seekable() const = 0;

    [[nodiscard]] virtual size_t
    read( char*  buffer,
          size_t nMaxBytesToRead ) = 0;

    virtual size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) = 0;

    size_t
    seekTo( uint64_t offset )
    {
        if ( offset > static_cast<uint64_t>( std::numeric_limits<long long int>::max() ) ) {
            throw std::invalid_argument( "Value " + std::to_string( offset ) + " out of range of long long int!" );
        }
        return seek( static_cast<long long int>( offset ) );
    }

    [[nodiscard]] virtual std::optional<size_t>
    size() const = 0;

    [[nodiscard]] virtual size_t
    tell() const = 0;

    virtual void
    clearerr() = 0;

protected:
    /**
     * Derived classes only need to construct a usable copy. Ideally, the state should be the same as the original.
     * The real @ref clone will call this and will ensure that the file position is transferred to the new object
     * if not already done.
     */
    [[nodiscard]] virtual UniqueFileReader
    cloneRaw() const
    {
        throw std::logic_error( "Not implemented!" );
    }

    [[nodiscard]] size_t
    effectiveOffset( long long int offset,
                     int           origin ) const
    {
        offset = [&] () {
            switch ( origin )
            {
            case SEEK_CUR:
                return saturatingAddition( static_cast<long long int>( tell() ), offset );
            case SEEK_SET:
                return offset;
            case SEEK_END:
                if ( const auto fileSize = size(); fileSize.has_value() ) {
                    return saturatingAddition( static_cast<long long int>( *fileSize ), offset );
                }
                throw std::logic_error( "File size is not available to seek from end!" );
            default:
                break;
            }
            throw std::invalid_argument( "Invalid seek origin supplied: " + std::to_string( origin ) );
        } ();

        const auto positiveOffset = static_cast<size_t>( std::max( offset, 0LL ) );
        /* Streaming file readers may return nullopt until EOF has been reached. */
        const auto fileSize = size();
        return fileSize.has_value() ? std::min( positiveOffset, *fileSize ) : positiveOffset;
    }
};
}  // namespace rapidgzip
