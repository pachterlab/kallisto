#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <map>

#include <core/FileUtils.hpp>                    // writeAll
#include <filereader/FileReader.hpp>

using std::size_t;


namespace indexed_bzip2
{
class BZ2ReaderInterface :
    public rapidgzip::FileReader
{
public:
    using WriteFunctor = std::function<void ( const void*, uint64_t )>;

public:
    ~BZ2ReaderInterface() override = default;

    [[nodiscard]] size_t
    read( char*  outputBuffer,
          size_t nBytesToRead ) final
    {
        return read( -1, outputBuffer, nBytesToRead );
    }

    /**
     * @param[out] outputBuffer should at least be large enough to hold @p nBytesToRead bytes
     * @return number of bytes written
     */
    size_t
    read( const int    outputFileDescriptor = -1,
          char* const  outputBuffer         = nullptr,
          const size_t nBytesToRead         = std::numeric_limits<size_t>::max() )
    {
        const auto writeFunctor =
            [nBytesDecoded = uint64_t( 0 ), outputFileDescriptor, outputBuffer]
            ( const void* const buffer,
              uint64_t    const size ) mutable
            {
                auto* const currentBufferPosition = outputBuffer == nullptr ? nullptr : outputBuffer + nBytesDecoded;
                const auto errorCode = rapidgzip::writeAll( outputFileDescriptor, currentBufferPosition, buffer, size );
                if ( errorCode != 0 ) {
                    std::stringstream message;
                    message << "Failed to write all bytes because of: " << strerror( errorCode )
                            << " (" << errorCode << ")";
                    throw std::runtime_error( std::move( message ).str() );
                }
                nBytesDecoded += size;
            };

        return read( writeFunctor, nBytesToRead );
    }

    virtual size_t
    read( const WriteFunctor& writeFunctor,
          size_t              nBytesToRead = std::numeric_limits<size_t>::max() ) = 0;

    /* BZip2 specific methods */

    [[nodiscard]] virtual bool
    blockOffsetsComplete() const = 0;

    /**
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    virtual std::map<size_t, size_t>
    blockOffsets() = 0;

    /**
     * Same as @ref blockOffsets but it won't force calculation of all blocks and simply returns
     * what is available at call time.
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    [[nodiscard]] virtual std::map<size_t, size_t>
    availableBlockOffsets() const = 0;

    virtual void
    setBlockOffsets( std::map<size_t, size_t> offsets ) = 0;

    /**
     * @return number of processed bits of compressed bzip2 input file stream
     * @note Bzip2 is block based and blocks are currently read fully, meaning that the granularity
     *       of the returned position is ~100-900kB. It's only useful for a rough estimate.
     */
    [[nodiscard]] virtual size_t
    tellCompressed() const = 0;

    void
    setShowProfileOnDestruction( bool showProfileOnDestruction )
    {
        m_showProfileOnDestruction = showProfileOnDestruction;
    }

protected:
    bool m_showProfileOnDestruction{ false };
};
}  // namespace indexed_bzip2
