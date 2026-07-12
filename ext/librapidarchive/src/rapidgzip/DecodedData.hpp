#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <core/FasterVector.hpp>
#include <core/VectorView.hpp>

#include "DecodedDataView.hpp"
#include "gzip/definitions.hpp"
#include "MarkerReplacement.hpp"


namespace rapidgzip::deflate
{
using MarkerVector = FasterVector<uint16_t>;
using DecodedVector = FasterVector<uint8_t>;


struct DecodedData
{
public:
    using WindowView = VectorView<uint8_t>;

    class Iterator
    {
    public:
        /**
         * This iterator provides a view of the decoded data as requested via an offset and a length.
         * If no relative offset or length is specified it will create a view of all of the data.
         * The interface will return subviews as pointer-length pairs because the data might not be
         * in one contiguous chunk internally.
         */
        explicit
        Iterator( const DecodedData& decodedData,
                  const size_t       offsetInChunk = 0,
                  const size_t       size = std::numeric_limits<size_t>::max() ) :
            m_data( decodedData ),
            m_size( size ),
            m_offsetInChunk( offsetInChunk )
        {
            /* Iterate over the chunks and decrease m_offsetInChunk by each chunk size until it validly points
             * into m_currentChunk, i.e., is smaller than its size. */
            for ( m_currentChunk = 0; m_currentChunk < m_data.data.size(); ++m_currentChunk ) {
                const auto& chunk = m_data.data[m_currentChunk];
                if ( ( m_offsetInChunk < chunk.size() ) && !chunk.empty() ) {
                    m_sizeInChunk = std::min( chunk.size() - m_offsetInChunk, m_size );
                    break;
                }
                m_offsetInChunk -= chunk.size();
            }
        }

        [[nodiscard]] operator bool() const noexcept
        {
            return ( m_currentChunk < m_data.data.size() ) && ( m_processedSize < m_size );
        }

        [[nodiscard]] std::pair<const void*, size_t>
        operator*() const
        {
            const auto& chunk = m_data.data[m_currentChunk];
            return { chunk.data() + m_offsetInChunk, m_sizeInChunk };
        }

        void
        operator++()
        {
            m_processedSize += m_sizeInChunk;
            m_offsetInChunk = 0;
            m_sizeInChunk = 0;

            if ( m_processedSize > m_size ) {
                throw std::logic_error( "Iterated over more bytes than was requested!" );
            }

            if ( !static_cast<bool>( *this ) ) {
                return;
            }

            ++m_currentChunk;
            for ( ; m_currentChunk < m_data.data.size(); ++m_currentChunk ) {
                const auto& chunk = m_data.data[m_currentChunk];
                if ( !chunk.empty() ) {
                    m_sizeInChunk = std::min( chunk.size(), m_size - m_processedSize );
                    break;
                }
            }
        }

    private:
        const DecodedData& m_data;
        const size_t m_size;

        size_t m_currentChunk{ 0 };
        size_t m_offsetInChunk{ 0 };
        size_t m_sizeInChunk{ 0 };
        size_t m_processedSize{ 0 };
    };

public:
    void
    append( DecodedVector&& toAppend )
    {
        if ( !toAppend.empty() ) {
            dataBuffers.emplace_back( std::move( toAppend ) );
            dataBuffers.back().shrink_to_fit();
            data.emplace_back( dataBuffers.back().data(), dataBuffers.back().size() );
        }
    }

    void
    append( DecodedDataView const& buffers );

    [[nodiscard]] size_t
    dataSize() const noexcept
    {
        const auto addSize = [] ( const size_t size, const auto& container ) { return size + container.size(); };
        return std::accumulate( data.begin(), data.end(), size_t( 0 ), addSize );
    }

    [[nodiscard]] size_t
    dataWithMarkersSize() const noexcept
    {
        const auto addSize = [] ( const size_t size, const auto& container ) { return size + container.size(); };
        return std::accumulate( dataWithMarkers.begin(), dataWithMarkers.end(), size_t( 0 ), addSize );
    }

    [[nodiscard]] size_t
    size() const noexcept
    {
        return dataSize() + dataWithMarkersSize();
    }

    [[nodiscard]] size_t
    sizeInBytes() const noexcept
    {
        return dataSize() * sizeof( uint8_t ) + dataWithMarkersSize() * sizeof( uint16_t );
    }

    /**
     * This is used to determine whether it is necessary to call applyWindow.
     * Testing for @ref dataWithMarkers.empty() is not sufficient because markers could be contained
     * in other members for derived classes! In that case @ref containsMarkers will be overridden.
     * @note Probably should not be called internally because it is allowed to be shadowed by a child class method.
     */
    [[nodiscard]] bool
    containsMarkers() const noexcept
    {
        return !dataWithMarkers.empty();
    }

    [[nodiscard]] size_t
    countMarkerSymbols() const;

    /**
     * Replaces all 16-bit wide marker symbols by looking up the referenced 8-bit symbols in @p window.
     * @note Probably should not be called internally because it is allowed to be shadowed by a child class method.
     */
    void
    applyWindow( WindowView const& window );

    /**
     * Returns the last 32 KiB decoded bytes. This can be called after decoding a block has finished
     * and then can be used to store and load it with deflate::Block::setInitialWindow to restart decoding
     * with the next block. Because this is not supposed to be called very often, it returns a copy of
     * the data instead of views.
     */
    [[nodiscard]] DecodedVector
    getLastWindow( WindowView const& previousWindow ) const;

    /**
     * @param skipBytes The number of bits to shift the previous window and fill it with new data.
     *        A value of 0 would simply return @p previousWindow while a value equal to size() would return
     *        the window as it would be after this whole block.
     * @note Should only be called after @ref applyWindow because @p skipBytes larger than @ref dataSize will throw.
     * @note Probably should not be called internally because it is allowed to be shadowed by a child class method.
     */
    [[nodiscard]] DecodedVector
    getWindowAt( WindowView const& previousWindow,
                 size_t            skipBytes ) const;

    void
    shrinkToFit()
    {
        for ( auto& container : dataBuffers ) {
            container.shrink_to_fit();
        }
        for ( auto& container : dataWithMarkers ) {
            container.shrink_to_fit();
        }
    }

    /**
     * Check decoded blocks that account for possible markers whether they actually contain markers and, if not so,
     * convert and move them to actual decoded data.
     */
    void
    cleanUnmarkedData();

#ifdef TEST_DECODED_DATA
    [[nodiscard]] const std::vector<MarkerVector>&
    getDataWithMarkers() const noexcept
    {
        return dataWithMarkers;
    }

    [[nodiscard]] const std::vector<VectorView<uint8_t> >&
    getData() const noexcept
    {
        return data;
    }
#endif

private:
    /**
     * Use vectors of vectors to avoid reallocations. The order of this data is:
     * - @ref dataWithMarkers (front to back)
     * - @ref data (front to back)
     * This order is fixed because there should be no reason for markers after we got enough data without markers!
     * There is no append( DecodedData ) method because this property might not be retained after using
     * @ref cleanUnmarkedData.
     */
    std::vector<MarkerVector> dataWithMarkers;
    std::vector<MarkerVector> reusedDataBuffers;
    std::vector<DecodedVector> dataBuffers;
    std::vector<VectorView<uint8_t> > data;
};


inline void
DecodedData::append( DecodedDataView const& buffers )
{
    static constexpr auto ALLOCATION_CHUNK_SIZE = 128_Ki;

    const auto& appendToEquallySizedChunks =
        [] ( auto&       targetChunks,
             const auto& buffer )
        {
            constexpr auto ALLOCATION_ELEMENT_COUNT = ALLOCATION_CHUNK_SIZE / sizeof( targetChunks[0][0] );

            if ( targetChunks.empty() ) {
                targetChunks.emplace_back().reserve( ALLOCATION_ELEMENT_COUNT );
            }

            for ( size_t nCopied = 0; nCopied < buffer.size(); ) {
                auto& copyTarget = targetChunks.back();
                const auto nFreeElements = copyTarget.capacity() - copyTarget.size();
                if ( nFreeElements == 0 ) {
                    targetChunks.emplace_back().reserve( ALLOCATION_ELEMENT_COUNT );
                    continue;
                }

                const auto nToCopy = std::min( nFreeElements, buffer.size() - nCopied );
                copyTarget.insert( copyTarget.end(), buffer.begin() + nCopied, buffer.begin() + nCopied + nToCopy );
                nCopied += nToCopy;
            }
        };

    if ( buffers.dataWithMarkersSize() > 0 ) {
        if ( !data.empty() ) {
            throw std::invalid_argument( "It is not allowed to append data with markers when fully decoded data "
                                         "has already been appended because the ordering will be wrong!" );
        }

        for ( const auto& buffer : buffers.dataWithMarkers ) {
            appendToEquallySizedChunks( dataWithMarkers, buffer );
        }
    }

    /* Add complexity to the already complex dataBuffers + data (views) structure by trying to force the
     * dataBuffer chunks to 128 KiB makes no sense because this method for appending views is only called
     * when decompressing with rapidgzip and as soon as we have 32 KiB of symbols, the decompression should
     * delegate to ISA-L except in pathological edge cases such as very large deflate blocks. */
    if ( buffers.dataSize() > 0 ) {
        auto& copied = dataBuffers.emplace_back();
        copied.reserve( buffers.dataSize() );
        for ( const auto& buffer : buffers.data ) {
            copied.insert( copied.end(), buffer.begin(), buffer.end() );
        }
        data.emplace_back( copied.data(), copied.size() );
    }
}


[[nodiscard]] inline size_t
DecodedData::countMarkerSymbols() const
{
    size_t result{ 0 };
    for ( const auto& chunk : dataWithMarkers ) {
        result += std::count_if( chunk.begin(), chunk.end(),
                                 [] ( const uint16_t symbol ) { return ( symbol & 0xFF00U ) != 0; } );
    }
    return result;
}


inline void
DecodedData::applyWindow( WindowView const& window )
{
    const auto markerCount = dataWithMarkersSize();
    if ( markerCount == 0 ) {
        dataWithMarkers.clear();
        return;
    }

    /* Because of the overhead of copying the window, avoid it for small replacements. */
    if ( markerCount >= 128_Ki ) {
        const std::array<uint8_t, 64_Ki> fullWindow =
            [&window] () noexcept
            {
                std::array<uint8_t, 64_Ki> result{};
                std::iota( result.begin(), result.begin() + 256, 0 );
                std::copy( window.begin(), window.end(), result.begin() + MAX_WINDOW_SIZE );
                return result;
            }();

        for ( auto& chunk : dataWithMarkers ) {
            auto* const target = reinterpret_cast<uint8_t*>( chunk.data() );
            /* Do not use std::transform because it allows out-of-order execution!
             * std::transform might be ~3% faster for FASTQ files btu it is hard to measure as timings seem
             * to vary a lot. Still, we need to be correct before being fast, but maybe something can be
             * done with inline Assembler or the intermediary "compressed" marker format might also help.
             * Note that these 3% shouldn't matter because we already are quite faster than before:
             *     rapidgzip-v0.9.0 : 2706.3 <= 2862.5 +- 2.3 <= 2979.9
             *     rapidgzip-v0.8.1 : 1988.2 <= 2067.8 +- 1.4 <= 2139.8
             */
            for ( size_t i = 0; i < chunk.size(); ++i ) {
                target[i] = fullWindow[chunk[i]];
            }
        }
    } else {
        /* For maximum size windows, we can skip one check because even UINT16_MAX is valid. */
        static_assert( std::numeric_limits<uint16_t>::max() - MAX_WINDOW_SIZE + 1U == MAX_WINDOW_SIZE );
        if ( window.size() >= MAX_WINDOW_SIZE ) {
            const MapMarkers<true> mapMarkers( window );
            for ( auto& chunk : dataWithMarkers ) {
                /* Do not use std::transform because it allows out-of-order execution and if a later i
                 * is computed first, it might overwrite values that are need for earlier i's because
                 * we are transforming in-place into a vector with smaller element size! */
                auto* const target = reinterpret_cast<uint8_t*>( chunk.data() );
                for ( size_t i = 0; i < chunk.size(); ++i ) {
                    target[i] = mapMarkers( chunk[i] );
                }
            }
        } else {
            const MapMarkers<false> mapMarkers( window );
            for ( auto& chunk : dataWithMarkers ) {
                auto* const target = reinterpret_cast<uint8_t*>( chunk.data() );
                /* Do not use std::transform because it allows out-of-order execution! */
                for ( size_t i = 0; i < chunk.size(); ++i ) {
                    target[i] = mapMarkers( chunk[i] );
                }
            }
        }
    }

    if ( !reusedDataBuffers.empty() ) {
        throw std::logic_error( "It seems like data already was replaced but we still got markers!" );
    }
    std::swap( reusedDataBuffers, dataWithMarkers );

    /* Prepend a VectorView to @ref data for each reused chunk buffer. */
    std::vector<VectorView<uint8_t> > dataViews;
    dataViews.reserve( reusedDataBuffers.size() + data.size() );
    for ( auto& chunk : reusedDataBuffers ) {
        /**
         * @todo Note that this leaves half of the chunk space unused because the number of elements stays
         *       the same while the element type size is halved. I assume that shrink_to_fit would be
         *       expensive. Therefore, it would be nice to simple join neihgboring chunks to fill all available
         *       space and free the rest of the chunks. But, depending on the individual chunk sizes,
         *       this can become complex.
         * @note It is kind of an exception that this works! Because reinterpret_cast with target types
         *       (unsigned) char* are excepted from the strict aliasing rules. For other types,
         *       it would be undefined behavior, e.g., trying to reuse a vector of uint32_t as uint16_t!
         * @see https://en.cppreference.com/w/cpp/language/reinterpret_cast#Type_aliasing
         */
        dataViews.emplace_back( reinterpret_cast<uint8_t*>( chunk.data() ), chunk.size() );
    }
    std::move( data.begin(), data.end(), std::back_inserter( dataViews ) );
    std::swap( data, dataViews );

    dataWithMarkers.clear();
}


[[nodiscard]] inline DecodedVector
DecodedData::getLastWindow( WindowView const& previousWindow ) const
{
    return getWindowAt( previousWindow, size() );
}


[[nodiscard]] inline DecodedVector
DecodedData::getWindowAt( WindowView const& previousWindow,
                          size_t const      skipBytes ) const
{
    if ( skipBytes > size() ) {
        throw std::invalid_argument( "Amount of bytes to skip is larger than this block!" );
    }

    DecodedVector window( MAX_WINDOW_SIZE );
    size_t prefilled{ 0 };
    if ( skipBytes < MAX_WINDOW_SIZE ) {
        const auto lastBytesToCopyFromPrevious = MAX_WINDOW_SIZE - skipBytes;
        if ( lastBytesToCopyFromPrevious <= previousWindow.size() ) {
            for ( size_t j = previousWindow.size() - lastBytesToCopyFromPrevious; j < previousWindow.size();
                  ++j, ++prefilled )
            {
                window[prefilled] = previousWindow[j];
            }
            // prefilled = lastBytesToCopyFromPrevious = MAX_WINDOW_SIZE - skipBytes
        } else {
            /* If previousWindow.size() < MAX_WINDOW_SIZE, which might happen at the start of streams,
             * then behave as if previousWindow was padded with leading zeros. */
            const auto zerosToFill = lastBytesToCopyFromPrevious - previousWindow.size();
            for ( ; prefilled < zerosToFill; ++prefilled ) {
                window[prefilled] = 0;
            }

            for ( size_t j = 0; j < previousWindow.size(); ++j, ++prefilled ) {
                window[prefilled] = previousWindow[j];
            }
            // prefilled = lastBytesToCopyFromPrevious - previousWindow.size() + previousWindow.size()
        }
        assert( prefilled == MAX_WINDOW_SIZE - skipBytes );
    }

    const auto remainingBytes = window.size() - prefilled;

    /* Skip over skipBytes in data and then copy the last remainingBytes before it. */
    auto offset = skipBytes - remainingBytes;
    /* if skipBytes < MAX_WINDOW_SIZE
     *     offset = skipBytes - ( window.size() - ( MAX_WINDOW_SIZE - skipBytes ) ) = 0
     * if skipBytes >= MAX_WINDOW_SIZE
     *     offset = skipBytes - ( window.size() - 0 ) */

    const auto copyFromDataWithMarkers =
        [this, &offset, &prefilled, &window] ( const auto& mapMarker )
        {
            for ( const auto& chunk : dataWithMarkers ) {
                if ( prefilled >= window.size() ) {
                    break;
                }

                if ( offset >= chunk.size() ) {
                    offset -= chunk.size();
                    continue;
                }

                for ( size_t i = offset; ( i < chunk.size() ) && ( prefilled < window.size() ); ++i, ++prefilled ) {
                    window[prefilled] = mapMarker( chunk[i] );
                }
                offset = 0;
            }
        };

    if ( previousWindow.size() >= MAX_WINDOW_SIZE ) {
        copyFromDataWithMarkers( MapMarkers</* full window */ true>( previousWindow ) );
    } else {
        copyFromDataWithMarkers( MapMarkers</* full window */ false>( previousWindow ) );
    }

    for ( const auto& chunk : data ) {
        if ( prefilled >= window.size() ) {
            break;
        }

        if ( offset >= chunk.size() ) {
            offset -= chunk.size();
            continue;
        }

        for ( size_t i = offset; ( i < chunk.size() ) && ( prefilled < window.size() ); ++i, ++prefilled ) {
            window[prefilled] = chunk[i];
        }
        offset = 0;
    }

    return window;
}


inline void
DecodedData::cleanUnmarkedData()
{
    while ( !dataWithMarkers.empty() ) {
        const auto& toDowncast = dataWithMarkers.back();
        /* Try to not only downcast whole chunks of data but also as many bytes as possible for the last chunk. */
        const auto marker = std::find_if(
            toDowncast.rbegin(), toDowncast.rend(),
            [] ( auto value ) { return value > std::numeric_limits<uint8_t>::max(); } );

        const auto sizeWithoutMarkers = static_cast<size_t>( std::distance( toDowncast.rbegin(), marker ) );
        auto downcasted = dataBuffers.emplace( dataBuffers.begin(), sizeWithoutMarkers );
        data.insert( data.begin(), VectorView<uint8_t>( downcasted->data(), downcasted->size() ) );
        std::transform( marker.base(), toDowncast.end(), downcasted->begin(),
                        [] ( auto symbol ) { return static_cast<uint8_t>( symbol ); } );

        if ( marker == toDowncast.rend() ) {
            dataWithMarkers.pop_back();
        } else {
            dataWithMarkers.back().resize( dataWithMarkers.back().size() - sizeWithoutMarkers );
            break;
        }
    }

    shrinkToFit();
}


/**
 * m rapidgzip && src/tools/rapidgzip -v -d -c -P 0 4GiB-base64.gz | wc -c
 * Non-polymorphic: Decompressed in total 4294967296 B in 1.49444 s -> 2873.96 MB/s
 * With virtual ~DecodedData() = default: Decompressed in total 4294967296 B in 3.58325 s -> 1198.62 MB/s
 * I don't know why it happens. Maybe it affects inline of function calls or moves of instances.
 */
static_assert( !std::is_polymorphic_v<DecodedData>, "Simply making it polymorphic halves performance!" );


#ifdef HAVE_IOVEC
[[nodiscard]] inline std::vector<::iovec>
toIoVec( const DecodedData& decodedData,
         const size_t       offsetInBlock,
         const size_t       dataToWriteSize )
{
    std::vector<::iovec> buffersToWrite;
    for ( auto it = rapidgzip::deflate::DecodedData::Iterator( decodedData, offsetInBlock, dataToWriteSize );
          static_cast<bool>( it ); ++it )
    {
        const auto& [data, size] = *it;
        ::iovec buffer{};
        /* The const_cast should be safe because vmsplice and writev should not modify the input data. */
        buffer.iov_base = const_cast<void*>( reinterpret_cast<const void*>( data ) );;
        buffer.iov_len = size;
        buffersToWrite.emplace_back( buffer );
    }
    return buffersToWrite;
}
#endif  // HAVE_IOVEC
}  // namespace rapidgzip::deflate
