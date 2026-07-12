#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <core/FasterVector.hpp>
#include <rapidgzip/CompressedVector.hpp>
#include <rapidgzip/DecodedData.hpp>
#include <rapidgzip/gzip/crc32.hpp>
#include <rapidgzip/gzip/gzip.hpp>


namespace rapidgzip
{
/**
 * Rpmalloc does worse than standard malloc (Clang 13) for the case when using 128 cores, chunk size 4 MiB
 * with imported index of Silesia (compression ratio ~3.1), i.e., the decompressed chunk sizes are ~12 MiB
 * and probably deviate wildly in size (4-100 MiB maybe?). I guess that this leads to overallocation and
 * memory slab reuse issues in rpmalloc.
 * Allocating memory chunks in much more deterministic sizes seems to alleviate this problem immensely!
 *
 * Problematic commit:
 *     dd678c7 2022-11-06 mxmlnkn [performance] Use rpmalloc to increase multi-threaded malloc performance
 *
 * Approximate benchmark results for:
 *     rapidgzip -P $( nproc ) -o /dev/null -f --export-index index silesia-256x.tar.pigz
 *     taskset --cpu-list 0-127 rapidgzip -P 128 -d -o /dev/null --import-index index silesia-256x.tar.pigz
 *
 * Commit:
 *     rapidgzip-v0.5.0   16 GB/s
 *     dd678c7~1        16 GB/s
 *     dd678c7           8 GB/s
 *
 * dd678c7 with ALLOCATION_CHUNK_SIZE:
 *     64  KiB       19.4 19.7       GB/s
 *     256 KiB       20.8 20.7       GB/s
 *     1   MiB       21.2 20.7 20.8  GB/s
 *     4   MiB       8.4 8.5 8.3     GB/s
 *
 * It seems to be pretty stable across magnitudes as long as the number of allocations doesn't get too
 * large and as long as the allocation chunk size is much smaller than the decompressed data chunk size.
 * 1 MiB seems like the natural choice because the optimum (compressed) chunk size is around 4 MiB
 * and it would also be exactly one hugepage if support for that would ever be added.
 * Beware that each chunk is only as large as one gzip stream. And bgzip creates gzip streams that are only
 * ~64 KiB each! Therefore, when decoding bgzip while importing the index, we need to account for this here
 * and avoid frequent overallocations and resizes, which slow down the program immensely!
 *
 * Test with:
 *     m rapidgzip && src/tools/rapidgzip -o /dev/null -d \
 *       --import-index silesia-32x.tar.bgzip.gzi silesia-32x.tar.bgzip
 *
 * ALLOCATION_CHUNK_SIZE  Bandwidth
 *        32  KiB         2.4 GB/s
 *        64  KiB         2.5 GB/s
 *        128 KiB         2.5 GB/s
 *        256 KiB         2.4 GB/s
 *        512 KiB         1.5 GB/s
 *        1   MiB         370 MB/s
 */
static constexpr size_t ALLOCATION_CHUNK_SIZE = 128_Ki;


/**
 * This class adds higher-level capabilities onto @ref deflate::DecodedData, which was only intended for
 * returning decompression results and aggregating them during decompression of a single deflate block.
 * This class instead is intended to aggregate results from multiple deflate blocks, possibly even multiple
 * gzip streams. It is used to hold the chunk data for parallel decompression.
 * It also adds some further metadata like deflate block and stream boundaries and helpers for creating
 * evenly distributed checkpoints for a gzip seek index.
 *
 * Specialized use cases can optimize memory usage or add post-processing steps by implementing the two
 * @ref append methods, @ref applyWindow, and @ref finalize. The shadowed methods in the base class
 * should be called from the reimplemented methods in order to keep default functionality. This call
 * can also be knowingly omitted, e.g., for only counting bytes instead of appending them.
 *
 * - @ref append is called by @ref GzipChunkFetcher after each deflate::Block call back, which could be
 *   every block or up to maximum 32 KiB of decompressed data.
 * - @ref finalize is called after the first stage of decompression has finished.
 *   At this point, the number of elements in the chunk is finalized. Elements can be 16-bit wide markers.
 * - @ref applyWindow is called during the second decompression stage and the ChunkData will hold the fully
 *   decompressed data after this call.
 */
struct ChunkData :
    public deflate::DecodedData
{
    using BaseType = deflate::DecodedData;

    using Window = CompressedVector<FasterVector<uint8_t> >;
    using SharedWindow = std::shared_ptr<const Window>;
    using WindowView = VectorView<uint8_t>;

    struct Configuration
    {
        size_t splitChunkSize{ std::numeric_limits<size_t>::max() };
        /** This should be used to decide what kind of footer to expect and what to do after the footer. */
        FileType fileType{ FileType::NONE };
        bool crc32Enabled{ true };
        std::optional<CompressionType> windowCompressionType;
        bool windowSparsity{ true };
        /**
         * This is used by the chunk decoding implementations, but it feels more correct to have this stored here,
         * because it affects the chunk configuration and in the future it might be cleaner to check for a maximum
         * size inside ChunkData::append instead, or offer a getter, which could be overwritten / replaced with
         * dynamic checks for stopping decompression preemptively, e.g., if the decompressed data is known to be wrong.
         */
        size_t maxDecompressedChunkSize{ std::numeric_limits<size_t>::max() };
        std::optional<char> newlineCharacter{};
    };

    struct Subchunk
    {
        [[nodiscard]] bool
        operator==( const Subchunk& other ) const
        {
            return ( encodedOffset == other.encodedOffset )
                   && ( decodedOffset == other.decodedOffset )
                   && ( encodedSize == other.encodedSize )
                   && ( decodedSize == other.decodedSize )
                   && ( newlineCount == other.newlineCount )
                   && ( static_cast<bool>( window ) == static_cast<bool>( other.window ) )
                   && ( !static_cast<bool>( window ) || !static_cast<bool>( other.window )
                        || ( *window == *other.window ) );
        }

        [[nodiscard]] bool
        hasBeenPostProcessed( const bool requireNewlineCount ) const
        {
            return static_cast<bool>( window ) && usedWindowSymbols.empty()
                   && ( newlineCount.has_value() || !requireNewlineCount );
        }

    public:
        size_t encodedOffset{ 0 };
        size_t decodedOffset{ 0 };
        size_t encodedSize{ 0 };
        size_t decodedSize{ 0 };
        std::optional<size_t> newlineCount{};
        SharedWindow window{};
        std::vector<bool> usedWindowSymbols{};
    };

    class Statistics
    {
    public:
        void
        merge( const Statistics& other )
        {
            falsePositiveCount           += other.falsePositiveCount;
            blockFinderDuration          += other.blockFinderDuration;
            decodeDuration               += other.decodeDuration;
            decodeDurationInflateWrapper += other.decodeDurationInflateWrapper;
            decodeDurationIsal           += other.decodeDurationIsal;
            appendDuration               += other.appendDuration;
            applyWindowDuration          += other.applyWindowDuration;
            computeChecksumDuration      += other.computeChecksumDuration;
            compressWindowDuration       += other.compressWindowDuration;
            markerCount                  += other.markerCount;
            nonMarkerCount               += other.nonMarkerCount;
            realMarkerCount              += other.realMarkerCount;
        }

    public:
        size_t falsePositiveCount{ 0 };
        double blockFinderDuration{ 0 };
        double decodeDuration{ 0 };
        double decodeDurationInflateWrapper{ 0 };
        double decodeDurationIsal{ 0 };
        double appendDuration{ 0 };
        double applyWindowDuration{ 0 };
        double computeChecksumDuration{ 0 };
        double compressWindowDuration{ 0 };
        uint64_t markerCount{ 0 };
        uint64_t nonMarkerCount{ 0 };
        uint64_t realMarkerCount{ 0 };
    };

public:
    explicit
    ChunkData( const Configuration& configurationToUse ) :
        configuration( configurationToUse )
    {
        setCRC32Enabled( configurationToUse.crc32Enabled );
    }

    ~ChunkData() = default;
    ChunkData() = default;
    ChunkData( ChunkData&& ) = default;
    ChunkData( const ChunkData& ) = delete;
    ChunkData& operator=( ChunkData&& ) = default;
    ChunkData& operator=( const ChunkData& ) = delete;

    [[nodiscard]] CompressionType
    windowCompressionType() const
    {
        if ( configuration.windowCompressionType ) {
            return *configuration.windowCompressionType;
        }
        /* Only bother with overhead-introducing compression for large chunk compression ratios. */
        return configuration.windowSparsity || ( decodedSizeInBytes * 8 > 2 * encodedSizeInBits )
               ? CompressionType::ZLIB
               : CompressionType::NONE;
    }

    void
    append( deflate::DecodedVector&& toAppend )
    {
        auto t0 = now();

        if ( crc32s.back().enabled() ) {
            crc32s.back().update( toAppend.data(), toAppend.size() );

            const auto t1 = now();
            statistics.computeChecksumDuration += duration( t0, t1 );
            t0 = t1;
        }

        BaseType::append( std::move( toAppend ) );
        statistics.appendDuration += duration( t0 );
    }

    void
    append( deflate::DecodedDataView const& toAppend )
    {
        auto t0 = now();

        /* Ignore data with markers. Those will be CRC32 computed inside @ref applyWindow. */
        if ( crc32s.back().enabled() ) {
            for ( const auto& buffer : toAppend.data ) {
                crc32s.back().update( buffer.data(), buffer.size() );
            }

            const auto t1 = now();
            statistics.computeChecksumDuration += duration( t0, t1 );
            t0 = t1;
        }

        BaseType::append( toAppend );
        statistics.appendDuration += duration( t0 );
    }

    void
    applyWindow( WindowView const&     window,
                 CompressionType const windowCompressionType )
    {
        const auto markerCount = dataWithMarkersSize();
        const auto tApplyStart = now();

        /* This is expensive! It adds 20-30% overhead for the FASTQ file! Therefore disable it.
         * The result for this statistics for:
         *     SRR22403185_2.fastq.gz
         *         Total decompressed bytes                 : 361'815'302
         *         Non-marker symbols                       :  16'270'407 (4.49688 %)
         *         Replaced marker symbol buffers           : 345'544'895 (95.5031 %)
         *         Actual marker symbol count in buffers    :  48'956'747 (14.168 %)
         *     silesia.tar.gz
         *         Total decompressed bytes                 : 211'957'760
         *         Non-marker symbols                       : 137'967'986 (65.0922 %)
         *         Replaced marker symbol buffers           :  73'989'774 (34.9078 %)
         *         Actual marker symbol count in buffers    :  22'555'742 (30.4849 %)
         *     4GiB-base64.gz
         *         Total decompressed bytes                 : 4'294'967'296
         *         Non-marker symbols                       : 4'272'350'980 (99.4734 %)
         *         Replaced marker symbol buffers           :    22'616'316 (0.526577 %)
         *         Actual marker symbol count in buffers    :       162'330 (0.717756 %)
         *     CTU-13-Dataset.tar.gz
         *         Total decompressed bytes                 : 79'747'543'040
         *         Non-marker symbols                       : 55'838'926'274 (70.0196 %)
         *         Replaced marker symbol buffers           : 23'908'616'766 (29.9804 %)
         *         Actual marker symbol count in buffers    :  2'868'357'239 (11.9972 %)
         *     wikidata-20220103-all.json.gz
         *         Total decompressed bytes                 : 1'428'353'996'731
         *         Non-marker symbols                       :    23'033'941'599 (1.61262 %)
         *         Replaced marker symbol buffers           : 1'405'320'055'132 (98.3874 %)
         *         Actual marker symbol count in buffers    :   863'915'563'663 (61.4746 %)
         *
         * -> An alternative format that uses a mix of 8-bit and 16-bit and maybe a separate 1-bit buffer
         *    to store which byte is which, would reduce memory usage, and therefore also allocation
         *    overhead by 80%! Or maybe run-time encode it a la: <n 8-bit values> <8-bit value> ... <m 16-bit values>
         *    This would hopefully speed up window applying because hopefully long runs of 8-bit values could
         *    simply be memcopied and even runs of 16-bit values could be processed in a loop.
         *    This kind of compression would also add overhead though and it proabably would be too difficult
         *    to do inside deflate::Block, so it should probably be applied in post in
         *    ChunkData::append( DecodedDataViews ). This might be something that could be optimized with SIMD,
         *    the same applies to the equally necessary new ChunkData::applyWindow method.
         *    -> The count could be 7-bit so that the 8-th bit can be used to store the 8/16-bit value flag.
         *       In the worst case: interleaved 8-bit and 16-bit values, this would add an overhead of 25%:
         *       <n><8><n><16hi><16lo> <n><8>...
         *    Ideally a format that has no overhead even in the worst-case would be nice.
         *    This would be possible by using 4-bit values for <n> but then the maximum runlength would be 3-bit -> 7,
         *    which seems insufficient as it might lead to lots of slow execution branching in the applyWindow method.
         */
        static constexpr bool ENABLE_REAL_MARKER_COUNT = false;
        if constexpr ( ENABLE_REAL_MARKER_COUNT ) {
            statistics.realMarkerCount += countMarkerSymbols();
        }

        BaseType::applyWindow( window );

        const auto tApplyEnd = now();
        if ( markerCount > 0 ) {
            statistics.markerCount += markerCount;
            statistics.applyWindowDuration += duration( tApplyStart, tApplyEnd );
        }

        const auto alreadyProcessedSize = std::accumulate(
            crc32s.begin(), crc32s.end(), size_t( 0 ),
            [] ( const auto sum, const auto& crc32 ) { return sum + crc32.streamSize(); } );
        if ( crc32s.front().enabled() && ( alreadyProcessedSize < BaseType::dataSize() ) ) {
            /* Markers should only appear up to the first gzip footer because otherwise a new gzip stream
             * would have started. A new gzip stream must not contain markers because there are no unresolvable
             * back-references! Because of this, it is safe to only update the first CRC32.
             * Beware that we do not only have to compute the CRC32 of markers but also for data that has been
             * been converted from dataWithMarkers inside DecodedData::cleanUnmarkedData. */
            const auto toProcessSize = BaseType::dataSize() - alreadyProcessedSize;
            CRC32Calculator crc32;
            for ( auto it = DecodedData::Iterator( *this, 0, toProcessSize ); static_cast<bool>( it ); ++it ) {
                const auto [buffer, size] = *it;
                crc32.update( buffer, size );
            }
            crc32s.front().prepend( crc32 );

            statistics.computeChecksumDuration += duration( tApplyEnd );
        }

        /* Replace markers in and compress the resulting fully-resolved window provided by each subchunk,
         * i.e., at the end of each subchunk. In benchmarks with random base64 data and ISA-L, this takes
         * roughly 0.5 ms per 32 KiB window (0.048s for 97 compressed windows). */
        const auto tWindowCompressionStart = now();
        size_t decodedOffsetInBlock{ 0 };
        for ( auto& subchunk : m_subchunks ) {
            decodedOffsetInBlock += subchunk.decodedSize;

            if ( !subchunk.window ) {
                auto subchunkWindow = getWindowAt( window, decodedOffsetInBlock );
                /* Set unused symbols to 0 to increase compressibility. */
                if ( subchunkWindow.size() == subchunk.usedWindowSymbols.size() ) {
                    for ( size_t i = 0; i < subchunkWindow.size(); ++i ) {
                        if ( !subchunk.usedWindowSymbols[i] ) {
                            subchunkWindow[i] = 0;
                        }
                    }
                }
                subchunk.window = std::make_shared<Window>( std::move( subchunkWindow ), windowCompressionType );
            }
            subchunk.usedWindowSymbols = std::vector<bool>();  // Free memory!

            /* Count lines if requested. */
            if ( configuration.newlineCharacter && !subchunk.newlineCount ) {
                size_t newlineCount = 0;
                using rapidgzip::deflate::DecodedData;
                for ( auto it = DecodedData::Iterator( *this, subchunk.decodedOffset, subchunk.decodedSize );
                      static_cast<bool>( it ); ++it )
                {
                    const auto& [buffer, size] = *it;
                    newlineCount += std::count( reinterpret_cast<const char*>( buffer ),
                                                reinterpret_cast<const char*>( buffer ) + size,
                                                configuration.newlineCharacter.value() );
                }
                subchunk.newlineCount = newlineCount;
            }
        }
        statistics.compressWindowDuration += duration( tWindowCompressionStart );

        /* Check that it counts as fully post-processed from here on. */
        if ( !hasBeenPostProcessed() ) {
            std::stringstream message;
            message << "[Info] Chunk is not recognized as post-processed even though it has been!\n"
                    << "[Info]    Subchunks : " << m_subchunks.size() << "\n"
                    << "[Info]    Contains markers : " << containsMarkers() << "\n";
            for ( auto& subchunk : m_subchunks ) {
                if ( subchunk.hasBeenPostProcessed( configuration.newlineCharacter.has_value() ) ) {
                    continue;
                }
                message << "[Info] Subchunk is not recognized as post-processed even though it has been!\n"
                        << "[Info]    Has window : " << static_cast<bool>( subchunk.window ) << "\n"
                        << "[Info]    Used window symbols empty : " << subchunk.usedWindowSymbols.empty() << "\n"
                        << "[Info]    Has newline count : " << subchunk.newlineCount.has_value() << "\n";
                if ( configuration.newlineCharacter.has_value() ) {
                    message << "[Info]    Newline character : " << static_cast<int>( subchunk.newlineCount.value() )
                            << "\n";
                }
            }
        #ifdef RAPIDGZIP_FATAL_PERFORMANCE_WARNINGS
            throw std::logic_error( std::move( message ).str() );
        #else
            std::cerr << std::move( message ).str();
        #endif
        }
    }

    [[nodiscard]] bool
    matchesEncodedOffset( size_t offset ) const noexcept
    {
        if ( maxEncodedOffsetInBits == std::numeric_limits<size_t>::max() ) {
            return offset == encodedOffsetInBits;
        }
        return ( encodedOffsetInBits <= offset ) && ( offset <= maxEncodedOffsetInBits );
    }

    void
    setEncodedOffset( size_t offset );

    void
    setSubchunks( std::vector<Subchunk>&& newSubchunks )
    {
        m_subchunks = std::move( newSubchunks );
    }

    /**
     * @note Probably should not be called internally because it is allowed to be shadowed by a child class method.
     */
    void
    finalize( size_t newEncodedEndOffsetInBits )
    {
        const auto oldMarkerSize = BaseType::dataWithMarkersSize();
        cleanUnmarkedData();
        const auto toProcessSize = oldMarkerSize - BaseType::dataWithMarkersSize();
        if ( toProcessSize > 0 ) {
            const auto tComputeHashStart = now();

            CRC32Calculator crc32;
            /* Iterate over contiguous chunks of memory. */
            for ( auto it = DecodedData::Iterator( *this, 0, toProcessSize ); static_cast<bool>( it ); ++it ) {
                const auto [buffer, size] = *it;
                crc32.update( buffer, size );
            }
            /* Note that the data with markers ought not cross footer boundaries because after a footer,
             * a new gzip stream begins, which should be known to not contain any unresolvable backreferences.
             * That's why we can simply merge the CRC32 for the cleaned data with the first CRC32. */
            crc32s.front().prepend( crc32 );

            statistics.computeChecksumDuration += duration( tComputeHashStart );
        }

        statistics.nonMarkerCount += dataSize();

        encodedEndOffsetInBits = newEncodedEndOffsetInBits;
        encodedSizeInBits = newEncodedEndOffsetInBits - encodedOffsetInBits;
        decodedSizeInBytes = BaseType::size();

        if ( m_subchunks.empty() ) {
            m_subchunks = split( configuration.splitChunkSize );
        }
    }

    /**
     * Appends a deflate block boundary.
     * @return true if it was appended, false if the last boundary is identical to the given one.
     */
    bool
    appendDeflateBlockBoundary( const size_t encodedOffset,
                                const size_t decodedOffset )
    {
        if ( blockBoundaries.empty()
             || ( blockBoundaries.back().encodedOffset != encodedOffset )
             || ( blockBoundaries.back().decodedOffset != decodedOffset ) )
        {
            blockBoundaries.emplace_back( BlockBoundary{ encodedOffset, decodedOffset } );
            return true;
        }
        return false;
    }

    /**
     * Appends generic footer information at the given offset.
     */
    void
    appendFooter( const Footer& footer )
    {
        footers.emplace_back( footer );

        const auto wasEnabled = crc32s.back().enabled();
        crc32s.emplace_back();
        crc32s.back().setEnabled( wasEnabled );
    }

    void
    setCRC32Enabled( bool enabled )
    {
        configuration.crc32Enabled = enabled;
        for ( auto& calculator : crc32s ) {
            calculator.setEnabled( enabled );
        }
    }

    /**
     * @return When true is returned, @ref GzipChunkFetcher will queue the call to @ref applyWindow in the thread pool.
     *         After the call to @ref applyWindow, this function must return true!
     */
    [[nodiscard]] bool
    hasBeenPostProcessed() const
    {
        return !m_subchunks.empty() && !containsMarkers()
               && std::all_of( m_subchunks.begin(), m_subchunks.end(), [this] ( const auto& subchunk ) {
                      return subchunk.hasBeenPostProcessed( configuration.newlineCharacter.has_value() );
                  } );
    }

    [[nodiscard]] const std::vector<Subchunk>&
    subchunks() const noexcept
    {
        return m_subchunks;
    }

    /**
     * Chunks smaller than the returned value should not be created. In practice, this currently means that
     * such small chunks are appended to the previous one. This means however that some chunks can grow
     * larger than configuration.splitChunkSize.
     */
    [[nodiscard]] constexpr size_t
    minimumSplitChunkSize() const
    {
        return configuration.splitChunkSize / 4U;
    }

    /**
     * Implement a kind of virtual method by using a std::function member because making ChunkData polymorphic
     * had catastrophic impact on the performance for unknown reason.
     */
    [[nodiscard]] deflate::DecodedVector
    getWindowAt( WindowView const& previousWindow,
                 size_t            skipBytes ) const
    {
        return m_getWindowAt
                ? m_getWindowAt( *this, previousWindow, skipBytes )
                : DecodedData::getWindowAt( previousWindow, skipBytes );
    }

protected:
    [[nodiscard]] std::vector<Subchunk>
    split( size_t spacing ) const;

public:
    size_t encodedOffsetInBits{ std::numeric_limits<size_t>::max() };
    size_t encodedSizeInBits{ 0 };

    /* This should only be evaluated when it is unequal std::numeric_limits<size_t>::max() and unequal
     * Base::encodedOffsetInBits. Then, [Base::encodedOffsetInBits, maxEncodedOffsetInBits] specifies a valid range
     * for the block offset. Such a range might happen for finding uncompressed deflate blocks because of the
     * byte-padding. */
    size_t maxEncodedOffsetInBits{ std::numeric_limits<size_t>::max() };
    /* Initialized with size() after thread has finished writing into ChunkData. Redundant but avoids a lock
     * because the marker replacement will momentarily lead to different results returned by size! */
    size_t decodedSizeInBytes{ 0 };
    /* This is currently only set in @ref finalize and used in @ref setEncodedOffset to initialize
     * @ref encodedSizeInBits. */
    size_t encodedEndOffsetInBits{ std::numeric_limits<size_t>::max() };

    Configuration configuration;

    /* Decoded offsets are relative to the decoded offset of this ChunkData because that might not be known
     * during first-pass decompression. */
    std::vector<BlockBoundary> blockBoundaries;
    std::vector<Footer> footers;
    /* There will be ( footers.size() + 1 ) CRC32 calculators. */
    std::vector<CRC32Calculator> crc32s{ std::vector<CRC32Calculator>( 1 ) };

    Statistics statistics{};

    bool stoppedPreemptively{ false };

protected:
    /**
     * Takes ChunkData& as first argument instead of capturing this in order to avoid having to implement
     * custom move and copy constructors.
     */
    std::function<deflate::DecodedVector( const ChunkData&, WindowView const&, size_t )> m_getWindowAt;
    std::vector<Subchunk> m_subchunks;
};


inline std::ostream&
operator<<( std::ostream&    out,
            const ChunkData& chunk )
{
    out << "ChunkData{\n";
    out << "  encodedOffsetInBits: " << chunk.encodedOffsetInBits << "\n";
    out << "  encodedSizeInBits: " << chunk.encodedSizeInBits << "\n";
    out << "  maxEncodedOffsetInBits: " << chunk.maxEncodedOffsetInBits << "\n";
    out << "  decodedSizeInBytes: " << chunk.decodedSizeInBytes << "\n";
    out << "  blockBoundaries: { ";
    for ( const auto& boundary : chunk.blockBoundaries ) {
        out << boundary.encodedOffset << ":" << boundary.decodedOffset << ", ";
    }
    out << "}\n";
    out << "  footers: { ";
    for ( const auto& footer : chunk.footers ) {
        out << footer.blockBoundary.encodedOffset << ":" << footer.blockBoundary.decodedOffset << ", ";
    }
    out << "}\n";
    out << "}\n";
    return out;
}


inline void
ChunkData::setEncodedOffset( size_t offset )
{
    if ( !matchesEncodedOffset( offset ) ) {
        throw std::invalid_argument( "The real offset to correct to should lie inside the offset range!" );
    }

    if ( encodedEndOffsetInBits == std::numeric_limits<size_t>::max() ) {
        throw std::invalid_argument( "Finalize must be called before setEncodedOffset!" );
    }

    if ( encodedEndOffsetInBits < offset ) {
        std::stringstream message;
        message << "The chunk start " << offset << " must not be after the chunk end " << encodedEndOffsetInBits << "!";
        throw std::invalid_argument( std::move( message ).str() );
    }

    encodedSizeInBits = encodedEndOffsetInBits - offset;
    encodedOffsetInBits = offset;
    maxEncodedOffsetInBits = offset;

    /* Adjust the encoded offset of the first subchunk because it may have been a range at the time of splitting. */
    if ( !m_subchunks.empty() ) {
        const auto nextSubchunk = std::next( m_subchunks.begin() );
        const auto nextOffset = nextSubchunk == m_subchunks.end() ? encodedEndOffsetInBits : nextSubchunk->encodedOffset;
        m_subchunks.front().encodedOffset = offset;
        m_subchunks.front().encodedSize = nextOffset - offset;
    }
}


[[nodiscard]] inline std::vector<ChunkData::Subchunk>
ChunkData::split( [[maybe_unused]] const size_t spacing ) const
{
    if ( encodedEndOffsetInBits == std::numeric_limits<size_t>::max() ) {
        throw std::invalid_argument( "Finalize must be called before splitting the chunk!" );
    }

    if ( spacing == 0 ) {
        throw std::invalid_argument( "Spacing must be a positive number of bytes." );
    }

    if ( ( encodedSizeInBits == 0 ) && ( decodedSizeInBytes == 0 ) ) {
        return {};
    }

    const auto nBlocks = static_cast<size_t>( std::round( static_cast<double>( decodedSizeInBytes )
                                                          / static_cast<double>( spacing ) ) );
    Subchunk wholeChunkAsSubchunk;
    wholeChunkAsSubchunk.encodedOffset = encodedOffsetInBits;
    wholeChunkAsSubchunk.decodedOffset = 0;
    wholeChunkAsSubchunk.encodedSize = encodedSizeInBits;
    wholeChunkAsSubchunk.decodedSize = decodedSizeInBytes;
    /* blockBoundaries does not contain the first block begin but all thereafter including the boundary after
     * the last block, i.e., the begin of the next deflate block not belonging to this ChunkData. */
    if ( ( nBlocks <= 1 ) || blockBoundaries.empty() ) {
        return { wholeChunkAsSubchunk };
    }

    /* The idea for partitioning is: Divide the size evenly and into subchunks and then choose the block boundary
     * that is closest to that value. */
    const auto perfectSpacing = static_cast<double>( decodedSizeInBytes ) / static_cast<double>( nBlocks );

    std::vector<Subchunk> result;
    result.reserve( nBlocks + 1 );

    BlockBoundary lastBoundary{ encodedOffsetInBits, 0 };
    /* The first and last boundaries are static, so we only need to find nBlocks - 1 further boundaries. */
    for ( size_t iSubchunk = 1; iSubchunk < nBlocks; ++iSubchunk ) {
        const auto perfectDecompressedOffset = static_cast<size_t>( static_cast<double>( iSubchunk ) * perfectSpacing );

        const auto isCloser =
            [perfectDecompressedOffset] ( const auto& b1, const auto& b2 )
            {
                return absDiff( b1.decodedOffset, perfectDecompressedOffset )
                       < absDiff( b2.decodedOffset, perfectDecompressedOffset );
            };
        auto closest = std::min_element( blockBoundaries.begin(), blockBoundaries.end(), isCloser );

        /* If there are duplicate decodedOffset, min_element returns the first, so skip over empty blocks (pigz).
         * Using the last block with the same decodedOffset makes handling the last block after this for-loop easier. */
        while ( ( std::next( closest ) != blockBoundaries.end() )
                && ( closest->decodedOffset == std::next( closest )->decodedOffset ) )
        {
            ++closest;
        }

        /* For very small spacings, the same boundary might be found twice. Avoid empty subchunks because of that. */
        if ( closest->decodedOffset <= lastBoundary.decodedOffset ) {
            continue;
        }

        if ( closest->encodedOffset <= lastBoundary.encodedOffset ) {
            throw std::logic_error( "If the decoded offset is strictly larger than so must be the encoded one!" );
        }

        Subchunk subchunk;
        subchunk.encodedOffset = lastBoundary.encodedOffset;
        subchunk.decodedOffset = result.empty() ? 0 : result.back().decodedOffset + result.back().decodedSize;
        subchunk.encodedSize = closest->encodedOffset - lastBoundary.encodedOffset;
        subchunk.decodedSize = closest->decodedOffset - lastBoundary.decodedOffset;
        result.emplace_back( subchunk );
        lastBoundary = *closest;
    }

    if ( lastBoundary.decodedOffset > decodedSizeInBytes ) {
        throw std::logic_error( "There should be no boundary outside of the chunk range!" );
    }
    if ( ( lastBoundary.decodedOffset < decodedSizeInBytes ) || result.empty() ) {
        /* Create the last subchunk from lastBoundary and the chunk end. */
        Subchunk subchunk;
        subchunk.encodedOffset = lastBoundary.encodedOffset,
        subchunk.decodedOffset = result.empty() ? 0 : result.back().decodedOffset + result.back().decodedSize;
        subchunk.encodedSize = encodedEndOffsetInBits - lastBoundary.encodedOffset,
        subchunk.decodedSize = decodedSizeInBytes - lastBoundary.decodedOffset,
        result.emplace_back( subchunk );
    } else if ( lastBoundary.decodedOffset == decodedSizeInBytes ) {
        /* Enlarge the last subchunk encoded size to also encompass the empty blocks before the chunk end.
         * Assuming that blockBoundaries contain the boundary at the chunk end and knowing that the loop
         * above always searches for the last boundary with the same decodedOffset, this branch shouldn't happen. */
        result.back().encodedSize = encodedEndOffsetInBits - result.back().encodedOffset;
    }

    if ( encodedEndOffsetInBits - encodedOffsetInBits != encodedSizeInBits ) {
        std::stringstream message;
        message << "The offset: " << encodedOffsetInBits << ", size: " << encodedSizeInBits << ", and end offset: "
                << encodedEndOffsetInBits << " are inconsistent!";
        throw std::logic_error( std::move( message ).str() );
    }

    const auto subchunkEncodedSizeSum =
        std::accumulate( result.begin(), result.end(), size_t( 0 ),
                         [] ( size_t sum, const auto& block ) { return sum + block.encodedSize; } );
    const auto subchunkDecodedSizeSum =
        std::accumulate( result.begin(), result.end(), size_t( 0 ),
                         [] ( size_t sum, const auto& block ) { return sum + block.decodedSize; } );
    if ( ( subchunkEncodedSizeSum != encodedSizeInBits ) || ( subchunkDecodedSizeSum != decodedSizeInBytes ) ) {
        std::stringstream message;
        message << "[Warning] Block splitting was unsuccessful. This might result in higher memory usage but is "
                << "otherwise harmless. Please report this performance bug with a reproducing example.\n"
                << "  subchunkEncodedSizeSum: " << subchunkEncodedSizeSum << "\n"
                << "  encodedSizeInBits     : " << encodedSizeInBits      << "\n"
                << "  subchunkDecodedSizeSum: " << subchunkDecodedSizeSum << "\n"
                << "  decodedSizeInBytes    : " << decodedSizeInBytes     << "\n";
    #ifdef RAPIDGZIP_FATAL_PERFORMANCE_WARNINGS
        throw std::logic_error( std::move( message ).str() );
    #else
        std::cerr << std::move( message ).str();
    #endif
        return { wholeChunkAsSubchunk };  // fallback without any splitting done at all
    }

    return result;
}


/**
 * m rapidgzip && src/tools/rapidgzip -v -d -c -P 0 4GiB-base64.gz | wc -c
 * Non-polymorphic: Decompressed in total 4294967296 B in 1.49444 s -> 2873.96 MB/s
 * With virtual ~DecodedData() = default: Decompressed in total 4294967296 B in 3.58325 s -> 1198.62 MB/s
 * I don't know why it happens. Maybe it affects inline of function calls or moves of instances.
 */
static_assert( !std::is_polymorphic_v<ChunkData>, "Simply making it polymorphic halves performance!" );


#if defined( HAVE_VMSPLICE )
/**
 * Tries to use writeAllSpliceUnsafe and, if successful, also extends lifetime by adding the block data
 * shared_ptr into a list.
 *
 * @note Limitations:
 *  - To avoid querying the pipe buffer size, it is only done once. This might introduce subtle errors when it is
 *    dynamically changed after this point.
 *  - The lifetime can only be extended on block granularity even though chunks would be more suited.
 *    This results in larger peak memory than strictly necessary.
 *  - In the worst case we would read only 1B out of each block, which would extend the lifetime
 *    of thousands of large blocks resulting in an out of memory issue.
 *    - This would only be triggerable by using the API. The current CLI and not even the Python
 *      interface would trigger this because either they don't splice to a pipe or only read
 *      sequentially.
 * @note It *does* account for pages to be spliced into yet another pipe buffer by waiting for buffer size
 *       amount of data being written before freeing, and likely reusing, the memory.
 */
[[nodiscard]] inline int
writeAllSplice( [[maybe_unused]] const int                         outputFileDescriptor,
                [[maybe_unused]] const std::shared_ptr<ChunkData>& chunkData,
                [[maybe_unused]] const std::vector<::iovec>&       buffersToWrite )
{
    return SpliceVault::getInstance( outputFileDescriptor ).first->splice( buffersToWrite, chunkData );
}
#endif  // HAVE_VMSPLICE


[[nodiscard]] inline int
writeAll( const std::shared_ptr<ChunkData>& chunkData,
          const int                         outputFileDescriptor,
          const size_t                      offsetInBlock,
          const size_t                      dataToWriteSize )
{
    if ( ( outputFileDescriptor < 0 ) || ( dataToWriteSize == 0 ) ) {
        return 0;
    }

#ifdef HAVE_VMSPLICE
    const auto buffersToWrite = toIoVec( *chunkData, offsetInBlock, dataToWriteSize );
    const auto errorCode = writeAllSplice( outputFileDescriptor, chunkData, buffersToWrite );
    if ( errorCode != 0 ) {
        return writeAllToFdVector( outputFileDescriptor, buffersToWrite );
    }
#else
    using rapidgzip::deflate::DecodedData;

    for ( auto it = DecodedData::Iterator( *chunkData, offsetInBlock, dataToWriteSize );
          static_cast<bool>( it ); ++it )
    {
        const auto& [buffer, size] = *it;
        const auto errorCode = writeAllToFd( outputFileDescriptor, buffer, size );
        if ( errorCode != 0 ) {
            return errorCode;
        }
    }
#endif

    return 0;
}


/**
 * This subclass of @ref ChunkData only counts the decompressed bytes and does not store them.
 */
struct ChunkDataCounter final :
    public ChunkData
{
    explicit
    ChunkDataCounter( const Configuration& configurationToUse ) :
        ChunkData( configurationToUse )
    {
        /**
         * The internal index will only contain the offsets and empty windows.
         * The index should not be exported when this is used.
         * Return a dummy window so that decoding can be resumed after stopping.
         */
        m_getWindowAt = [] ( const ChunkData&, WindowView const&, size_t ) {
            return deflate::DecodedVector( deflate::MAX_WINDOW_SIZE, 0 );
        };
    }

    ChunkDataCounter()
    {
        m_getWindowAt = [] ( const ChunkData&, WindowView const&, size_t ) {
            return deflate::DecodedVector( deflate::MAX_WINDOW_SIZE, 0 );
        };
    }

    ~ChunkDataCounter() = default;
    ChunkDataCounter( ChunkDataCounter&& ) = default;
    ChunkDataCounter( const ChunkDataCounter& ) = delete;
    ChunkDataCounter& operator=( ChunkDataCounter&& ) = default;
    ChunkDataCounter& operator=( const ChunkDataCounter& ) = delete;

    void
    append( deflate::DecodedVector&& toAppend )
    {
        decodedSizeInBytes += toAppend.size();
    }

    void
    append( deflate::DecodedDataView const& toAppend )
    {
        decodedSizeInBytes += toAppend.size();
    }

    void
    finalize( size_t newEncodedEndOffsetInBits )
    {
        encodedEndOffsetInBits = newEncodedEndOffsetInBits;
        encodedSizeInBits = encodedEndOffsetInBits - encodedOffsetInBits;
        /* Do not overwrite decodedSizeInBytes like is done in the base class
         * because DecodedData::size() would return 0! Instead, it is updated inside append. */

        m_subchunks = split( configuration.splitChunkSize );
    }

    /**
     * No splitting necessary for memory reduction because we don't hold the results anyway.
     */
    [[nodiscard]] std::vector<Subchunk>
    split( [[maybe_unused]] const size_t spacing ) const
    {
        return ChunkData::split( std::numeric_limits<size_t>::max() );
    }
};
}  // namespace rapidgzip
