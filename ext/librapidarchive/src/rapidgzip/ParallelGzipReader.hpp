#pragma once

#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <core/AffinityHelpers.hpp>
#include <core/BlockMap.hpp>
#include <core/common.hpp>
#include <filereader/FileReader.hpp>
#include <filereader/Shared.hpp>

#ifdef WITH_PYTHON_SUPPORT
    #include <filereader/Python.hpp>
    #include <filereader/Standard.hpp>
#endif

#include "gzip/crc32.hpp"
#include "gzip/gzip.hpp"
#include "GzipChunkFetcher.hpp"
#include "GzipBlockFinder.hpp"
#include "IndexFileFormat.hpp"


namespace rapidgzip
{
#ifdef WITH_PYTHON_SUPPORT
enum class IOReadMethod : uint8_t
{
    SEQUENTIAL = 0,
    PREAD = 1,
    LOCKED_READ_AND_SEEK = 2,
};

[[nodiscard]] static UniqueFileReader
wrapFileReader( UniqueFileReader&& fileReader,
                IOReadMethod       ioReadMethod )
{
    switch ( ioReadMethod )
    {
    case IOReadMethod::SEQUENTIAL:
        return std::make_unique<SinglePassFileReader>( std::move( fileReader ) );
    case IOReadMethod::PREAD:
    case IOReadMethod::LOCKED_READ_AND_SEEK:
    {
        auto sharedFile = ensureSharedFileReader( std::move( fileReader ) );
        sharedFile->setUsePread( ioReadMethod == IOReadMethod::PREAD );
        return sharedFile;
    }
    }
    return std::move( fileReader );
}
#endif  // WITH_PYTHON_SUPPORT


/**
 * @note Calls to this class are not thread-safe! Even though they use threads to evaluate them in parallel.
 */
template<typename T_ChunkData = ChunkData>
class ParallelGzipReader final :
    public FileReader
{
public:
    using ChunkData = T_ChunkData;
    using ChunkConfiguration = typename ChunkData::Configuration;
    /**
     * The fetching strategy should support parallelization via prefetching for sequential accesses while
     * avoiding a lot of useless prefetches for random or multi-stream sequential accesses like those occurring
     * via ratarmount.
     * The fetching strategy does not have to and also should not account for backward and strided accesses
     * because the prefetch and cache units are very large and striding or backward accessing over multiple
     * megabytes should be extremely rare.
     */
    using ChunkFetcher = rapidgzip::GzipChunkFetcher<FetchingStrategy::FetchMultiStream, ChunkData>;
    using BlockFinder = typename ChunkFetcher::BlockFinder;
    using BitReader = gzip::BitReader;
    using WriteFunctor = std::function<void ( const std::shared_ptr<ChunkData>&, size_t, size_t )>;
    using Window = WindowMap::Window;

    struct NewlineOffset
    {
        uint64_t lineOffset{ 0 };
        uint64_t uncompressedOffsetInBytes{ 0 };
    };

public:
    /**
     * Quick benchmarks for spacing on AMD Ryzen 3900X 12-core.
     *
     * @verbatim
     * base64 /dev/urandom | head -c $(( 4 * 1024 * 1024 * 1024 )) > 4GiB-base64
     * gzip 4GiB-base64
     *
     * function benchmarkWc()
     * {
     *     for chunkSize in 128 $(( 1*1024 )) $(( 2*1024 )) $(( 4*1024 )) $(( 8*1024 )) $(( 16*1024 )) $(( 32*1024 )); do
     *         echo "Chunk Size: $chunkSize KiB"
     *         for i in $( seq 5 ); do
     *             src/tools/rapidgzip --chunk-size $chunkSize -v -P 0 -d -c "$1" 2>rapidgzip.log | wc -c
     *             grep "Decompressed in total" rapidgzip.log
     *         done
     *     done
     * }
     *
     * m rapidgzip
     * benchmarkWc 4GiB-base64.gz
     *
     *
     * spacing | bandwidth / (MB/s) | file read multiplier
     * --------+--------------------+----------------------
     * 128 KiB | ~1250              | 2.08337
     *   1 MiB | ~3500              | 1.13272
     *   2 MiB | ~3800              | 1.06601
     *   4 MiB | ~4000              | 1.03457
     *   8 MiB | ~4200              | 1.0169
     *  16 MiB | ~4400              | 1.00799
     *  32 MiB | ~4100              | 1.00429
     * @endverbatim
     *
     * For higher chunk sizes, the bandwidths become very unstable,
     * probably because even work division becomes a problem relative to the file size.
     * Furthermore, caching behavior might worsen for larger chunk sizes.
     *
     * @verbatim
     * wget https://sun.aei.polsl.pl/~sdeor/corpus/silesia.zip
     * mkdir -p silesia && ( cd silesia && unzip ../silesia.zip )
     * tar -cf silesia.tar silesia/  # 211957760 B -> 212 MB, 203 MiB, gzip 66 MiB -> compression factor: 3.08
     * for (( i=0; i<40; ++i )); do cat 'silesia.tar'; done | pigz > 40xsilesia.tar.gz
     * m rapidgzip
     * benchmarkWc 40xsilesia.tar.gz
     *
     * spacing | bandwidth / (MB/s)
     * --------+--------------------
     * 128 KiB | ~1450
     *   1 MiB | ~2500
     *   2 MiB | ~2800
     *   4 MiB | ~3400
     *   8 MiB | ~3800
     *  16 MiB | ~4100
     *  32 MiB | ~4100
     * @endverbatim
     *
     * Beware, on 2xAMD EPYC CPU 7702, when decoding with more than 64 cores, the optimum is
     * at 2 MiB instead of 4-8 MiB! Maybe these are NUMA domain + caching issues combined?
     *
     * AMD Ryzen 3900X Caches:
     *  - L1: 64 kiB (50:50 instruction:cache) per core -> 768 kiB
     *  - L2: 512 kiB per core -> 6 MiB
     *  - L3: 64 MiB shared (~5.3 MiB per core)
     *  - RAM: 2x16GiB DIMM DDR4 3600 MHz (0.3 ns), 2x32GiB DIMM DDR4 3600 MHz (0.3 ns)
     *
     * AMD EPYC CPU 7702:
     *  - L1: 64 kiB (50:50 instruction:cache) per core -> 4 MiB
     *  - L2: 512 kiB per core -> 32 MiB
     *  - L3: 256 MiB shared (4 MiB per core)
     *
     * -> That EPYC processor is the same generation Zen 2 and therefore has identical L1 and L2 caches
     *    and the L3 cache size is even higher, so it must be a NUMA issue.
     *
     * Non-compressible data is a special case because it only needs to do a memcpy.
     *
     * @verbatim
     * head -c $(( 4 * 1024 * 1024 * 1024 )) /dev/urandom | gzip > 4GiB-random.gz
     * m rapidgzip
     * benchmarkWc 4GiB-random.gz
     *
     * spacing | bandwidth / (MB/s) | file read multiplier
     * --------+--------------------+----------------------
     * 128 KiB | ~1300              | 2.00049
     *   1 MiB | ~3400              | 1.12502
     *   2 MiB | ~3900              | 1.06253
     *   4 MiB | ~4000              | 1.03129
     *   8 MiB | ~4100              | 1.01567
     *  16 MiB | ~4200              | 1.00786
     *  32 MiB | ~4200              | 1.00396
     * @endverbatim
     *
     * Another set of benchmarks that exclude the bottleneck for writing the results to a pipe by
     * using the rapidgzip option --count-lines. Note that in contrast to pugz, the decompressed
     * blocks are still processed in sequential order. Processing them out of order by providing
     * a map-reduce like interface might accomplish even more speedups.
     *
     * @verbatim
     * m rapidgzip
     * for chunkSize in 128 $(( 1*1024 )) $(( 2*1024 )) $(( 4*1024 )) $(( 8*1024 )) $(( 16*1024 )) $(( 32*1024 )); do
     *     echo "Chunk Size: $chunkSize KiB"
     *     for i in $( seq 5 ); do
     *         src/tools/rapidgzip -v --chunk-size $chunkSize -P 0 --count-lines 4GiB-base64.gz 2>rapidgzip.log
     *         grep "Decompressed in total" rapidgzip.log
     *     done
     * done
     *
     * spacing | bandwidth / (MB/s)
     * --------+--------------------
     * 128 KiB | ~1500
     *   1 MiB | ~4600
     *   2 MiB | ~5000
     *   4 MiB | ~5400
     *   8 MiB | ~5400
     *  16 MiB | ~5100
     *  32 MiB | ~4900
     * @endverbatim
     *
     * The factor 2 amount of read data can be explained with the BitReader always buffering 128 KiB!
     * Therefore if the work chunk is too small, it leads to this problem.
     * @note Beware the actual result of wc -l! With the wrong vmsplice usage, it returned random results
     *       for chunk sizes smaller than 4 MiB or even for higher chunk sizes with alternative malloc
     *       implementations like mimalloc.
     * @note The optimum at ~8 MiB for incompressible data vs ~4 MiB for base64 data with a compression
     *       ratio ~1.3 might be explainable with a roughly equal decompressed block size. In general,
     *       we would like the chunk size to be measured in decompressed data because the decompressed
     *       bandwidth is much more stable than the compressed bandwidth over a variety of data.
     * @todo We might be able to reduce this overhead by buffering up to untilOffset and then
     *       only increase the buffer in much smaller steps, e.g., 8 KiB.
     *       This might actually be easy to implement by making the BitReader chunk size adjustable.
     * @todo Possibly increase the chunk size to 4 or 8 MiB again after implementing an out of memory
     *       guard for high compression ratios so that CTU-13-Dataset.tar.gz can be decompressed with
     *       less than 30 GB of RAM!
     *       Rebenchmark of course whether it makes sense or not anymore. E.g., speeding up the block
     *       finder might enable smaller chunk sizes.
     *
     * @verbatim
     * for (( i=0; i<10; ++i )); do cat 'silesia.tar'; done | lbzip2 > 10xsilesia.tar.bz2
     * stat --format=%s -L 10xsilesia.tar.bz2
     *     546 315 457
     * benchmarkWc 10xsilesia.tar.bz2
     *
     * spacing | bandwidth / (MB/s)
     * --------+--------------------
     * 128 KiB | ~370
     *   1 MiB | ~410
     *   2 MiB | ~510
     *   4 MiB | ~600 <-
     *   8 MiB | ~560
     *  16 MiB | ~540
     *  32 MiB | ~550
     * @endverbatim
     *
     * @verbatim
     * benchmarkWc silesia.tar.bz2
     * stat --format=%s -L silesia.tar.bz2
     *     54 591 465 = 52.06 MiB
     *
     * spacing | bandwidth / (MB/s)
     * --------+--------------------
     * 128 KiB | ~340
     *   1 MiB | ~400 <-
     *   2 MiB | ~400
     *   4 MiB | ~400
     *   8 MiB | ~400
     *  16 MiB | ~400
     *  32 MiB | ~400
     * @endverbatim
     *
     * There simply is not enough work to distribute. That's why it is slow for larger chunk sizes.
     * For smaller chunk sizes it becomes slow because some chunks won't find anything to decode but they
     * still count towards the maximum cached chunk size.
     * @todo They shouldn't count towards that limit because they don't consume much memory anyway.
     *       Maybe test those somehow and move them into a different lookup cache, or simply don't count them.
     *       The latter might be expensive if they become too many and if it isn't a simply bool check.
     *       Unfortunately, it isn't even easily possible to check for exception. We would have to call future::get()
     *       in a try-catch-block and repackage the result thereafter or change the ChunkFetcher interface to not
     *       return blocks. Lots of work. Or simply don't use chunk sizes smaller than 1 MiB because compressed bzip2
     *       should become much larger than 900kB.
     */
    explicit
    ParallelGzipReader( UniqueFileReader fileReader,
                        size_t           parallelization = 0,
                        uint64_t         chunkSizeInBytes = 4_Mi ) :
        m_chunkSizeInBytes( std::max( 8_Ki, chunkSizeInBytes ) ),
        m_sharedFileReader( ensureSharedFileReader( std::move( fileReader ) ) ),
        m_fetcherParallelization( parallelization == 0 ? availableCores() : parallelization ),
        m_startBlockFinder(
            [this] () {
                return std::make_unique<BlockFinder>(
                    UniqueFileReader( m_sharedFileReader->clone() ),
                    /* spacing in bytes */ m_chunkSizeInBytes );
            }
        )
    {
        setMaxDecompressedChunkSize( 20U * m_chunkSizeInBytes );

        const auto fileSize = m_sharedFileReader->size();
        if ( fileSize && ( m_chunkSizeInBytes * 2U * parallelization > *fileSize ) ) {
            /* Use roughly as many chunks as there is parallelization.
             * Multiply a factor of two, to give the thread pool more time to be filled out.
             * Bound the minimum chunk size because of the block finder overhead for gzip,
             * because <900kB chunks might not have any real work to do, and to avoid many threads
             * being started for very small files.
             * This formula is mostly optimized for silesia.tar.bz2.
             * Speed isn't that important for small gzip files because it decompresses many times faster.
             * In the first place, this implementation is intended towards very large files not small files. */
            m_chunkSizeInBytes =
                std::max( 512_Ki, ceilDiv( ceilDiv( *fileSize, 3U * parallelization ), 512_Ki ) * 512_Ki );
        }

        m_sharedFileReader->setStatisticsEnabled( m_statisticsEnabled );
        if ( !m_sharedFileReader->seekable() ) {
            /* The ensureSharedFileReader helper should wrap non-seekable file readers inside SinglePassFileReader. */
            throw std::logic_error( "BitReader should always be seekable even if the underlying file is not!" );
        }

        const auto& [lock, file] = m_sharedFileReader->underlyingFile();
        // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
        auto* const singlePassFileReader = dynamic_cast<SinglePassFileReader*>( file );
        if ( singlePassFileReader != nullptr ) {
            singlePassFileReader->setMaxReusableChunkCount(
                static_cast<size_t>(
                    std::ceil( static_cast<double>( parallelization ) * static_cast<double>( m_chunkSizeInBytes )
                               / static_cast<double>( SinglePassFileReader::CHUNK_SIZE ) ) ) );
            setKeepIndex( false );
        }
    }

#ifdef WITH_PYTHON_SUPPORT
    /* These constructor overloads are for easier construction in the Cython-interface.
     * For C++, the FileReader constructor would have been sufficient. */

    explicit
    ParallelGzipReader( int          fileDescriptor,
                        size_t       parallelization,
                        uint64_t     chunkSizeInBytes,
                        IOReadMethod ioReadMethod ) :
        ParallelGzipReader( wrapFileReader( std::make_unique<StandardFileReader>( fileDescriptor ), ioReadMethod ),
                            parallelization, chunkSizeInBytes )
    {}

    explicit
    ParallelGzipReader( const std::string& filePath,
                        size_t             parallelization,
                        uint64_t           chunkSizeInBytes,
                        IOReadMethod       ioReadMethod ) :
        ParallelGzipReader( wrapFileReader( std::make_unique<StandardFileReader>( filePath ), ioReadMethod ),
                            parallelization, chunkSizeInBytes )
    {}

    explicit
    ParallelGzipReader( PyObject*    pythonObject,
                        size_t       parallelization,
                        uint64_t     chunkSizeInBytes,
                        IOReadMethod ioReadMethod ) :
        ParallelGzipReader( wrapFileReader( std::make_unique<PythonFileReader>( pythonObject ), ioReadMethod ),
                            parallelization, chunkSizeInBytes )
    {}
#endif

    ~ParallelGzipReader() override
    {
        if ( m_showProfileOnDestruction && m_statisticsEnabled ) {
            std::stringstream out;
            out << std::boolalpha;
            out << "[ParallelGzipReader] Time spent:\n";
            out << "    Writing to output         : " << m_writeOutputTime << " s\n";
            out << "    Computing CRC32           : " << m_crc32Time << " s\n";
            out << "    Number of verified CRC32s : " << m_verifiedCRC32Count << "\n";
            out << "\nChunk Configuration:\n";
            out << "    CRC32 enabled      : " << m_crc32.enabled() << "\n";
            out << "    Window compression : " << ( m_windowCompressionType
                                                    ? toString( *m_windowCompressionType )
                                                    : std::string( "Default" ) ) << "\n";
            out << "    Window sparsity    : " << m_windowSparsity << "\n";
            std::cerr << std::endl;
        }
    }

    void
    setStatisticsEnabled( bool enabled )
    {
        m_statisticsEnabled = enabled;
        if ( m_chunkFetcher ) {
            m_chunkFetcher->setStatisticsEnabled( m_statisticsEnabled );
        }
        if ( m_sharedFileReader ) {
            m_sharedFileReader->setStatisticsEnabled( m_statisticsEnabled );
        }
    }

    /**
     * @note Only will work if m_statisticsEnabled is true.
     */
    void
    setShowProfileOnDestruction( bool showProfileOnDestruction )
    {
        m_showProfileOnDestruction = showProfileOnDestruction;
        if ( m_chunkFetcher ) {
            m_chunkFetcher->setShowProfileOnDestruction( m_showProfileOnDestruction );
        }
        if ( m_sharedFileReader ) {
            m_sharedFileReader->setShowProfileOnDestruction( m_showProfileOnDestruction );
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
        if ( !m_sharedFileReader || !m_sharedFileReader->seekable() ) {
            return false;
        }

        const auto& [lock, file] = m_sharedFileReader->underlyingFile();
        // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
        auto* const singlePassFileReader = dynamic_cast<SinglePassFileReader*>( file );
        return singlePassFileReader == nullptr;
    }

    void
    close() override
    {
        m_chunkFetcher.reset();
        m_blockFinder.reset();
        m_sharedFileReader.reset();
    }

    [[nodiscard]] bool
    closed() const override
    {
        return !m_sharedFileReader || m_sharedFileReader->closed();
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
        if ( !m_blockMap->finalized() ) {
            return std::nullopt;
        }
        return m_blockMap->back().second;
    }

    void
    clearerr() override
    {
        if ( m_sharedFileReader ) {
            m_sharedFileReader->clearerr();
        }
        m_atEndOfFile = false;
        throw std::invalid_argument( "Not fully tested!" );
    }

    /**
     * @param nBytesToRead This parameter can be performance-critical! Very small calls must be avoided because
     *        lots of checks such as for closed() can become expensive as they might require locking mutexes!
     *        Optimal inclusive ranges for number of bytes per call on Ryzen 3900X 12-core:
     *        - parallelization=1            : [8 KiB, 256 MiB]  no threading and therefore any value is alright!
     *        - parallelization=24           : [8 MiB, 256 MiB] smaller 128 B becomes unusably slow!
     *        - parallelization=24 with index: [32 KiB, 4 MiB]  smaller 128 B becomes unusably slow!
     *        Therefore, my recommendation would be to simply use 4 or 8 MiB,
     *        but [128 KiB, 256 MiB] is generally fine if you can live with ~20 % slowdowns.
     *        See also the comment for DECOMPRESSION_BUFFER_SIZE in benchmarkGzip.cpp!
     */
    [[nodiscard]] size_t
    read( char*  outputBuffer,
          size_t nBytesToRead ) override
    {
        return read( -1, outputBuffer, nBytesToRead );
    }

    /* Simpler file reader interface for Python-interfacing */

    // NOLINTBEGIN(misc-no-recursion)

    size_t
    read( const int    outputFileDescriptor = -1,
          char* const  outputBuffer         = nullptr,
          const size_t nBytesToRead         = std::numeric_limits<size_t>::max() )
    {
        const auto writeFunctor =
            [nBytesDecoded = uint64_t( 0 ), outputFileDescriptor, outputBuffer]
            ( const std::shared_ptr<ChunkData>& chunkData,
              size_t const                      offsetInBlock,
              size_t const                      dataToWriteSize ) mutable
            {
                if ( dataToWriteSize == 0 ) {
                    return;
                }

                const auto errorCode = writeAll( chunkData, outputFileDescriptor, offsetInBlock, dataToWriteSize );
                if ( errorCode != 0 ) {
                    std::stringstream message;
                    message << "Failed to write all bytes because of: " << strerror( errorCode )
                            << " (" << errorCode << ")";
                    throw std::runtime_error( std::move( message ).str() );
                }

                if ( outputBuffer != nullptr ) {
                    using rapidgzip::deflate::DecodedData;

                    size_t nBytesCopied{ 0 };
                    for ( auto it = DecodedData::Iterator( *chunkData, offsetInBlock, dataToWriteSize );
                          static_cast<bool>( it ); ++it )
                    {
                        const auto& [buffer, bufferSize] = *it;
                        auto* const currentBufferPosition = outputBuffer + nBytesDecoded + nBytesCopied;
                        std::memcpy( currentBufferPosition, buffer, bufferSize );
                        nBytesCopied += bufferSize;
                    }
                }

                nBytesDecoded += dataToWriteSize;
            };

        if ( ( outputFileDescriptor == -1 ) && ( outputBuffer == nullptr ) ) {
            /* An empty std::function gives that read method options to optimize, e.g., via seeking. */
            return read( WriteFunctor{}, nBytesToRead );
        }
        return read( writeFunctor, nBytesToRead );
    }

    size_t
    read( const WriteFunctor& writeFunctor,
          const size_t        nBytesToRead = std::numeric_limits<size_t>::max() )
    {
        if ( !writeFunctor && m_blockMap->finalized() ) {
            const auto oldOffset = tell();
            const auto newOffset = seek( nBytesToRead > static_cast<size_t>( std::numeric_limits<long long int>::max() )
                                         ? std::numeric_limits<long long int>::max()
                                         : nBytesToRead,
                                         SEEK_CUR );
            return newOffset - oldOffset;
        }

        if ( closed() ) {
            throw std::invalid_argument( "You may not call read on closed ParallelGzipReader!" );
        }

        if ( eof() || ( nBytesToRead == 0 ) ) {
            return 0;
        }

        size_t nBytesDecoded = 0;
        while ( ( nBytesDecoded < nBytesToRead ) && !eof() ) {
        #ifdef WITH_PYTHON_SUPPORT
            checkPythonSignalHandlers();
            const ScopedGILUnlock unlockedGIL;
        #endif

            const auto blockResult = chunkFetcher().get( m_currentPosition );
            if ( !blockResult ) {
                m_atEndOfFile = true;
                break;
            }
            const auto& [decodedOffsetInBytes, chunkData] = *blockResult;

            if ( chunkData->containsMarkers() ) {
                throw std::logic_error( "Did not expect to get results with markers!" );
            }

            /* Copy data from fetched block to output. */

            const auto offsetInBlock = m_currentPosition - decodedOffsetInBytes;
            const auto blockSize = chunkData->decodedSizeInBytes;
            if ( offsetInBlock >= blockSize ) {
                std::stringstream message;
                message << "[ParallelGzipReader] Block does not contain the requested offset! "
                        << "Requested offset from chunk fetcher: " << m_currentPosition
                        << " (" << formatBytes( m_currentPosition ) << ")"
                        << ", decoded offset: " << decodedOffsetInBytes
                        << " (" << formatBytes( decodedOffsetInBytes ) << ")"
                        << ", block data encoded offset: " << formatBits( chunkData->encodedOffsetInBits )
                        << ", block data encoded size: " << formatBits( chunkData->encodedSizeInBits )
                        << ", block data size: " << chunkData->decodedSizeInBytes
                        << " (" << formatBytes( chunkData->decodedSizeInBytes ) << ")"
                        << " markers: " << chunkData->dataWithMarkersSize();
                throw std::logic_error( std::move( message ).str() );
            }

            const auto nBytesToDecode = std::min( blockSize - offsetInBlock, nBytesToRead - nBytesDecoded );

            [[maybe_unused]] const auto tCRC32Start = now();
            processCRC32( chunkData, offsetInBlock, nBytesToDecode );
            if ( m_statisticsEnabled ) {
                m_crc32Time += duration( tCRC32Start );
            }

            if ( writeFunctor ) {
                [[maybe_unused]] const auto tWriteStart = now();
                writeFunctor( chunkData, offsetInBlock, nBytesToDecode );
                if ( m_statisticsEnabled ) {
                    m_writeOutputTime += duration( tWriteStart );
                }
            }

            nBytesDecoded += nBytesToDecode;
            m_currentPosition += nBytesToDecode;

            const auto& [lock, file] = m_sharedFileReader->underlyingFile();
            // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
            auto* const singlePassFileReader = dynamic_cast<SinglePassFileReader*>( file );
            if ( singlePassFileReader != nullptr ) {
                /* Release only up to the beginning of the currently used chunk in order to theoretically enable
                 * to clear the full cache and then continue again. This effectively require a recomputation of
                 * the current chunk if it was not fully read yet. */
                singlePassFileReader->releaseUpTo( /* floor int division */ chunkData->encodedOffsetInBits / CHAR_BIT );
            }

            if ( !m_keepIndex && m_windowMap ) {
                m_windowMap->releaseUpTo( chunkData->encodedOffsetInBits );
            }
        }

        return nBytesDecoded;
    }

    size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) override
    {
        if ( closed() ) {
            throw std::invalid_argument( "You may not call seek on closed ParallelGzipReader!" );
        }

        if ( origin == SEEK_END ) {
            /* size() requires the block offsets to be available! */
            if ( !m_blockMap->finalized() ) {
                read();
            }
        }
        const auto positiveOffset = effectiveOffset( offset, origin );

        if ( positiveOffset == tell() ) {
            /* This extra check for EOF is necessary for empty files! */
            m_atEndOfFile = m_blockMap->finalized() && ( m_currentPosition >= m_blockMap->back().second );
            return positiveOffset;
        }

        /* Backward seeking is no problem at all! 'tell' may only return <= size()
         * as value meaning we are now < size() and therefore EOF can be cleared! */
        if ( positiveOffset < tell() ) {
            if ( !m_keepIndex ) {
                throw std::invalid_argument( "Seeking (back) not supported when index-keeping has been disabled!" );
            }

            if ( !seekable() ) {
                throw std::invalid_argument( "Cannot seek backwards with non-seekable input!" );
            }
            m_atEndOfFile = false;
            m_currentPosition = positiveOffset;
            return positiveOffset;
        }

        /* m_blockMap is only accessed by read and seek, which are not to be called from different threads,
         * so we do not have to lock it. */
        const auto blockInfo = m_blockMap->findDataOffset( positiveOffset );
        if ( positiveOffset < blockInfo.decodedOffsetInBytes ) {
            throw std::logic_error( "Block map returned unwanted block!" );
        }

        if ( blockInfo.contains( positiveOffset ) ) {
            m_currentPosition = positiveOffset;
            m_atEndOfFile = m_blockMap->finalized() && ( m_currentPosition >= m_blockMap->back().second );
            return tell();
        }

        if ( m_blockMap->finalized() ) {
            m_atEndOfFile = true;
            m_currentPosition = m_blockMap->back().second;
            return tell();
        }

        /* Jump to furthest known point as performance optimization. Note that even if that is right after
         * the last byte, i.e., offset == size(), then no eofbit is set even in ifstream! In ifstream you
         * can even seek to after the file end with no fail bits being set in my tests! */
        m_atEndOfFile = false;
        m_currentPosition = blockInfo.decodedOffsetInBytes + blockInfo.decodedSizeInBytes;
        read( -1, nullptr, positiveOffset - tell() );
        return tell();
    }

    // NOLINTEND(misc-no-recursion)

    /* Block compression specific methods */

    [[nodiscard]] bool
    blockOffsetsComplete() const
    {
        return m_blockMap->finalized();
    }

    /**
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    [[nodiscard]] std::map<size_t, size_t>
    blockOffsets()
    {
        if ( !m_blockMap->finalized() ) {
            read();
            if ( !m_blockMap->finalized() || !blockFinder().finalized() ) {
                throw std::logic_error( "Reading everything should have finalized the block map!" );
            }
        }

        return m_blockMap->blockOffsets();
    }

    /**
     * This is the first instance for me where returning a const value makes sense because it contains
     * a shared pointer to the WindowMap, which is not to be modified. Making GzipIndex const forces
     * the caller to deep clone the index and WindowMap for, e.g., the setBlockOffsets API, which
     * destructively moves from the WindowMap.
     */
    [[nodiscard]] const GzipIndex  // NOLINT(readability-const-return-type)
    gzipIndex( bool withLineOffsets = false )
    {
        const auto offsets = blockOffsets();  // Also finalizes reading implicitly.
        if ( offsets.empty() || !m_windowMap ) {
            return {};
        }

        const auto archiveSize = m_sharedFileReader->size();
        if ( !archiveSize && !m_indexIsImported ) {
            /* If the index was imported, then this warning is moot and using the last chunk offset is sufficient. */
            std::cerr << "[Warning] The input file size should have become available after finalizing the index!\n";
            std::cerr << "[Warning] Will use the last chunk end offset as size. This might lead to errors on import!\n";
        }

        GzipIndex index;
        index.compressedSizeInBytes = archiveSize ? *archiveSize : ceilDiv( offsets.rbegin()->first, 8U );
        index.uncompressedSizeInBytes = offsets.rbegin()->second;
        index.windowSizeInBytes = 32_Ki;

        if ( withLineOffsets ) {
            if ( !m_newlineCharacter ) {
                throw std::runtime_error( "Cannot add line offsets to index when they were not gathered!" );
            }
            index.hasLineOffsets = true;

            switch ( m_newlineCharacter.value() )
            {
            case '\n':
                index.newlineFormat = NewlineFormat::LINE_FEED;
                break;
            case '\r':
                index.newlineFormat = NewlineFormat::CARRIAGE_RETURN;
                break;
            default:
                throw std::runtime_error(
                    "Cannot add line offsets to index when the gathered line offsets gathered are something "
                    "other than newline or carriage return!" );
            }
        }

        /* Heuristically determine a checkpoint spacing from the existing checkpoints. */
        size_t maximumDecompressedSpacing{ 0 };
        for ( auto it = offsets.begin(), nit = std::next( offsets.begin() ); nit != offsets.end(); ++it, ++nit ) {
            maximumDecompressedSpacing = std::max( maximumDecompressedSpacing, nit->second - it->second );
        }
        index.checkpointSpacing = maximumDecompressedSpacing / 32_Ki * 32_Ki;

        auto lineOffset = m_newlineOffsets.begin();
        for ( const auto& [compressedOffsetInBits, uncompressedOffsetInBytes] : offsets ) {
            Checkpoint checkpoint;
            checkpoint.compressedOffsetInBits = compressedOffsetInBits;
            checkpoint.uncompressedOffsetInBytes = uncompressedOffsetInBytes;

            if ( index.hasLineOffsets ) {
                while ( ( lineOffset != m_newlineOffsets.end() )
                        && ( lineOffset->uncompressedOffsetInBytes < uncompressedOffsetInBytes ) )
                {
                    ++lineOffset;
                }

                if ( lineOffset == m_newlineOffsets.end() ) {
                    std::stringstream message;
                    message << "Failed to find line offset for uncompressed offset: "
                            << formatBytes( uncompressedOffsetInBytes ) << ", number of line offsets to stored: "
                            << m_newlineOffsets.size();
                    throw std::runtime_error( std::move( message ).str() );
                }

                if ( lineOffset->uncompressedOffsetInBytes != uncompressedOffsetInBytes ) {
                    throw std::logic_error( "Line offset not found for uncompressed offset "
                                            + std::to_string( uncompressedOffsetInBytes ) + "!" );
                }

                checkpoint.lineOffset = lineOffset->lineOffset;
            }

            index.checkpoints.emplace_back( checkpoint );
        }

        index.windows = m_windowMap;

        return index;
    }

    /**
     * Same as @ref blockOffsets but it won't force calculation of all blocks and simply returns
     * what is available at call time.
     * @return vectors of block data: offset in file, offset in decoded data
     *         (cumulative size of all prior decoded blocks).
     */
    [[nodiscard]] std::map<size_t, size_t>
    availableBlockOffsets() const
    {
        return m_blockMap->blockOffsets();
    }

    [[nodiscard]] auto
    statistics() const
    {
        if ( !m_chunkFetcher ) {
            throw std::invalid_argument( "No chunk fetcher initialized!" );
        }
        return m_chunkFetcher->statistics();
    }

    void
    setCRC32Enabled( bool enabled )
    {
        if ( m_crc32.enabled() == enabled ) {
            return;
        }

        m_crc32.setEnabled( enabled && ( tell() == 0 ) );
        applyChunkDataConfiguration();
    }

    void
    setMaxDecompressedChunkSize( uint64_t maxDecompressedChunkSize )
    {
        /* Anything smaller than the chunk size doesn't make much sense. Even that would be questionable.
         * as it would lead to slow downs in almost every case. */
        m_chunkConfiguration.maxDecompressedChunkSize = std::max( m_chunkSizeInBytes, maxDecompressedChunkSize );
        applyChunkDataConfiguration();
    }

    [[nodiscard]] uint64_t
    maxDecompressedChunkSize() const noexcept
    {
        return m_chunkConfiguration.maxDecompressedChunkSize;
    }

    /**
     * Must only be changed before doing any read call! Else, some of the chunks will already have been processed
     * with the existing newline format.
     *
     * @param newlineCharacter If nullopt, newline counting will be disabled.
     */
    void
    setNewlineCharacter( const std::optional<char>& newlineCharacter )
    {
        if ( newlineCharacter == m_newlineCharacter ) {
            return;
        }

        /* The check could be improved here, e.g., check for queued futures. */
        if ( !m_newlineOffsets.empty() || ( m_blockMap && !m_blockMap->empty() ) ) {
            throw std::invalid_argument( "May not change newline counting behavior after some chunks have been read!" );
        }
        m_newlineCharacter = newlineCharacter;
        applyChunkDataConfiguration();
    }

    [[nodiscard]] std::optional<char>
    newlineCharacter() const noexcept
    {
        return m_newlineCharacter;
    }

    [[nodiscard]] const std::vector<NewlineOffset>&
    newlineOffsets() const noexcept
    {
        return m_newlineOffsets;
    }

private:
    void
    setBlockOffsets( const std::map<size_t, size_t>& offsets )
    {
        /**
         * @todo Join very small consecutive block offsets until it roughly reflects the chunk size?
         * Because currently, the version with the BGZI index is slower than without!
         * rapidgzip -d -o /dev/null 4GiB-base64.bgz
         * > Decompressed in total 4294967296 B in 0.565016 s -> 7601.49 MB/s
         * rapidgzip -d -o /dev/null --import-index 4GiB-base64.bgz
         * > Decompressed in total 4294901760 B in 1.22275 s -> 3512.5 MB/s
         */

        if ( offsets.empty() ) {
            if ( m_blockMap->dataBlockCount() == 0 ) {
                return;
            }
            throw std::invalid_argument( "May not clear offsets. Construct a new ParallelGzipReader instead!" );
        }

        setBlockFinderOffsets( offsets );

        if ( offsets.size() < 2 ) {
            throw std::invalid_argument( "Block offset map must contain at least one valid block and one EOS block!" );
        }
        m_blockMap->setBlockOffsets( offsets );
    }

public:
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

        m_indexIsImported = true;
        m_keepIndex = true;

        std::optional<char> newlineCharacter;
        switch ( index.newlineFormat )
        {
        case NewlineFormat::LINE_FEED:
            newlineCharacter = '\n';
            break;
        case NewlineFormat::CARRIAGE_RETURN:
            newlineCharacter = '\r';
            break;
        }

        if ( index.hasLineOffsets && newlineCharacter ) {
            m_newlineCharacter = newlineCharacter;
            m_newlineOffsets.resize( index.checkpoints.size() );
            std::transform( index.checkpoints.begin(), index.checkpoints.end(), m_newlineOffsets.begin(),
                            [] ( const auto& checkpoint ) {
                NewlineOffset newlineOffset;
                newlineOffset.lineOffset = checkpoint.lineOffset;
                newlineOffset.uncompressedOffsetInBytes = checkpoint.uncompressedOffsetInBytes;
                return newlineOffset;
            } );

            /* Checkpoints should already be sorted and therefore also the newline offsets. Check just to be sure.
             * We are not sorting here because it may be impossible to sort by line offsets and uncompressed offsets
             * for inconsistent index data! */
            const auto lessLineOffset = [] ( const auto& a, const auto& b ) { return a.lineOffset < b.lineOffset; };
            if ( !std::is_sorted( m_newlineOffsets.begin(), m_newlineOffsets.end(), lessLineOffset ) ) {
                throw std::invalid_argument( "Index checkpoints must be sorted by line offsets!" );
            }
        }

        /* Generate simple compressed to uncompressed offset map from index. */
        std::map<size_t, size_t> newBlockOffsets;
        for ( size_t i = 0; i < index.checkpoints.size(); ++i ) {
            const auto& checkpoint = index.checkpoints[i];

            /**
             * Skip emission of an index, if the next checkpoint would still let us be below the chunk size.
             * Always copy the zeroth index as is necessary for a valid index!
             *
             * This is for merging very small index points as might happen when importing BGZF indexes.
             * Small index will lead to relatively larger overhead for the threading and will degrade performance:
             * Merge n blocks:
             *     0 -> ~3.3 GB/s, Total existing blocks: 65793
             *     2 -> ~5.0 GB/s, Total existing blocks: 32897
             *     4 -> ~5.7 GB/s, Total existing blocks: 16449
             *     8 -> ~6.0 GB/s, Total existing blocks: 8225
             *    16 -> ~6.8 GB/s, Total existing blocks: 4113
             *    32 -> ~6.8 GB/s, Total existing blocks: 2057
             *    64 -> ~7.2 GB/s, Total existing blocks: 1029
             *   128 -> ~6.9 GB/s, Total existing blocks: 515
             *
             * Without index import (chunk size 4 MiB):
             * src/tools/rapidgzip -d -o /dev/null 4GiB-base64.bgz
             *   Total existing blocks: 766 blocks
             *   Index reading took: 0.00259098 s
             *
             *   Decompressed in total 4294967296 B in 0.580731 s -> 7395.79 MB/s
             *   Decompressed in total 4294967296 B in 0.576022 s -> 7456.26 MB/s
             *   Decompressed in total 4294967296 B in 0.597594 s -> 7187.1 MB/s
             */
            if ( !newBlockOffsets.empty() && ( i + 1 < index.checkpoints.size() )
                 && ( index.checkpoints[i + 1].uncompressedOffsetInBytes - newBlockOffsets.rbegin()->second
                      <= m_chunkSizeInBytes ) )
            {
                continue;
            }

            newBlockOffsets.emplace( checkpoint.compressedOffsetInBits, checkpoint.uncompressedOffsetInBytes );

            /* Copy window data.
             * For some reason, indexed_gzip also stores windows for the very last checkpoint at the end of the file,
             * which is useless because there is nothing thereafter. But, don't filter it here so that exportIndex
             * mirrors importIndex better.
             * Bgzip indexes do not have windows because they are not needed, so we do not have to insert anything
             * into the WindowMap in that case. Bgzip indexes will be detected by the magic bytes and in that case
             * windows should never be looked up in the WindowMap in the first place. */
            if ( const auto window = lockedWindows.second->find( checkpoint.compressedOffsetInBits );
                 window != lockedWindows.second->end() )
            {
                m_windowMap->emplaceShared( checkpoint.compressedOffsetInBits, window->second );
            }
        }

        /* Input file-end offset if not included in checkpoints. */
        if ( const auto fileEndOffset = newBlockOffsets.find( index.compressedSizeInBytes * 8 );
             fileEndOffset == newBlockOffsets.end() )
        {
            newBlockOffsets.emplace( index.compressedSizeInBytes * 8, index.uncompressedSizeInBytes );
            m_windowMap->emplace( index.compressedSizeInBytes * 8, {}, CompressionType::NONE );
        } else if ( fileEndOffset->second != index.uncompressedSizeInBytes ) {
            throw std::invalid_argument( "Index has contradicting information for the file end information!" );
        }

        setBlockOffsets( newBlockOffsets );

        chunkFetcher().clearCache();
    }

    void
    importIndex( UniqueFileReader indexFile )
    {
        const auto t0 = now();
        setBlockOffsets( readGzipIndex( std::move( indexFile ), m_sharedFileReader->clone(),
                                        m_fetcherParallelization ) );
        if ( m_showProfileOnDestruction ) {
            std::cerr << "[ParallelGzipReader::importIndex] Took " << duration( t0 ) << " s\n";
        }
    }

    void
    exportIndex( const std::function<void( const void* buffer, size_t size )>& checkedWrite,
                 const IndexFormat                                             indexFormat = IndexFormat::INDEXED_GZIP )
    {
        const auto t0 = now();

        if ( !m_keepIndex ) {
            throw std::invalid_argument( "Exporting index not supported when index-keeping has been disabled!" );
        }

        switch ( indexFormat )
        {
        case IndexFormat::INDEXED_GZIP:
            indexed_gzip::writeGzipIndex( gzipIndex(), checkedWrite );
            break;
        case IndexFormat::GZTOOL:
            gztool::writeGzipIndex( gzipIndex( false ), checkedWrite );
            break;
        case IndexFormat::GZTOOL_WITH_LINES:
            gztool::writeGzipIndex( gzipIndex( true ), checkedWrite );
            break;
        }

        if ( m_showProfileOnDestruction ) {
            std::cerr << "[ParallelGzipReader::exportIndex] Took " << duration( t0 ) << " s\n";
        }
    }

#ifdef WITH_PYTHON_SUPPORT
    void
    importIndex( PyObject* pythonObject )
    {
        importIndex( std::make_unique<PythonFileReader>( pythonObject ) );
    }

    void
    exportIndex( PyObject*         pythonObject,
                 const IndexFormat indexFormat = IndexFormat::INDEXED_GZIP )
    {
        const auto file = std::make_unique<PythonFileReader>( pythonObject );
        const auto checkedWrite =
            [&file] ( const void* buffer, size_t size )
            {
                if ( file->write( reinterpret_cast<const char*>( buffer ), size ) != size ) {
                    throw std::runtime_error( "Failed to write data to index!" );
                }
            };

        exportIndex( checkedWrite, indexFormat );
    }
#endif

    void
    gatherLineOffsets()
    {
        /* Check whether the newline information has already been collected from an imported index or earlier call. */
        if ( !m_newlineCharacter ) {
            return;
        }

        const Finally restorePosition{ [this, oldOffset = tell()] () { seek( oldOffset ); } };

        /* If it was already toggled on, simply read until the end to gather all offsets. */
        if ( !blockOffsetsComplete() ) {
            read();
            return;
        }

        /* Block offset is complete, check if line offsets are complete by checking the last one. */
        uint64_t processedBytes{ m_newlineOffsets.empty() ? 0 : m_newlineOffsets.back().uncompressedOffsetInBytes };
        if ( const auto fileSize = size(); fileSize && !m_newlineOffsets.empty() && ( processedBytes >= fileSize ) ) {
            return;
        }

        /* This may be necessary when the block map has been finalized because an index without line information
         * has been imported! In that case, we need to gather line information manually like a user would. */

        /* Collect line offsets until the next chunk offset has been added to the map. Then, we can look for the
         * line number at that exact chunk offset and insert it and clear our temporary results. */
        uint64_t processedLines{ m_newlineOffsets.empty() ? 0 : m_newlineOffsets.back().lineOffset };
        /** Index i stores the byte offset for the (processedLines + i)-th line. */
        std::vector<uint64_t> newlineOffsets;

        const auto collectLineOffsets =
            [this, &processedLines, &newlineOffsets, &processedBytes, newlineCharacter = m_newlineCharacter.value()]
            ( const std::shared_ptr<rapidgzip::ChunkData>& chunkData,
              const size_t                                 offsetInChunk,
              const size_t                                 dataToWriteSize )
            {
                /* Iterate over the requested data range of the chunk and collect byte offsets for every newline. */
                using rapidgzip::deflate::DecodedData;
                for ( auto it = DecodedData::Iterator( *chunkData, offsetInChunk, dataToWriteSize );
                      static_cast<bool>( it ); ++it )
                {
                    const auto& [buffer, size] = *it;

                    const std::string_view view{ reinterpret_cast<const char*>( buffer ), size };
                    for ( auto position = view.find( newlineCharacter, 0 );
                          position != std::string_view::npos;
                          position = view.find( newlineCharacter, position + 1 ) )
                    {
                        newlineOffsets.emplace_back( processedBytes + position );
                    }

                    processedBytes += size;
                }

                /* Iterate over all found newline offsets and start inserting an actual byte->newline offset pair
                 * but only once per chunk to reduce the index size. */
                auto it = newlineOffsets.begin();
                while ( it != newlineOffsets.end() ) {
                    const auto chunkInfo = m_blockMap->findDataOffset( *it );
                    if ( !chunkInfo.contains( *it ) ) {
                        /* I don't think this can happen. It happens when the currently processed chunk
-                        * is not yet registered in the chunk map. */
                        std::cerr << "[Warning] Offset in processed chunk was not found in chunk map!\n";
                        break;
                    }

                    if ( m_newlineOffsets.empty() || ( m_newlineOffsets.back().uncompressedOffsetInBytes != *it ) ) {
                        NewlineOffset newlineOffset;
                        newlineOffset.lineOffset = static_cast<uint64_t>( std::distance( newlineOffsets.begin(), it ) )
                                                   + processedLines;
                        newlineOffset.uncompressedOffsetInBytes = chunkInfo.decodedOffsetInBytes;

                        if ( !m_newlineOffsets.empty() ) {
                            if ( m_newlineOffsets.back().uncompressedOffsetInBytes >= *it ) {
                                std::stringstream message;
                                message << "Got earlier or equal chunk offset than the last processed one! "
                                        << "Last newline byte offset: "
                                        << m_newlineOffsets.back().uncompressedOffsetInBytes
                                        << ", found newline byte offset: " << *it;
                                throw std::logic_error( std::move( message ).str() );
                            }

                            if ( m_newlineOffsets.back().lineOffset > newlineOffset.lineOffset ) {
                                throw std::logic_error( "Got earlier line offset than the last processed one!" );
                            }
                        }

                        m_newlineOffsets.emplace_back( newlineOffset );
                    }

                    /* Skip over all newlines still in the last processed chunk. */
                    while ( ( it != newlineOffsets.end() ) && chunkInfo.contains( *it ) ) {
                        ++it;
                    }
                }

                processedLines += static_cast<uint64_t>( std::distance( newlineOffsets.begin(), it ) );
                newlineOffsets.erase( newlineOffsets.begin(), it );
            };

        seekTo( processedBytes );
        read( collectLineOffsets );

        /* Insert information for the end-of-file offset. */
        if ( m_newlineOffsets.empty() || ( processedBytes > m_newlineOffsets.back().uncompressedOffsetInBytes ) ) {
            NewlineOffset newlineOffset;
            newlineOffset.uncompressedOffsetInBytes = processedBytes;
            newlineOffset.lineOffset = processedLines + newlineOffsets.size();
            m_newlineOffsets.emplace_back( newlineOffset );
        }
    }

    /**
     * @return number of processed bits of compressed bzip2 input file stream
     * @note Bzip2 is block based and blocks are currently read fully, meaning that the granularity
     *       of the returned position is ~100-900kB. It's only useful for a rough estimate.
     */
    [[nodiscard]] size_t
    tellCompressed() const
    {
        if ( !m_blockMap || m_blockMap->empty() ) {
            return 0;
        }

        const auto blockInfo = m_blockMap->findDataOffset( m_currentPosition );
        if ( blockInfo.contains( m_currentPosition ) ) {
            return blockInfo.encodedOffsetInBits;
        }
        return m_blockMap->back().first;
    }

    /**
     * Closes all threads and saves the work. They will be restarted when needed again, e.g., on seek or read.
     * This is intended for use with fusepy. You can start a ParallelGzipReader use it to create the block map
     * and print out user output and then you join all threads before FUSE forks the process. FUSE requires
     * threads to be created after it forks, it seems:
     * @see https://github.com/libfuse/libfuse/wiki/FAQ#how-should-threads-be-started
     * Personally, the only problem I observed was background process not finishing even after unmounting,
     * however, contrary to the FAQ it seems that threads were not joined because the file system seemed to work.
     */
    void
    joinThreads()
    {
        m_chunkFetcher.reset();
        m_blockFinder.reset();
    }

    /**
     * Index-keeping can be disabled as a memory usage optimization when it will never be needed.
     * Currently, this will clear windows for chunks that have been fully decompressed once.
     * Trying to seek in the file with this option enabled will throw an error!
     */
    void
    setKeepIndex( bool keep )
    {
        m_keepIndex = keep;
        applyChunkDataConfiguration();
    }

    void
    setWindowSparsity( bool useSparseWindows )
    {
        m_windowSparsity = useSparseWindows;
        applyChunkDataConfiguration();
    }

    void
    setWindowCompressionType( CompressionType windowCompressionType )
    {
        m_windowCompressionType = windowCompressionType;
        applyChunkDataConfiguration();
    }

    [[nodiscard]] std::string
    fileTypeAsString()
    {
        return toString( blockFinder().fileType() );
    }

    void
    setDeflateStreamCRC32s( std::unordered_map<size_t, uint32_t> crc32s )
    {
        m_deflateStreamCRC32s = std::move( crc32s );
    }

    void
    addDeflateStreamCRC32( size_t   endOfStreamOffsetInBytes,
                           uint32_t crc32 )
    {
        m_deflateStreamCRC32s.insert_or_assign( endOfStreamOffsetInBytes, crc32 );
    }

private:
    void
    applyChunkDataConfiguration()
    {
        if ( !m_chunkFetcher ) {
            return;
        }

        m_chunkConfiguration.crc32Enabled = m_crc32.enabled();
        m_chunkConfiguration.windowCompressionType =
            m_keepIndex ? m_windowCompressionType : std::make_optional( CompressionType::NONE );
        /* Window sparsity only makes sense when keeping the index. */
        m_chunkConfiguration.windowSparsity = m_keepIndex && m_windowSparsity;
        m_chunkConfiguration.newlineCharacter = m_newlineCharacter;

        m_chunkFetcher->setChunkConfiguration( m_chunkConfiguration );
    }

    BlockFinder&
    blockFinder()  // NOLINT(misc-no-recursion)
    {
        /* This guard makes the warned-about recursion via setBlockFinderOffsets safe. */
        if ( m_blockFinder ) {
            return *m_blockFinder;
        }

        if ( !m_startBlockFinder ) {
            throw std::logic_error( "Block finder creator was not initialized correctly!" );
        }

        m_blockFinder = m_startBlockFinder();
        if ( !m_blockFinder ) {
            throw std::logic_error( "Block finder creator failed to create new block finder!" );
        }

        if ( m_blockMap->finalized() ) {
            setBlockFinderOffsets( m_blockMap->blockOffsets() );
        }

        return *m_blockFinder;
    }

    ChunkFetcher&
    chunkFetcher()
    {
        if ( m_chunkFetcher ) {
            return *m_chunkFetcher;
        }

        /* As a side effect, blockFinder() creates m_blockFinder if not already initialized! */
        blockFinder();

        m_chunkFetcher = std::make_unique<ChunkFetcher>( ensureSharedFileReader( m_sharedFileReader->clone() ),
                                                         m_blockFinder, m_blockMap, m_windowMap,
                                                         m_fetcherParallelization );
        if ( !m_chunkFetcher ) {
            throw std::logic_error( "Block fetcher should have been initialized!" );
        }

        m_chunkFetcher->setShowProfileOnDestruction( m_showProfileOnDestruction );
        m_chunkFetcher->setStatisticsEnabled( m_statisticsEnabled );
        m_chunkFetcher->addChunkIndexingCallback(
            [this] ( const auto& chunk, auto ) { this->gatherLineOffsets( chunk ); } );
        applyChunkDataConfiguration();

        return *m_chunkFetcher;
    }

    void
    gatherLineOffsets( const std::shared_ptr<const ChunkData>& chunk )
    {
        if ( !m_newlineCharacter.has_value() ) {
            return;
        }

        if ( !chunk ) {
            throw std::logic_error( "ParallelGzipReader::gatherLineOffsets should only be called with valid chunk!" );
        }

        for ( const auto& subchunk : chunk->subchunks() ) {
            if ( !subchunk.newlineCount ) {
                throw std::logic_error( "Newline count in subchunk is missing!" );
            }
            if ( chunk->configuration.newlineCharacter != m_newlineCharacter ) {
                throw std::logic_error( "Newline character in subchunk does not match the configured one!" );
            }

            const auto blockInfo = m_blockMap->getEncodedOffset( subchunk.encodedOffset );
            if ( !blockInfo ) {
                std::stringstream message;
                message << "Failed to find subchunk offset: " << formatBits( subchunk.encodedOffset )
                        << "even though it should have been inserted at the top of this method!";
                throw std::logic_error( std::move( message ).str() );
            }

            if ( m_newlineOffsets.empty() ) {
                m_newlineOffsets.emplace_back( NewlineOffset{ 0, 0 } );
            }

            const auto& lastLineCount = m_newlineOffsets.back();
            if ( lastLineCount.uncompressedOffsetInBytes != blockInfo->decodedOffsetInBytes ) {
                std::stringstream message;
                message << "Did not find line count for preceding decompressed offset: "
                        << formatBytes( blockInfo->decodedOffsetInBytes );
                throw std::logic_error( std::move( message ).str() );
            }

            m_newlineOffsets.emplace_back( NewlineOffset{
                lastLineCount.lineOffset + subchunk.newlineCount.value(),
                blockInfo->decodedOffsetInBytes + subchunk.decodedSize
            } );
        }
    }

    void
    setBlockFinderOffsets( const std::map<size_t, size_t>& offsets )  // NOLINT(misc-no-recursion)
    {
        if ( offsets.empty() ) {
            throw std::invalid_argument( "A non-empty list of block offsets is required!" );
        }

        typename BlockFinder::BlockOffsets encodedBlockOffsets;
        for ( auto it = offsets.begin(), nit = std::next( offsets.begin() ); nit != offsets.end(); ++it, ++nit )
        {
            /* Ignore blocks with no data, i.e., EOS blocks. */
            if ( it->second != nit->second ) {
                encodedBlockOffsets.push_back( it->first );
            }
        }
        /* The last block is not pushed because "std::next( it )" is end but last block must be EOS anyways
         * or else BlockMap will not work correctly because the implied size of that last block is 0! */

        blockFinder().setBlockOffsets( std::move( encodedBlockOffsets ) );
    }

    void
    processCRC32( const std::shared_ptr<ChunkData>& chunkData,
                  [[maybe_unused]] size_t const     offsetInBlock,
                  [[maybe_unused]] size_t const     dataToWriteSize )
    {
        if ( ( m_nextCRC32ChunkOffset == 0 ) && m_blockFinder ) {
            const auto [offset, errorCode] = m_blockFinder->get( /* block index */ 0, /* timeout */ 0 );
            if ( offset && ( errorCode == BlockFinder::GetReturnCode::SUCCESS ) ) {
                m_nextCRC32ChunkOffset = *offset;
            }
        }

        if ( !m_crc32.enabled()
             || ( m_nextCRC32ChunkOffset != chunkData->encodedOffsetInBits )
             || chunkData->crc32s.empty() ) {
            return;
        }

        m_nextCRC32ChunkOffset = chunkData->encodedOffsetInBits + chunkData->encodedSizeInBits;

        /* As long as CRC32 is enabled, this should not happen and we filter above for !m_crc32.enabled(). */
        if ( chunkData->crc32s.size() != chunkData->footers.size() + 1 ) {
            throw std::logic_error( "Fewer CRC32s in chunk than expected based on the gzip footers!" );
        }

        const auto totalCRC32StreamSize = std::accumulate(
            chunkData->crc32s.begin(), chunkData->crc32s.end(), size_t( 0 ),
            [] ( size_t sum, const auto& calculator ) { return sum + calculator.streamSize(); } );
        if ( totalCRC32StreamSize != chunkData->decodedSizeInBytes ) {
            std::stringstream message;
            message << "CRC32 computation stream size (" << formatBytes( totalCRC32StreamSize ) << ") differs from "
                    << "chunk size: " << formatBytes( chunkData->decodedSizeInBytes ) << "!\n"
                    << "Please open an issue or disable integrated CRC32 verification as a quick workaround.";
            throw std::logic_error( std::move( message ).str() );
        }

        /* Process CRC32 of chunk. */
        m_crc32.append( chunkData->crc32s.front() );
        for ( size_t i = 0; i < chunkData->footers.size(); ++i ) {
            const auto& footer = chunkData->footers[i];
            const auto footerByteOffset = ceilDiv( footer.blockBoundary.encodedOffset, CHAR_BIT );
            if ( const auto externalCRC32 = m_deflateStreamCRC32s.find( footerByteOffset );
                 externalCRC32 != m_deflateStreamCRC32s.end() )
            {
                m_crc32.verify( m_crc32.crc32() );
            } else if ( hasCRC32( chunkData->configuration.fileType ) && m_crc32.verify( footer.gzipFooter.crc32 ) ) {
                m_verifiedCRC32Count++;
            }
            m_crc32 = chunkData->crc32s.at( i + 1 );
        }
    }

private:
    uint64_t m_chunkSizeInBytes{ 4_Mi };
    ChunkConfiguration m_chunkConfiguration;

    std::unique_ptr<SharedFileReader> m_sharedFileReader;

    size_t m_currentPosition = 0;  /**< the current position as can only be modified with read or seek calls. */
    bool m_atEndOfFile = false;

    /** Benchmarking */
    bool m_statisticsEnabled{ false };
    bool m_showProfileOnDestruction{ false };
    double m_writeOutputTime{ 0 };
    double m_crc32Time{ 0 };
    uint64_t m_verifiedCRC32Count{ 0 };

    size_t const m_fetcherParallelization;

    std::function<std::shared_ptr<BlockFinder>( void )> const m_startBlockFinder;

    /** Necessary for prefetching decoded blocks in parallel. */
    std::shared_ptr<BlockFinder>     m_blockFinder;
    std::shared_ptr<BlockMap>  const m_blockMap{ std::make_shared<BlockMap>() };
    /**
     * The window map should contain windows to all encoded block offsets inside @ref m_blockMap.
     * The windows are stored in a separate map even though all keys should be identical because BlockMap is
     * too "finished". I don't see how to generically and readably add generic user data / windows to it.
     * Furthermore, the windows might potentially be written out-of-order while block offsets should be inserted
     * in order into @ref m_blockMap.
     */
    std::shared_ptr<WindowMap> const m_windowMap{ std::make_shared<WindowMap>() };
    bool m_keepIndex{ true };
    bool m_windowSparsity{ true };
    std::optional<CompressionType> m_windowCompressionType;
    std::unique_ptr<ChunkFetcher> m_chunkFetcher;
    /**
     * Note that the uncompressed offset can point to any byte offset inside the line depending on how the chunks
     * are split. Only the offset to the 0-th line is exact of course. To get any other line beginning exactly,
     * you need to start from the previous line and search for the next newline.
     * Note also that not all line offsets have to be in this vector. That's why it is a vector of pairs and not
     * a simply vector of values. Line offsets are only available at spacings. To get an exact line offset, you
     * need to start reading from the next smaller one and skip over as many newline characters as necessary.
     */
    std::vector<NewlineOffset> m_newlineOffsets;
    std::optional<char> m_newlineCharacter;

    CRC32Calculator m_crc32;
    uint64_t m_nextCRC32ChunkOffset{ 0 };
    std::unordered_map<size_t, uint32_t> m_deflateStreamCRC32s;

    bool m_indexIsImported{ false };
};
}  // namespace rapidgzip
