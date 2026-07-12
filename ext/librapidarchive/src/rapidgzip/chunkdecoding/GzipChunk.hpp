#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <zlib.h>

#include <core/common.hpp>
#include <rapidgzip/blockfinder/DynamicHuffman.hpp>
#include <rapidgzip/blockfinder/Uncompressed.hpp>
#include <rapidgzip/ChunkData.hpp>
#include <rapidgzip/gzip/definitions.hpp>
#include <rapidgzip/gzip/deflate.hpp>
#include <rapidgzip/gzip/gzip.hpp>
#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    #include <rapidgzip/gzip/isal.hpp>
#endif
#include <rapidgzip/gzip/zlib.hpp>

#include "DecompressionError.hpp"


namespace rapidgzip
{
template<typename T_ChunkData>
class GzipChunk final
{
public:
    using ChunkData = T_ChunkData;
    using ChunkConfiguration = typename ChunkData::Configuration;
    using Subchunk = typename ChunkData::Subchunk;
    using WindowView = typename ChunkData::WindowView;


    static void
    startNewSubchunk( std::vector<Subchunk>& subchunks,
                      const size_t           encodedOffset )
    {
        const auto nextDecodedOffset = subchunks.empty()
                                       ? 0
                                       : subchunks.back().decodedOffset + subchunks.back().decodedSize;

        auto& subchunk = subchunks.emplace_back();
        subchunk.encodedOffset = encodedOffset;
        subchunk.decodedSize = 0;
        subchunk.decodedOffset = nextDecodedOffset;
    }

    static void
    determineUsedWindowSymbolsForLastSubchunk( std::vector<Subchunk>& subchunks,
                                               gzip::BitReader&       bitReader )
    {
        if ( subchunks.empty() || ( subchunks.back().encodedSize == 0 ) ) {
            return;
        }

        auto& subchunk = subchunks.back();

        /* Only gather sparsity information when it is necessary (non-empty window) or may become necessary
         * (no window yet). The window may already be initialized and empty for deflate bocks after gzip headers. */
        if ( subchunk.window && subchunk.window->empty() ) {
            return;
        }

        try {
            /* Get the window as soon as possible to avoid costly long seeks back outside the BitReader buffer.
             * Especially, don't do it during chunk splitting because it would be too late in general. */
            const Finally seekBack{ [&bitReader, oldOffset = bitReader.tell()] () {
                bitReader.seekTo( oldOffset ); }
            };
            bitReader.seek( subchunk.encodedOffset + subchunk.encodedSize );
            subchunk.usedWindowSymbols = deflate::getUsedWindowSymbols( bitReader );
        } catch ( const std::exception& ) {
            /* Ignore errors such as EOF and even decompression errors because we are only collecting extra
             * data and might already be at the end of the given chunk size, so shouldn't return errors for
             * data thereafter. */
        }

        /* Check whether the no window is needed at all. This may happen when analyzing the very first deflate
         * block and it is at the start of a gzip stream or if the subchunk starts with a non-compressed block. */
        const auto& usedWindowSymbols = subchunk.usedWindowSymbols;
        if ( std::all_of( usedWindowSymbols.begin(), usedWindowSymbols.end(), [] ( bool p ) { return !p; } ) ) {
            subchunk.usedWindowSymbols = std::vector<bool>();  // Free memory!
            subchunk.window = std::make_shared<typename ChunkData::Window>();
        }
    }

    static void
    finalizeWindowForLastSubchunk( ChunkData&             chunk,
                                   std::vector<Subchunk>& subchunks,
                                   gzip::BitReader&       bitReader )
    {
        if ( subchunks.empty() ) {
            return;
        }

        /* Finalize the window of the previous subchunk. Either initialize it to be empty because it is at the
         * start of a new gzip stream and does not need a window, or determine the sparsity. Note that the very
         * first subchunk at offset 0 cannot have a corresponding footer! */
        bool subchunkRequiresWindow{ true };
        const auto nextWindowOffset = subchunks.back().decodedOffset + subchunks.back().decodedSize;
        for ( auto footer = chunk.footers.rbegin(); footer != chunk.footers.rend(); ++footer ) {
            if ( footer->blockBoundary.decodedOffset == nextWindowOffset ) {
                subchunkRequiresWindow = false;
                break;
            }
            /* Footer are sorted ascending and we iterate in reverse order, so we can preemptively quit this
             * search when we find a smaller offset than wanted. This improves performance for many footers
             * as basically only the newly added ones since the last subchunk are checked, resulting in an
             * overall O(n) complexity instead of O(n^2) where n is the number of footers. This is why std::find
             * is not used. */
            if ( footer->blockBoundary.decodedOffset < nextWindowOffset ) {
                break;
            }
        }

        if ( !subchunkRequiresWindow ) {
            subchunks.back().window = std::make_shared<typename ChunkData::Window>();
        } else if ( chunk.configuration.windowSparsity ) {
            determineUsedWindowSymbolsForLastSubchunk( subchunks, bitReader );
        }
    }

    static void
    finalizeChunk( ChunkData&              chunk,
                   std::vector<Subchunk>&& subchunks,
                   gzip::BitReader&        bitReader,
                   const size_t            nextBlockOffset )
    {
        /* Finalize started subchunk. Merge with previous one if it is very small. */
        subchunks.back().encodedSize = nextBlockOffset - subchunks.back().encodedOffset;
        if ( ( subchunks.size() >= 2 )
             && ( subchunks.back().decodedSize < chunk.minimumSplitChunkSize() ) )
        {
            const auto lastSubchunk = subchunks.back();
            subchunks.pop_back();

            subchunks.back().encodedSize += lastSubchunk.encodedSize;
            subchunks.back().decodedSize += lastSubchunk.decodedSize;
            subchunks.back().usedWindowSymbols.clear();
            subchunks.back().window.reset();
        }

        finalizeWindowForLastSubchunk( chunk, subchunks, bitReader );

        chunk.setSubchunks( std::move( subchunks ) );
        chunk.finalize( nextBlockOffset );
    }

    static void
    appendDeflateBlockBoundary( ChunkData&             chunk,
                                std::vector<Subchunk>& subchunks,
                                gzip::BitReader&       bitReader,
                                const size_t           encodedOffset,
                                const size_t           decodedOffset )
    {
        /* Try to append the block boundary, and preemptively return if the boundary already exists. */
        if ( !chunk.appendDeflateBlockBoundary( encodedOffset, decodedOffset ) ) {
            return;
        }

        if ( subchunks.empty() ) {
            return;
        }

        /* Do on-the-fly chunk splitting. */
        if ( subchunks.back().decodedSize >= chunk.configuration.splitChunkSize ) {
            subchunks.back().encodedSize = encodedOffset - subchunks.back().encodedOffset;
            finalizeWindowForLastSubchunk( chunk, subchunks, bitReader );
            startNewSubchunk( subchunks, encodedOffset );
        }
    }

    /**
     * @param decodedSize If given, it is used to avoid overallocations. It is NOT used as a stop condition.
     * @param exactUntilOffset Decompress until this known bit offset in the encoded stream. It must lie on
     *                         a deflate block boundary.
     */
    template<typename InflateWrapper>
    [[nodiscard]] static ChunkData
    decodeChunkWithInflateWrapper( UniqueFileReader          && sharedFileReader,
                                   size_t                 const encodedOffsetInBits,
                                   size_t                 const exactUntilOffset,
                                   WindowView             const initialWindow,
                                   std::optional<size_t>  const decodedSize,
                                   ChunkConfiguration    const& chunkDataConfiguration )
    {
        const auto tStart = now();

        ChunkData result{ chunkDataConfiguration };
        result.encodedOffsetInBits = encodedOffsetInBits;

        gzip::BitReader bitReader( std::move( sharedFileReader ) );
        bitReader.seek( result.encodedOffsetInBits );
        InflateWrapper inflateWrapper( std::move( bitReader ), exactUntilOffset );
        inflateWrapper.setWindow( initialWindow );
        inflateWrapper.setFileType( result.configuration.fileType );

        size_t alreadyDecoded{ 0 };
        while ( true ) {
            const auto suggestedDecodeSize = decodedSize.value_or( ALLOCATION_CHUNK_SIZE );
            deflate::DecodedVector subchunk( suggestedDecodeSize > alreadyDecoded
                                             ? std::min( ALLOCATION_CHUNK_SIZE, suggestedDecodeSize - alreadyDecoded )
                                             : ALLOCATION_CHUNK_SIZE );
            std::optional<Footer> footer;

            /* In order for CRC32 verification to work, we have to append at most one gzip stream per subchunk
             * because the CRC32 calculator is swapped inside ChunkData::append. That's why the stop condition
             * tests for footer.has_value(). */
            size_t nBytesRead = 0;
            size_t nBytesReadPerCall{ 0 };
            for ( ; ( nBytesRead < subchunk.size() ) && !footer; nBytesRead += nBytesReadPerCall ) {
                std::tie( nBytesReadPerCall, footer ) = inflateWrapper.readStream( subchunk.data() + nBytesRead,
                                                                                   subchunk.size() - nBytesRead );
                if ( ( nBytesReadPerCall == 0 ) && !footer ) {
                    break;
                }
            }

            alreadyDecoded += nBytesRead;

            subchunk.resize( nBytesRead );
            result.append( std::move( subchunk ) );
            if ( footer ) {
                footer->blockBoundary.decodedOffset = alreadyDecoded;
                result.appendFooter( *footer );
            }

            if ( ( nBytesReadPerCall == 0 ) && !footer ) {
                break;
            }
        }

        uint8_t dummy{ 0 };
        auto [nBytesReadPerCall, footer] = inflateWrapper.readStream( &dummy, 1 );
        if ( ( nBytesReadPerCall == 0 ) && footer ) {
            footer->blockBoundary.decodedOffset = alreadyDecoded;
            result.appendFooter( *footer );
        }

        if ( exactUntilOffset != inflateWrapper.tellCompressed() ) {
            std::stringstream message;
            message << "The inflate wrapper offset (" << inflateWrapper.tellCompressed() << ") "
                    << "does not match the requested exact stop offset: " << exactUntilOffset << ". "
                    << "The archive or the index may be corrupted or the stop condition might contain a bug. "
                    << "Decoded: " << alreadyDecoded << " B";
            if ( decodedSize ) {
                message << " out of requested " << *decodedSize << " B";
            }
            message << ", started at offset: " << result.encodedOffsetInBits << ".";
            throw std::runtime_error( std::move( message ).str() );
        }

        result.finalize( exactUntilOffset );
        result.statistics.decodeDurationInflateWrapper = duration( tStart );
        return result;
    }


    /**
     * This is called from @ref decodeChunkWithRapidgzip in case the window has been fully resolved so that
     * normal decompression instead of two-staged one becomes possible.
     *
     * @param untilOffset In contrast to @ref decodeChunkWithInflateWrapper, this may be an inexact guess
     *                    from which another thread starts decoding!
     * @note This code is copy-pasted from decodeChunkWithInflateWrapper and adjusted to use the stopping
     *       points and deflate block properties as stop criterion.
     */
    template<typename InflateWrapper>
    [[nodiscard]] static ChunkData
    finishDecodeChunkWithInexactOffset( gzip::BitReader* const  bitReader,
                                        size_t           const  untilOffset,
                                        WindowView       const  initialWindow,
                                        size_t           const  maxDecompressedChunkSize,
                                        ChunkData&&             result,
                                        std::vector<Subchunk>&& subchunks )
    {
        if ( bitReader == nullptr ) {
            throw std::invalid_argument( "BitReader may not be nullptr!" );
        }

        const auto tStart = now();
        auto nextBlockOffset = bitReader->tell();
        bool stoppingPointReached{ false };
        auto alreadyDecoded = result.size();

        if ( ( alreadyDecoded > 0 ) && !bitReader->eof() ) {
            appendDeflateBlockBoundary( result, subchunks, *bitReader, nextBlockOffset, alreadyDecoded );
        }

        InflateWrapper inflateWrapper{ BitReader( *bitReader ) };
        inflateWrapper.setFileType( result.configuration.fileType );
        inflateWrapper.setWindow( initialWindow );
        inflateWrapper.setStoppingPoints( static_cast<StoppingPoint>( StoppingPoint::END_OF_BLOCK |
                                                                      StoppingPoint::END_OF_BLOCK_HEADER |
                                                                      StoppingPoint::END_OF_STREAM_HEADER ) );

        while ( !stoppingPointReached ) {
            deflate::DecodedVector buffer( ALLOCATION_CHUNK_SIZE );
            std::optional<Footer> footer;

            /* In order for CRC32 verification to work, we have to append at most one gzip stream per subchunk
             * because the CRC32 calculator is swapped inside ChunkData::append. */
            size_t nBytesRead = 0;
            size_t nBytesReadPerCall{ 0 };
            while ( ( nBytesRead < buffer.size() ) && !footer && !stoppingPointReached ) {
                std::tie( nBytesReadPerCall, footer ) = inflateWrapper.readStream( buffer.data() + nBytesRead,
                                                                                   buffer.size() - nBytesRead );
                nBytesRead += nBytesReadPerCall;
                subchunks.back().decodedSize += nBytesReadPerCall;

                /* We cannot stop decoding after a final block because the following decoder does not
                 * expect to start a gzip footer. Put another way, we are interested in START_OF_BLOCK
                 * not END_OF_BLOCK and therefore we have to infer one from the other. */
                bool isBlockStart{ false };

                switch ( inflateWrapper.stoppedAt() )
                {
                case StoppingPoint::END_OF_STREAM_HEADER:
                    isBlockStart = true;
                    break;

                case StoppingPoint::END_OF_BLOCK:
                    isBlockStart = !inflateWrapper.isFinalBlock();
                    break;

                case StoppingPoint::END_OF_BLOCK_HEADER:
                    if ( ( ( nextBlockOffset >= untilOffset )
                           && !inflateWrapper.isFinalBlock()
                           && ( inflateWrapper.compressionType() != deflate::CompressionType::FIXED_HUFFMAN ) )
                         || ( nextBlockOffset == untilOffset ) ) {
                        stoppingPointReached = true;
                    }
                    break;

                case StoppingPoint::NONE:
                    if ( ( nBytesReadPerCall == 0 ) && !footer ) {
                        stoppingPointReached = true;
                    }
                    break;

                default:
                    throw std::logic_error( "Got stopping point of a type that was not requested!" );
                }

                if ( isBlockStart ) {
                    nextBlockOffset = inflateWrapper.tellCompressed();

                    /* Do not push back the first boundary because it is redundant as it should contain the same encoded
                     * offset as @ref result and it also would have the same problem that the real offset is ambiguous
                     * for non-compressed blocks. */
                    if ( alreadyDecoded + nBytesRead > 0 ) {
                        appendDeflateBlockBoundary( result, subchunks, *bitReader,
                                                    nextBlockOffset, alreadyDecoded + nBytesRead );
                    }

                    if ( alreadyDecoded >= maxDecompressedChunkSize ) {
                        stoppingPointReached = true;
                        result.stoppedPreemptively = true;
                        break;
                    }
                }
            }

            alreadyDecoded += nBytesRead;

            buffer.resize( nBytesRead );
            result.append( std::move( buffer ) );
            if ( footer ) {
                nextBlockOffset = inflateWrapper.tellCompressed();
                footer->blockBoundary.decodedOffset = alreadyDecoded;
                result.appendFooter( *footer );
            }

            if ( ( inflateWrapper.stoppedAt() == StoppingPoint::NONE )
                 && ( nBytesReadPerCall == 0 ) && !footer ) {
                break;
            }
        }

        uint8_t dummy{ 0 };
        auto [nBytesReadPerCall, footer] = inflateWrapper.readStream( &dummy, 1 );
        if ( ( inflateWrapper.stoppedAt() == StoppingPoint::NONE ) && ( nBytesReadPerCall == 0 ) && footer ) {
            nextBlockOffset = inflateWrapper.tellCompressed();
            footer->blockBoundary.decodedOffset = alreadyDecoded;
            result.appendFooter( *footer );
        }

        finalizeChunk( result, std::move( subchunks ), *bitReader, nextBlockOffset );
        result.statistics.decodeDurationIsal = duration( tStart );
        /**
         * Without the std::move, performance is halved! It seems like copy elision on return does not work
         * with function arguments! @see https://en.cppreference.com/w/cpp/language/copy_elision
         * > In a return statement, when the operand is the name of a non-volatile object with automatic
         * > storage duration, **which isn't a function parameter** [...]
         * And that is only the non-mandatory copy elision, which isn't even guaranteed in C++17!
         */
        return std::move( result );
    }


    [[nodiscard]] static ChunkData
    decodeChunkWithRapidgzip( gzip::BitReader*          const bitReader,
                              size_t                    const untilOffset,
                              std::optional<WindowView> const initialWindow,
                              ChunkConfiguration       const& chunkDataConfiguration )
    {
        if ( bitReader == nullptr ) {
            throw std::invalid_argument( "BitReader must be non-null!" );
        }

        const auto maxDecompressedChunkSize = chunkDataConfiguration.maxDecompressedChunkSize;
        ChunkData result{ chunkDataConfiguration };
        result.encodedOffsetInBits = bitReader->tell();

        /* Initialize metadata for chunk splitting.
         * We will create a new subchunk if the decodedSize exceeds a threshold. */
        std::vector<Subchunk> subchunks;
        startNewSubchunk( subchunks, result.encodedOffsetInBits );

    #ifdef LIBRAPIDARCHIVE_WITH_ISAL
        /** @todo ZlibInflateWrapper is missing the API methods setStoppingPoints, stoppedAt, isFinalBlock, and
         *        and compressionType. If we had that, we could improve performance for -P 1 without ISA-L, by
         *        using zlib-ng! See zlib's inflate's flush parameter:
         *        > The flush parameter of inflate() can be Z_NO_FLUSH, Z_SYNC_FLUSH, Z_FINISH, Z_BLOCK, or Z_TREES.
         * Implementing compressedType would be another beast because we would have to parse the block type ourselves
         * as there seems to be no API. Should be doable though, especially as we need to store the bit offset anyway.
         */
        if ( initialWindow ) {
            return finishDecodeChunkWithInexactOffset<IsalInflateWrapper>(
                bitReader, untilOffset, *initialWindow, maxDecompressedChunkSize,
                std::move( result ), std::move( subchunks ) );
        }
    #endif

        /* If true, then read the gzip header. We cannot simply check the gzipHeader optional because we might
         * start reading in the middle of a gzip stream and will not meet the gzip header for a while or never. */
        bool isAtStreamEnd = false;
        size_t streamBytesRead = 0;
        size_t totalBytesRead = 0;
        bool didReadHeader{ false };

        /* Allocate on heap because it is ~217 kB large!
         * Allocating it once for this whole chunk should be negligible overhead. */
        auto block = std::make_unique<deflate::Block<> >();
        if ( initialWindow ) {
            block->setInitialWindow( *initialWindow );
        }

        /* Loop over possibly gzip streams and deflate blocks. We cannot use GzipReader even though it does
         * something very similar because GzipReader only works with fully decodable streams but we
         * might want to return buffer with placeholders in case we don't know the initial window, yet! */
        size_t nextBlockOffset{ 0 };
    #ifdef LIBRAPIDARCHIVE_WITH_ISAL
        size_t cleanDataCount{ 0 };
    #endif
        while ( true )
        {
            if ( isAtStreamEnd ) {
                const auto headerOffset = bitReader->tell();
                auto error = Error::NONE;

                switch ( result.configuration.fileType )
                {
                case FileType::NONE:
                case FileType::BZIP2:
                    throw std::logic_error( "[GzipChunkFetcher::decodeChunkWithRapidgzip] Invalid file type!" );
                case FileType::BGZF:
                case FileType::GZIP:
                    error = gzip::readHeader( *bitReader ).second;
                    break;
                case FileType::ZLIB:
                    error = zlib::readHeader( *bitReader ).second;
                    break;
                case FileType::DEFLATE:
                    error = Error::NONE;
                    break;
                }

                if ( error != Error::NONE ) {
                    if ( error == Error::END_OF_FILE ) {
                        break;
                    }
                    std::stringstream message;
                    message << "Failed to read gzip/zlib header at offset " << formatBits( headerOffset )
                            << " because of error: " << toString( error );
                    throw std::domain_error( std::move( message ).str() );
                }

            #ifdef LIBRAPIDARCHIVE_WITH_ISAL
                return finishDecodeChunkWithInexactOffset<IsalInflateWrapper>(
                    bitReader, untilOffset, /* initialWindow */ {}, maxDecompressedChunkSize, std::move( result ),
                    std::move( subchunks ) );
            #endif

                didReadHeader = true;
                block->reset( /* initial window */ VectorView<uint8_t>{} );

                isAtStreamEnd = false;
            }

            nextBlockOffset = bitReader->tell();

            if ( totalBytesRead >= maxDecompressedChunkSize ) {
                result.stoppedPreemptively = true;
                break;
            }

        #ifdef LIBRAPIDARCHIVE_WITH_ISAL
            if ( cleanDataCount >= deflate::MAX_WINDOW_SIZE ) {
                return finishDecodeChunkWithInexactOffset<IsalInflateWrapper>(
                    bitReader, untilOffset, result.getLastWindow( {} ), maxDecompressedChunkSize, std::move( result ),
                    std::move( subchunks ) );
            }
        #endif

            if ( auto error = block->readHeader( *bitReader ); error != Error::NONE ) {
                /* Encountering EOF while reading the (first bit for the) deflate block header is only
                 * valid if it is the very first deflate block given to us. Else, it should not happen
                 * because the final block bit should be set in the previous deflate block. */
                if ( ( error == Error::END_OF_FILE ) && ( bitReader->tell() == result.encodedOffsetInBits ) ) {
                    break;
                }

                std::stringstream message;
                message << "Failed to read deflate block header at offset " << formatBits( result.encodedOffsetInBits )
                        << " (position after trying: " << formatBits( bitReader->tell() ) << ": "
                        << toString( error );
                throw std::domain_error( std::move( message ).str() );
            }

            /**
             * Preemptive Stop Condition.
             * @note It is only important for performance that the deflate blocks we are matching here are the same
             *       as the block finder will find.
             * @note We do not have to check for an uncompressed block padding of zero because the deflate decoder
             *       counts that as an error anyway!
             */
            if ( ( ( nextBlockOffset >= untilOffset )
                   && !block->isLastBlock()
                   && ( block->compressionType() != deflate::CompressionType::FIXED_HUFFMAN ) )
                 || ( nextBlockOffset == untilOffset ) ) {
                break;
            }

            /* Do not push back the first boundary because it is redundant as it should contain the same encoded
             * offset as @ref result and it also would have the same problem that the real offset is ambiguous
             * for non-compressed blocks. */
            if ( totalBytesRead > 0 ) {
                appendDeflateBlockBoundary( result, subchunks, *bitReader, nextBlockOffset, totalBytesRead );
            }

            /* Loop until we have read the full contents of the current deflate block-> */
            size_t blockBytesRead{ 0 };
            while ( !block->eob() )
            {
                const auto [bufferViews, error] = block->read( *bitReader, std::numeric_limits<size_t>::max() );
                if ( error != Error::NONE ) {
                    std::stringstream message;
                    message << "Failed to decode deflate block at " << formatBits( result.encodedOffsetInBits )
                            << " because of: " << toString( error );
                    throw std::domain_error( std::move( message ).str() );
                }

            #ifdef LIBRAPIDARCHIVE_WITH_ISAL
                cleanDataCount += bufferViews.dataSize();
            #endif

                result.append( bufferViews );
                blockBytesRead += bufferViews.size();

                /* Non-compressed deflate blocks are limited to 64 KiB and the largest Dynamic Huffman Coding
                 * deflate blocks I have seen were 128 KiB in compressed size. With a maximum compression
                 * ratio of 1032, this would result in ~128 MiB. Fortunately, simple runs of zeros did compress
                 * to only 8 KiB blocks, i.e., ~8 MiB decompressed.
                 * However, igzip -0 can compress the whole file in a single deflate block.
                 * Decompressing such a file is not supported (yet). It would require some heavy
                 * refactoring of the ChunkData class to support resuming the decompression so that
                 * we can simply break and return here instead of throwing an exception. This would basically
                 * require putting a whole GzipReader in the ChunkData so that even random access is supported
                 * in an emulated manner. */
                if ( blockBytesRead > 256_Mi ) {
                    throw std::runtime_error( "A single deflate block that decompresses to more than 256 MiB was "
                                              "encountered. This is not supported to avoid out-of-memory errors." );
                }
            }
            streamBytesRead += blockBytesRead;
            totalBytesRead += blockBytesRead;
            subchunks.back().decodedSize += blockBytesRead;

            if ( block->isLastBlock() ) {
                Footer footer;

                switch ( result.configuration.fileType )
                {
                case FileType::NONE:
                case FileType::BZIP2:
                    throw std::logic_error( "Cannot decode stream if the file type is not specified!" );

                case FileType::DEFLATE:
                    if ( bitReader->tell() % BYTE_SIZE != 0 ) {
                        bitReader->read( BYTE_SIZE - bitReader->tell() % BYTE_SIZE );
                    }
                    break;

                case FileType::ZLIB:
                {
                    footer.zlibFooter = zlib::readFooter( *bitReader );
                    /** @todo check Adler32 checksum when computation has been implemented. */
                    break;
                }

                case FileType::BGZF:
                case FileType::GZIP:
                {
                    footer.gzipFooter = gzip::readFooter( *bitReader );
                    /* We only check for the stream size and CRC32 if we have read the whole stream including the header! */
                    if ( didReadHeader ) {
                        if ( streamBytesRead != footer.gzipFooter.uncompressedSize ) {
                            std::stringstream message;
                            message << "Mismatching size (" << streamBytesRead << " <-> footer: "
                                    << footer.gzipFooter.uncompressedSize << ") for gzip stream!";
                            throw std::runtime_error( std::move( message ).str() );
                        }
                    }
                    break;
                }
                }

                footer.blockBoundary.decodedOffset = totalBytesRead;
                footer.blockBoundary.encodedOffset = bitReader->tell();  // End-of-footer offset for now!
                result.appendFooter( footer );

                isAtStreamEnd = true;
                didReadHeader = false;
                streamBytesRead = 0;

                if ( bitReader->eof() ) {
                    nextBlockOffset = bitReader->tell();
                    break;
                }
            }
        }

        finalizeChunk( result, std::move( subchunks ), *bitReader, nextBlockOffset );
        return result;
    }


    [[nodiscard]] static ChunkData
    decodeChunk( UniqueFileReader                      && sharedFileReader,
                 size_t                             const blockOffset,
                 size_t                             const untilOffset,
                 typename ChunkData::SharedWindow   const initialWindow,
                 std::optional<size_t>              const decodedSize,
                 std::atomic<bool>                 const& cancelThreads,
                 ChunkConfiguration                const& chunkDataConfiguration,
                 bool                               const untilOffsetIsExact = false )
    {
        if ( initialWindow && untilOffsetIsExact ) {
        #ifdef LIBRAPIDARCHIVE_WITH_ISAL
            using InflateWrapper = IsalInflateWrapper;
        #else
            using InflateWrapper = ZlibInflateWrapper;
        #endif

            const auto fileSize = sharedFileReader->size();
            const auto window = initialWindow->decompress();
            auto result = decodeChunkWithInflateWrapper<InflateWrapper>(
                std::move( sharedFileReader ),
                blockOffset,
                fileSize ? std::min( untilOffset, *fileSize * BYTE_SIZE ) : untilOffset,
                *window,
                decodedSize,
                chunkDataConfiguration );

            if ( decodedSize && ( result.decodedSizeInBytes != *decodedSize ) ) {
                std::stringstream message;
                message << "Decoded chunk size does not match the requested decoded size!\n"
                        << "  Block offset          : " << blockOffset << " b\n"
                        << "  Until offset          : " << untilOffset << " b\n"
                        << "  Encoded size          : " << ( untilOffset - blockOffset ) << " b\n"
                        << "  Actual encoded size   : " << result.encodedSizeInBits << " b\n"
                        << "  Decoded size          : " << result.decodedSizeInBytes << " B\n"
                        << "  Expected decoded size : " << *decodedSize << " B\n"
                        << "  Until offset is exact : " << untilOffsetIsExact << "\n"
                        << "  Initial Window        : " << std::to_string( window->size() ) << "\n";
                throw std::runtime_error( std::move( message ).str() );
            }

            return result;
        }

        gzip::BitReader bitReader( std::move( sharedFileReader ) );
        if ( initialWindow ) {
            bitReader.seekTo( blockOffset );
            const auto window = initialWindow->decompress();
            return decodeChunkWithRapidgzip( &bitReader, untilOffset, *window, chunkDataConfiguration );
        }

        const auto tryToDecode =
            [&] ( const std::pair<size_t, size_t>& offset ) -> std::optional<ChunkData>
            {
                try {
                    /* For decoding, it does not matter whether we seek to offset.first or offset.second but it did
                     * matter a lot for interpreting and correcting the encodedSizeInBits in GzipBlockFetcer::get! */
                    bitReader.seekTo( offset.second );
                    auto result = decodeChunkWithRapidgzip(
                        &bitReader, untilOffset, /* initialWindow */ std::nullopt, chunkDataConfiguration );
                    result.encodedOffsetInBits = offset.first;
                    result.maxEncodedOffsetInBits = offset.second;
                    result.encodedSizeInBits = result.encodedEndOffsetInBits - result.encodedOffsetInBits;
                    /** @todo Avoid out of memory issues for very large compression ratios by using a simple runtime
                     *        length encoding or by only undoing the Huffman coding in parallel and the LZ77 serially,
                     *        or by stopping decoding at a threshold and fall back to serial decoding in that case? */
                    return result;
                } catch ( const std::exception& exception ) {
                    /* Ignore errors and try next block candidate. This is very likely to happen if @ref blockOffset
                     * is only an estimated offset! If it happens because decodeChunkWithRapidgzip has a bug, then it
                     * might indirectly trigger an exception when the next required block offset cannot be found. */
                }
                return std::nullopt;
            };

        /* First simply try to decode at the current position to avoid expensive block finders in the case
         * that for some reason the @ref blockOffset guess was perfect. Note that this has to be added as
         * a separate stop condition for decoding the previous block! */
        if ( auto result = tryToDecode( { blockOffset, blockOffset } ); result ) {
            return *std::move( result );
        }

        /**
         * Performance when looking for dynamic and uncompressed blocks:
         * @verbatim
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c 4GiB-base64.gz | wc -c
         *    read in total 5.25998e+09 B out of 3263906195 B, i.e., read the file 1.61156 times
         *    Decompressed in total 4294967296 B in 4.35718 s -> 985.722 MB/s
         *
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c random.bin.gz | wc -c
         *    read in total 2.41669e+09 B out of 2148139037 B, i.e., read the file 1.12502 times
         *    Decompressed in total 2147483648 B in 1.40283 s -> 1530.82 MB/s
         * @endverbatim
         *
         * Performance when looking only for dynamic blocks:
         * @verbatim
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c 4GiB-base64.gz | wc -c
         *    read in total 3.67191e+09 B out of 3263906195 B, i.e., read the file 1.12501 times
         *    Decompressed in total 4294967296 B in 3.0489 s -> 1408.69 MB/s
         *  -> Almost 50% faster! And only 12% file read overhead instead of 61%!
         *
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c random.bin.gz | wc -c
         *   -> LONGER THAN 3 MIN!
         * @endverbatim
         *
         * Performance after implementing the chunked alternating search:
         * @verbatim
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c 4GiB-base64.gz | wc -c
         *    read in total 3.68686e+09 B out of 3263906195 B, i.e., read the file 1.12958 times
         *    Decompressed in total 4294967296 B in 3.06287 s -> 1402.27 MB/s
         *
         * m rapidgzip && time src/tools/rapidgzip -P 0 -d -c random.bin.gz | wc -c
         *    read in total 2.41669e+09 B out of 2148139037 B, i.e., read the file 1.12502 times
         *    Decompressed in total 2147483648 B in 1.31924 s -> 1627.82 MB/s
         * @endverbatim
         */

        const auto findNextDynamic =
            [&] ( size_t beginOffset, size_t endOffset ) {
                if ( beginOffset >= endOffset ) {
                    return std::numeric_limits<size_t>::max();
                }
                bitReader.seekTo( beginOffset );
                return blockfinder::seekToNonFinalDynamicDeflateBlock( bitReader, endOffset );
            };

        const auto findNextUncompressed =
            [&] ( size_t beginOffset, size_t endOffset ) {
                if ( beginOffset >= endOffset ) {
                    return std::make_pair( std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() );
                }
                bitReader.seekTo( beginOffset );
                return blockfinder::seekToNonFinalUncompressedDeflateBlock( bitReader, endOffset );
            };

        /**
         * 1. Repeat for each chunk:
         *    1. Initialize both offsets with possible matches inside the current chunk.
         *    2. Repeat until both offsets are invalid because no further matches were found in the chunk:
         *       1. Try decoding the earlier offset.
         *       2. Update the used offset by searching from last position + 1 until the chunk end.
         */
        const auto tBlockFinderStart = now();
        static constexpr auto CHUNK_SIZE = 8_Ki * BYTE_SIZE;
        size_t falsePositiveCount{ 0 };
        for ( auto chunkBegin = blockOffset; chunkBegin < untilOffset; chunkBegin += CHUNK_SIZE ) {
            /* Only look in the first 512 KiB of data. If nothing can be found there, then something is likely
             * to be wrong with the file and looking in the rest will also likely fail. And looking in the whole
             * range to be decompressed is multiple times slower than decompression because of the slower
             * block finder and because it requires intensive seeking for false non-compressed block positives. */
            if ( cancelThreads || ( chunkBegin - blockOffset >= 512_Ki * BYTE_SIZE ) ) {
                break;
            }

            const auto chunkEnd = std::min( static_cast<size_t>( chunkBegin + CHUNK_SIZE ), untilOffset );

            auto uncompressedOffsetRange = findNextUncompressed( chunkBegin, chunkEnd );
            auto dynamicHuffmanOffset = findNextDynamic( chunkBegin, chunkEnd );

            while ( ( uncompressedOffsetRange.first < chunkEnd ) || ( dynamicHuffmanOffset < chunkEnd ) ) {
                if ( cancelThreads ) {
                    break;
                }

                /* Choose the lower offset to test next. */
                std::pair<size_t, size_t> offsetToTest;
                if ( dynamicHuffmanOffset < uncompressedOffsetRange.first ) {
                    offsetToTest = { dynamicHuffmanOffset, dynamicHuffmanOffset };
                    dynamicHuffmanOffset = findNextDynamic( dynamicHuffmanOffset + 1, chunkEnd );
                } else {
                    offsetToTest = uncompressedOffsetRange;
                    uncompressedOffsetRange = findNextUncompressed( uncompressedOffsetRange.second + 1, chunkEnd );
                }

                /* Try decoding and measure the time. */
                const auto tBlockFinderStop = now();
                if ( auto result = tryToDecode( offsetToTest ); result ) {
                    result->statistics.blockFinderDuration = duration( tBlockFinderStart, tBlockFinderStop );
                    result->statistics.decodeDuration = duration( tBlockFinderStop );
                    result->statistics.falsePositiveCount = falsePositiveCount;
                    return *std::move( result );
                }

                falsePositiveCount++;
            }
        }

        std::stringstream message;
        message << "Failed to find any valid deflate block in [" << formatBits( blockOffset )
                << ", " << formatBits( untilOffset ) << ")";
        throw NoBlockInRange( std::move( message ).str() );
    }
};
}  // namespace rapidgzip
