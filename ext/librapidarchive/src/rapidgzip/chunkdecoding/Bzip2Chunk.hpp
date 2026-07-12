#pragma once

#include <algorithm>
#include <atomic>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <core/BitStringFinder.hpp>
#include <indexed_bzip2/bzip2.hpp>
#include <rapidgzip/ChunkData.hpp>
#include <rapidgzip/gzip/definitions.hpp>

#include "DecompressionError.hpp"


namespace rapidgzip
{
using Configuration = ChunkData::Configuration;


template<typename T_ChunkData>
class Bzip2Chunk
{
public:
    using ChunkData = T_ChunkData;
    using ChunkConfiguration = typename ChunkData::Configuration;
    using Window = typename ChunkData::Window;
    using Subchunk = typename ChunkData::Subchunk;


    [[nodiscard]] static ChunkData
    decodeUnknownBzip2Chunk( bzip2::BitReader*      const bitReader,
                             size_t                 const untilOffset,
                             std::optional<size_t>  const decodedSize,
                             ChunkConfiguration    const& chunkDataConfiguration )
    {
        if ( bitReader == nullptr ) {
            throw std::invalid_argument( "BitReader must be non-null!" );
        }

        ChunkData result{ chunkDataConfiguration };
        result.encodedOffsetInBits = bitReader->tell();
        result.maxEncodedOffsetInBits = result.encodedOffsetInBits;

        /* Initialize metadata for chunk splitting.
         * We will create a new subchunk if the decodedSize exceeds a threshold. */
        std::vector<Subchunk> subchunks;
        subchunks.emplace_back();
        subchunks.back().encodedOffset = result.encodedOffsetInBits;
        subchunks.back().decodedOffset = 0;
        subchunks.back().decodedSize = 0;

        /* If true, then read the gzip header. We cannot simply check the gzipHeader optional because we might
         * start reading in the middle of a gzip stream and will not meet the gzip header for a while or never. */
        bool isAtStreamEnd = false;
        size_t totalBytesRead = 0;

        /* Loop over possibly gzip streams and deflate blocks. We cannot use GzipReader even though it does
         * something very similar because GzipReader only works with fully decodable streams but we
         * might want to return buffer with placeholders in case we don't know the initial window, yet! */
        size_t nextBlockOffset{ 0 };
        while ( true )
        {
            if ( isAtStreamEnd ) {
                bzip2::readBzip2Header( *bitReader );
                isAtStreamEnd = false;
            }

            nextBlockOffset = bitReader->tell();

            /* Do on-the-fly chunk splitting. */
            if ( subchunks.back().decodedSize >= result.configuration.splitChunkSize ) {
                subchunks.back().encodedSize = nextBlockOffset - subchunks.back().encodedOffset;

                const auto nextdecodedOffset = subchunks.back().decodedOffset + subchunks.back().decodedSize;

                auto& subchunk = subchunks.emplace_back();
                subchunk.encodedOffset = nextBlockOffset;
                subchunk.decodedOffset = nextdecodedOffset;
                subchunk.decodedSize = 0;
            }

            /** @todo does this work when quitting on an empty block, i.e., if the next block is an end-of-stream one?
             *        test decodeUnknownBzip2Chunk with all block offsets */
            if ( totalBytesRead >= chunkDataConfiguration.maxDecompressedChunkSize ) {
                result.stoppedPreemptively = true;
                break;
            }

            /* This also reads the header and will throw on errors. */
            std::optional<bzip2::Block> block;
            try {
                block.emplace( *bitReader );
                if ( !block->eos() ) {
                    block->readBlockData();
                }
            } catch ( const bzip2::BitReader::EndOfFileReached& exception ) {
                /* Encountering EOF while reading the (first bit for the) deflate block header is only
                 * valid if it is the very first deflate block given to us. Else, it should not happen
                 * because the final block bit should be set in the previous deflate block. */
                if ( bitReader->tell() == result.encodedOffsetInBits ) {
                    break;
                }
                throw exception;
            }

            /* Preemptive Stop Condition. */
            if ( ( ( nextBlockOffset >= untilOffset ) && !block->eos() ) || ( nextBlockOffset == untilOffset ) ) {
                break;
            }

            /* Do not push back the first boundary because it is redundant as it should contain the same encoded
             * offset as @ref result and it also would have the same problem that the real offset is ambiguous
             * for non-compressed blocks. */
            if ( totalBytesRead > 0 ) {
                result.appendDeflateBlockBoundary( nextBlockOffset, totalBytesRead );
            }

            /* In contrast to deflate, bzip2 has dedicated end-of-stream blocks, which do not contain any data.
             * Therefore, we need to check for it before calling block.decodeBlock. */
            if ( block->eos() ) {
                const auto footerOffset = bitReader->tell();
                Footer footer;
                footer.blockBoundary = { footerOffset, totalBytesRead };
                result.appendFooter( std::move( footer ) );

                isAtStreamEnd = true;

                if ( bitReader->eof() ) {
                    nextBlockOffset = bitReader->tell();
                    break;
                }
                continue;
            }

            /* Loop until we have read the full contents of the current block. */
            size_t blockBytesRead{ 0 };
            while ( true )
            {
                const auto suggestedDecodeSize = decodedSize.value_or( ALLOCATION_CHUNK_SIZE );
                deflate::DecodedVector subchunk(
                    suggestedDecodeSize > totalBytesRead
                    ? std::min( ALLOCATION_CHUNK_SIZE, suggestedDecodeSize - totalBytesRead )
                    : ALLOCATION_CHUNK_SIZE );

                size_t nBytesRead = 0;
                while ( nBytesRead < subchunk.size() ) {
                    const auto nBytesReadPerCall = block->read(
                        subchunk.size() - nBytesRead,
                        reinterpret_cast<char*>( subchunk.data() ) + nBytesRead );
                    if ( nBytesReadPerCall == 0 ) {
                        break;
                    }
                    nBytesRead += nBytesReadPerCall;
                }

                subchunk.resize( nBytesRead );
                result.append( std::move( subchunk ) );

                if ( nBytesRead == 0 ) {
                    break;
                }

                blockBytesRead += nBytesRead;
                totalBytesRead += nBytesRead;
                subchunks.back().decodedSize += nBytesRead;

                /**
                 * > Because of the first-stage RLE compression (see above), the maximum length of plaintext
                 * > that a single 900 kB bzip2 block can contain is around 46 MB (45,899,236 bytes).
                 * @see https://en.wikipedia.org/wiki/Bzip2
                 * Note that @ref maxDecompressedChunkSize is still necessary because this only limits the
                 * decoded size of a single bzip2 block, while a chunk can contain multiple such blocks.
                 * > An even smaller file of 40 bytes can be achieved by using an input containing entirely
                 * > values of 251, an apparent compression ratio of 1147480.9:1.
                 * This makes chunk splitting and @ref maxDecompressedChunkSize still a requirement.
                 */
                if ( blockBytesRead > 64_Mi ) {
                    throw std::runtime_error( "A single bzip2 block that decompresses to more than 64 MiB was "
                                              "encountered. This is not supported to avoid out-of-memory errors." );
                }
            }
        }

        /* Finalize started subchunk. Merge with previous one if it is very small. */
        subchunks.back().encodedSize = nextBlockOffset - subchunks.back().encodedOffset;
        if ( ( subchunks.size() >= 2 )
             && ( subchunks.back().decodedSize < result.minimumSplitChunkSize() ) )
        {
            const auto lastSubchunk = subchunks.back();
            subchunks.pop_back();

            subchunks.back().encodedSize += lastSubchunk.encodedSize;
            subchunks.back().decodedSize += lastSubchunk.decodedSize;
        }

        /* Ensure that all subchunks have empty windows to avoid them being filled as windows are not necessary. */
        for ( auto& subchunk : subchunks ) {
            subchunk.window = std::make_shared<Window>();
        }

        result.setSubchunks( std::move( subchunks ) );
        result.finalize( nextBlockOffset );
        return result;
    }


    [[nodiscard]] static ChunkData
    decodeChunk( UniqueFileReader       && sharedFileReader,
                 size_t              const chunkOffset,
                 size_t              const untilOffset,
                 std::atomic<bool>  const& cancelThreads,
                 ChunkConfiguration const& chunkDataConfiguration )
    {
        bzip2::BitReader bitReader( sharedFileReader->clone() );

        const auto tryToDecode =
            [&] ( const size_t offset ) -> std::optional<ChunkData>
            {
                try {
                    bitReader.seekTo( offset );
                    auto result = decodeUnknownBzip2Chunk( &bitReader, untilOffset, /* decodedSize */ std::nullopt,
                                                           chunkDataConfiguration );
                    return result;
                } catch ( const std::exception& ) {
                    /* Ignore errors and try next block candidate. This is very likely to happen if @ref blockOffset
                     * is only an estimated offset! */
                }
                return std::nullopt;
            };

        if ( auto result = tryToDecode( chunkOffset ); result ) {
            return *std::move( result );
        }

        const auto blockFinderOffsetInBytes = chunkOffset / BYTE_SIZE;
        sharedFileReader->seekTo( blockFinderOffsetInBytes );
        BitStringFinder<bzip2::MAGIC_BITS_SIZE> blockFinder{
            std::move( sharedFileReader ), bzip2::MAGIC_BITS_BLOCK, /* fileBufferSizeBytes */ 64_Ki };
        while ( !cancelThreads ) {
            const auto foundRelativeOffset = blockFinder.find();
            if ( foundRelativeOffset == std::numeric_limits<size_t>::max() ) {
                break;
            }

            const auto blockOffset = blockFinderOffsetInBytes * BYTE_SIZE + foundRelativeOffset;
            if ( blockOffset >= untilOffset ) {
                break;
            }

            if ( blockOffset >= chunkOffset ) {
                auto result = tryToDecode( blockOffset );
                if ( result ) {
                    return *std::move( result );
                }
            }
        }

        std::stringstream message;
        message << "Failed to find any valid bzip2 block in [" << formatBits( chunkOffset )
                << ", " << formatBits( untilOffset ) << ")";
        throw NoBlockInRange( std::move( message ).str() );
    }
};
}  // namespace rapidgzip
