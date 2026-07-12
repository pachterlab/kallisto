#pragma once

#include <algorithm>
#include <cassert>
#include <deque>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>
#include <vector>

#include <core/BlockFinderInterface.hpp>
#include <core/common.hpp>
#include <indexed_bzip2/bzip2.hpp>

#include "blockfinder/Bgzf.hpp"
#include "gzip/deflate.hpp"
#include "gzip/format.hpp"


namespace rapidgzip
{
/**
 * This is a much more lean version of core/BlockFinder. It does not do any actual work aside from finding the first
 * deflate block. Instead, it is mostly doing bookkeeping and simple partitioning using @ref m_spacingInBits to generate
 * guesses beyond the known block offsets and inside the file range.
 *
 * Block offsets can be confirmed, in which case those will be returned. This is important for performant
 * prefetching and is hard to let the BlockMap do.
 * However, care has to be taken in its usage because block confirmation effectively invalidates block
 * previous indexes!
 */
class GzipBlockFinder final :
    public BlockFinderInterface
{
public:
    using BlockOffsets = std::vector<size_t>;

public:
    explicit
    GzipBlockFinder( UniqueFileReader fileReader,
                     size_t           spacing ) :
        m_file( std::move( fileReader ) ),
        m_fileSizeInBits( m_file->size()
                          ? std::make_optional( *m_file->size() * CHAR_BIT )
                          : std::nullopt ),
        m_spacingInBits( spacing * CHAR_BIT )
    {
        if ( m_spacingInBits < 32_Ki ) {
            /* Well, actually, it could make sense because this is about the spacing in the compressed data but
             * then even more! A spacing of 32 KiB in uncompressed data can lead to index sizes up to the
             * decompressed file. A spacing of 32 KiB in the compressed data can only lead to an index equal that
             * of the compressed file, so it behaves much more reasonable! */
            throw std::invalid_argument( "A spacing smaller than the window size makes no sense!" );
        }

        const auto detectedFormat = determineFileTypeAndOffset( m_file );
        if ( !detectedFormat ) {
            throw std::invalid_argument( "Failed to detect a valid file format." );
        }

        m_fileType = detectedFormat->first;
        if ( m_fileType == FileType::BGZF ) {
            m_bgzfBlockFinder = std::make_unique<blockfinder::Bgzf>( m_file->clone() );
        }

        m_blockOffsets.push_back( detectedFormat->second );
    }

    /**
     * @return number of block offsets. This number may increase as long as it is not finalized yet.
     */
    [[nodiscard]] size_t
    size() const override
    {
        const std::scoped_lock lock( m_mutex );
        return m_blockOffsets.size();
    }

    void
    finalize()
    {
        const std::scoped_lock lock( m_mutex );
        m_finalized = true;
    }

    [[nodiscard]] bool
    finalized() const override
    {
        const std::scoped_lock lock( m_mutex );
        return m_finalized;
    }

    [[nodiscard]] FileType
    fileType() const noexcept
    {
        return m_fileType;
    }

    /**
     * Insert a known to be exact block offset. They should in general be inserted in sequence because no
     * partitioning will be done before the largest inserted block offset.
     */
    void
    insert( size_t blockOffset )
    {
        const std::scoped_lock lock( m_mutex );
        insertUnsafe( blockOffset );
    }

    using BlockFinderInterface::get;

    /**
     * @return The block offset to the given block index or nothing when the block finder is finalized and the
     *         requested block out of range. When the requested block index is not a known one, a guess will
     *         be returned based on @ref m_spacingInBits.
     * @todo ADD TESTS FOR THIS
     */
    [[nodiscard]] std::pair<std::optional<size_t>, GetReturnCode>
    get( size_t                  blockIndex,
         [[maybe_unused]] double timeoutInSeconds ) override
    {
        const std::scoped_lock lock( m_mutex );

        if ( m_fileType == FileType::BGZF ) {
            return getBgzfBlock( blockIndex );
        }

        if ( blockIndex < m_blockOffsets.size() ) {
            return { m_blockOffsets[blockIndex], GetReturnCode::SUCCESS };
        }

        assert( !m_blockOffsets.empty() );
        const auto blockIndexOutside = blockIndex - m_blockOffsets.size();  // >= 0
        const auto partitionIndex = firstPartitionIndex() + blockIndexOutside;
        const auto blockOffset = partitionIndex * m_spacingInBits;

        if ( !m_fileSizeInBits ) {
            if ( const auto fileSize = m_file->size() ) {
                m_fileSizeInBits = *fileSize * CHAR_BIT;
            }
        }
        if ( !m_fileSizeInBits || ( blockOffset < *m_fileSizeInBits ) ) {
            return { blockOffset, GetReturnCode::SUCCESS };
        }

        /* Return the file size as offset for all indexes past the file.
         * This avoids:
         *  - the BlockFetcher waiting until this index becomes "available"
         *  - the previous index offset not being used because there is no untilOffset for it */
        if ( partitionIndex > 0 ) {
            return { *m_fileSizeInBits, GetReturnCode::FAILURE };
        }

        /* This shouldn't happen. */
        return { 0, GetReturnCode::FAILURE };
    }

    /**
     * @return Index for the block at the requested offset.
     */
    [[nodiscard]] size_t
    find( size_t encodedBlockOffsetInBits ) const override
    {
        const std::scoped_lock lock( m_mutex );

        /* Find in sorted vector by bisection. */
        const auto match = std::lower_bound( m_blockOffsets.begin(), m_blockOffsets.end(), encodedBlockOffsetInBits );
        if ( ( match != m_blockOffsets.end() ) && ( *match == encodedBlockOffsetInBits ) ) {
            return std::distance( m_blockOffsets.begin(), match );
        }

        if ( ( encodedBlockOffsetInBits > m_blockOffsets.back() )
             && ( encodedBlockOffsetInBits % m_spacingInBits == 0 ) )
        {
            const auto blockIndex = m_blockOffsets.size()
                                    + ( encodedBlockOffsetInBits / m_spacingInBits - firstPartitionIndex() );
            assert( ( firstPartitionIndex() + ( blockIndex - m_blockOffsets.size() ) ) * m_spacingInBits
                    == encodedBlockOffsetInBits /* see get for the inverse calculation this is taken from. */ );
            return blockIndex;
        }

        throw std::out_of_range( "No block with the specified offset " + std::to_string( encodedBlockOffsetInBits )
                                 + " exists in the block finder map!" );
    }

    void
    setBlockOffsets( const std::vector<size_t>& blockOffsets )
    {
        m_blockOffsets.assign( blockOffsets.begin(), blockOffsets.end() );
        finalize();
    }

    [[nodiscard]] size_t
    partitionOffsetContainingOffset( size_t blockOffset ) const
    {
        /* Round down to m_spacingInBits grid. */
        return ( blockOffset / m_spacingInBits ) * m_spacingInBits;
    }

    [[nodiscard]] constexpr size_t
    spacingInBits() const noexcept
    {
        return m_spacingInBits;
    }

private:
    [[nodiscard]] std::optional<size_t>
    fileSize()
    {
        if ( m_fileSizeInBits ) {
            return m_fileSizeInBits;
        }

        const auto fileSize = m_file->size();
        if ( fileSize ) {
            m_fileSizeInBits = *fileSize * CHAR_BIT;
            return m_fileSizeInBits;
        }

        return std::nullopt;
    }

    bool
    insertUnsafe( size_t blockOffset )
    {
        const auto size = fileSize();
        if ( size.has_value() && ( blockOffset >= *size ) ) {
            return false;
        }

        const auto match = std::lower_bound( m_blockOffsets.begin(), m_blockOffsets.end(), blockOffset );
        if ( ( match == m_blockOffsets.end() ) || ( *match != blockOffset ) ) {
            if ( m_finalized ) {
                throw std::invalid_argument( "Already finalized, may not insert further block offsets!" );
            }
            m_blockOffsets.insert( match, blockOffset );
            assert( std::is_sorted( m_blockOffsets.begin(), m_blockOffsets.end() ) );
        }

        return true;
    }

    void
    gatherMoreBgzfBlocks( size_t blockIndex )
    {
        while ( blockIndex + m_batchFetchCount >= m_blockOffsets.size() ) {
            const auto nextOffset = m_bgzfBlockFinder->find();
            if ( nextOffset < m_blockOffsets.back() + m_spacingInBits ) {
                continue;
            }
            if ( !insertUnsafe( nextOffset ) ) {
                break;
            }
        }
    }

    [[nodiscard]] std::pair<std::optional<size_t>, GetReturnCode>
    getBgzfBlock( size_t blockIndex )
    {
        if ( m_bgzfBlockFinder && !m_finalized ) {
            gatherMoreBgzfBlocks( blockIndex );
        }

        if ( blockIndex < m_blockOffsets.size() ) {
            return { m_blockOffsets[blockIndex], GetReturnCode::SUCCESS };
        }

        /* Size should be available at this point because EOF should be the only cause
         * for gatherMoreBgzfBlocks not gathering up to the specified index. */
        return { fileSize().value_or( std::numeric_limits<size_t>::max() ), GetReturnCode::FAILURE };
    }

    /**
     * @return the "index" corresponding to the first "guessed" block offset given by the formula i * m_spacingInBits
     *         for i in N_0 with the requirement that it must be larger (not equal) than the last confirmed offset.
     */
    [[nodiscard]] size_t
    firstPartitionIndex() const
    {
        /* Consider a spacing of 2. The guesses would return offsets at 0, 2, 4, 6, ...
         * If the last confirmed offset was 0 or 1 , then the next partition offset would be 2, i.e.,
         * we should return the index 1. If the last confirmed offset was 2 or 3, we should return 2 and so on.
         * This means we want to divide by the spacing and round the result down and add plus 1 to that. */
        return m_blockOffsets.back() / m_spacingInBits + 1;
    }

private:
    mutable std::mutex m_mutex;

    const UniqueFileReader m_file;
    std::optional<size_t> m_fileSizeInBits;
    bool m_finalized{ false };
    const size_t m_spacingInBits;

    /**
     * These should only contain confirmed block offsets in order. Use a deque to avoid having to move all
     * subsequent elements when inserting into the sorted container.
     */
    std::deque<size_t> m_blockOffsets;

    /** Only used for Bgzf files in which case it will gather offsets in chunks of this. */
    FileType m_fileType{ FileType::NONE };
    std::unique_ptr<blockfinder::Bgzf> m_bgzfBlockFinder;
    const size_t m_batchFetchCount = std::max<size_t>( 16, 3ULL * std::thread::hardware_concurrency() );
};
}  // namespace rapidgzip
