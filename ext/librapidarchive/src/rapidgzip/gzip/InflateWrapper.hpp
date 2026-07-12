#pragma once

#include <memory>
#include <optional>
#include <utility>

#include <core/VectorView.hpp>
#include <filereader/BufferView.hpp>

#include "definitions.hpp"
#include "gzip.hpp"
#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    #include "isal.hpp"
#endif

namespace rapidgzip
{
template<typename InflateWrapper,
         typename Container>
[[nodiscard]] Container
inflateWithWrapper( const Container&            toDecompress,
                    const std::optional<size_t> decompressedSize,
                    const VectorView<uint8_t>   dictionary = {},
                    const rapidgzip::FileType   fileType = rapidgzip::FileType::GZIP )
{
    if ( ( decompressedSize && ( *decompressedSize == 0 ) ) || toDecompress.empty() ) {
        return {};
    }

#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    if constexpr ( std::is_same_v<InflateWrapper, rapidgzip::IsalInflateWrapper> ) {
        if ( decompressedSize && dictionary.empty() ) {
            return rapidgzip::inflateWithIsal( toDecompress, *decompressedSize, fileType );
        }
    }
#endif

    gzip::BitReader bitReader(
        std::make_unique<BufferViewFileReader>( toDecompress.data(), toDecompress.size() ) );

    InflateWrapper inflateWrapper( std::move( bitReader ) );

    switch ( fileType )
    {
    case FileType::DEFLATE:
        inflateWrapper.setFileType( fileType );
        break;
    case FileType::BGZF:
    case FileType::GZIP:
    case FileType::ZLIB:
        inflateWrapper.setFileType( fileType );
        inflateWrapper.setStartWithHeader( true );
        break;
    default:
        throw std::invalid_argument( std::string( "Unsupported file type: " ) + toString( fileType ) );
    }

    if ( !dictionary.empty() ) {
        inflateWrapper.setWindow( dictionary );
    }

    static constexpr auto CHUNK_SIZE = 4_Ki;
    Container result;
    while ( true ) {
        const auto oldSize = result.size();
        if ( ( oldSize == 0 ) && decompressedSize && ( *decompressedSize > 0 ) ) {
            result.resize( *decompressedSize );
        } else {
            result.resize( oldSize + CHUNK_SIZE );
        }

        const auto [nBytesRead, footer] = inflateWrapper.readStream( result.data() + oldSize, result.size() - oldSize );
        result.resize( oldSize + nBytesRead );
        if ( ( nBytesRead == 0 ) && !footer ) {
            break;
        }
    }
    return result;
}
}  // namespace rapidgzip
