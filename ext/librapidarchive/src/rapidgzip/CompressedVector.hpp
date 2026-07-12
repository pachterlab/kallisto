#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <core/FasterVector.hpp>
#include <core/VectorView.hpp>

#include "gzip/InflateWrapper.hpp"
#ifdef LIBRAPIDARCHIVE_WITH_ISAL
    #include "gzip/isal.hpp"
#endif
#include "gzip/zlib.hpp"


namespace rapidgzip
{
enum class CompressionType : uint8_t
{
    NONE      =  0,
    DEFLATE   =  1,
    ZLIB      =  2,
    GZIP      =  3,
    BZIP2     =  4,
    LZ4       =  5,
    ZSTANDARD =  6,
    LZMA      =  7,
    XZ        =  8,
    BROTLI    =  9,
    LZIP      = 10,
    LZOP      = 11,
};


[[nodiscard]] constexpr const char*
toString( const CompressionType compressionType )
{
    switch ( compressionType )
    {
    case CompressionType::NONE     : return "NONE";
    case CompressionType::DEFLATE  : return "Deflate";
    case CompressionType::ZLIB     : return "ZLIB";
    case CompressionType::GZIP     : return "GZIP";
    case CompressionType::BZIP2    : return "BZIP2";
    case CompressionType::LZ4      : return "LZ4";
    case CompressionType::ZSTANDARD: return "ZStandard";
    case CompressionType::LZMA     : return "LZMA";
    case CompressionType::XZ       : return "XZ";
    case CompressionType::BROTLI   : return "Brotli";
    case CompressionType::LZIP     : return "LZIP";
    case CompressionType::LZOP     : return "LZOP";
    }
    return "Unknown";
}


/**
 * @p compressionType May also be NONE but in order to avoid unnecessary copies, it should be avoided.
 */
template<typename Container = FasterVector<uint8_t> >
[[nodiscard]] Container
compress( const VectorView<typename Container::value_type> toCompress,
          const CompressionType                            compressionType )
{
    using namespace rapidgzip;

    switch ( compressionType )
    {
    case CompressionType::GZIP:
    #ifdef LIBRAPIDARCHIVE_WITH_ISAL
        try {
            return rapidgzip::compressWithIsal<Container>( toCompress );
        } catch ( const std::runtime_error& exception ) {
            std::stringstream message;
            message << "[Warning] Compression with ISA-L failed unexpectedly with: " << exception.what() << "\n";
            message << "[Warning] Will use zlib as a fallback. Please report this bug anyway.\n";
        #ifdef RAPIDGZIP_FATAL_PERFORMANCE_WARNINGS
            throw std::logic_error( std::move( message ).str() );
        #else
            std::cerr << std::move( message ).str();
        #endif
            return rapidgzip::compressWithZlib<Container>( toCompress );
        }
    #else
        return rapidgzip::compressWithZlib<Container>( toCompress );
    #endif
        break;

    case CompressionType::ZLIB:
        return rapidgzip::compressWithZlib<Container>( toCompress, CompressionStrategy::DEFAULT, /* dictionary */ {},
                                                       ContainerFormat::ZLIB );

    case CompressionType::NONE:
        return Container( toCompress.begin(), toCompress.end() );
        break;

    default:
        break;
    }

    throw std::invalid_argument( std::string( "Only gzip compression and none are currently supported" )
                                 + ", but got: " + toString( compressionType ) );
}


/**
 * The methods by design are not called simply "data"/"size" to avoid it being used the wrong way
 * when replacing normal containers for this one.
 */
template<typename Container = FasterVector<uint8_t> >
class CompressedVector
{
public:
    using container_type = Container;

public:
    CompressedVector() = default;

    explicit
    CompressedVector( Container&&           toCompress,
                      const CompressionType compressionType ) :
        m_compressionType( compressionType ),
        m_decompressedSize( toCompress.size() ),
        m_data( compressionType == CompressionType::NONE
                ? std::make_shared<Container>( std::move( toCompress ) )
                : std::make_shared<Container>( compress<Container>( toCompress, compressionType ) ) )
    {}

    explicit
    CompressedVector( const VectorView<typename Container::value_type> toCompress,
                      const CompressionType                            compressionType ) :
        m_compressionType( compressionType ),
        m_decompressedSize( toCompress.size() ),
        m_data( std::make_shared<Container>( compress<Container>( toCompress, compressionType ) ) )
    {}

    explicit
    CompressedVector( Container&&           compressedData,
                      const size_t          decompressedSize,
                      const CompressionType compressionType ) :
        m_compressionType( compressionType ),
        m_decompressedSize( decompressedSize ),
        m_data( std::make_shared<Container>( std::move( compressedData ) ) )
    {}

    [[nodiscard]] CompressionType
    compressionType() const noexcept
    {
        return m_compressionType;
    }

    [[nodiscard]] std::shared_ptr<const Container>
    compressedData() const noexcept
    {
        return m_data ? m_data : std::make_shared<Container>();
    }

    [[nodiscard]] size_t
    compressedSize() const noexcept
    {
        return m_data ? m_data->size() : 0;
    }

    /**
     * @return a non-null shared pointer to the decompressed data!
     */
    [[nodiscard]] std::shared_ptr<const Container>
    decompress() const
    {
        if ( !m_data || empty() ) {
            return std::make_shared<Container>();
        }

        #ifdef LIBRAPIDARCHIVE_WITH_ISAL
            using InflateWrapper = rapidgzip::IsalInflateWrapper;
        #else
            using InflateWrapper = rapidgzip::ZlibInflateWrapper;
        #endif

        const auto decompressWithWrapper =
            [this] ( FileType fileType ) {
                return std::make_shared<Container>(
                    inflateWithWrapper<InflateWrapper>( *m_data, m_decompressedSize, /* window */ {}, fileType ) );
            };

        switch ( m_compressionType )
        {
        case CompressionType::DEFLATE:
            return decompressWithWrapper( FileType::DEFLATE );
        case CompressionType::GZIP:
            return decompressWithWrapper( FileType::GZIP );
        case CompressionType::ZLIB:
            return decompressWithWrapper( FileType::ZLIB );
        case CompressionType::NONE:
            return m_data;
        default:
            break;
        }

        throw std::invalid_argument( std::string( "Only gzip compression and none are currently supported" )
                                     + ", but got: " + toString( m_compressionType ) );
    }

    [[nodiscard]] size_t
    decompressedSize() const noexcept
    {
        return m_decompressedSize;
    }

    void
    clear()
    {
        m_data.reset();
        m_decompressedSize = 0;
    }

    [[nodiscard]] bool
    empty() const noexcept
    {
        return m_decompressedSize == 0;
    }

    [[nodiscard]] bool
    operator==( const CompressedVector& other ) const
    {
        return ( m_compressionType == other.m_compressionType )
               && ( static_cast<bool>( m_data ) == static_cast<bool>( other.m_data ) )
               && ( !m_data || !other.m_data || ( *m_data == *other.m_data ) )
               && ( m_decompressedSize == other.m_decompressedSize );
    }

    [[nodiscard]] bool
    operator!=( const CompressedVector& other ) const
    {
        return !( *this == other );
    }

private:
    CompressionType m_compressionType{ CompressionType::GZIP };
    size_t m_decompressedSize{ 0 };
    std::shared_ptr<const Container> m_data;
};
}  // namespace rapidgzip
