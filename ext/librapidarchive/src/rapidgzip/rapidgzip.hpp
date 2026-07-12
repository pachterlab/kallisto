#pragma once

#include <array>
#include <cstdint>
#include <string>

#include "gzip/deflate.hpp"
#include "gzip/gzip.hpp"
#include "gzip/GzipReader.hpp"
#include "IndexFileFormat.hpp"
#include "ParallelGzipReader.hpp"


static constexpr uint32_t RAPIDGZIP_VERSION_MAJOR{ 0 };
static constexpr uint32_t RAPIDGZIP_VERSION_MINOR{ 16 };
static constexpr uint32_t RAPIDGZIP_VERSION_PATCH{ 0 };
static constexpr uint32_t RAPIDGZIP_VERSION{
    RAPIDGZIP_VERSION_MAJOR * 0x10000UL + RAPIDGZIP_VERSION_MINOR * 0x100UL + RAPIDGZIP_VERSION_PATCH
};


namespace rapidgzip
{
static constexpr std::array<uint8_t, 3> VERSION = {
    RAPIDGZIP_VERSION_MAJOR,
    RAPIDGZIP_VERSION_MINOR,
    RAPIDGZIP_VERSION_PATCH,
};


static const std::string VERSION_STRING{
    std::to_string( RAPIDGZIP_VERSION_MAJOR ) + '.' +
    std::to_string( RAPIDGZIP_VERSION_MINOR ) + '.' +
    std::to_string( RAPIDGZIP_VERSION_PATCH )
};
}  // namespace rapidgzip
