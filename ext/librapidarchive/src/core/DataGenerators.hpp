#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include "common.hpp"

namespace rapidgzip
{
inline void
createRandomTextFile( const std::string& path,
                      const uint64_t     size )
{
    std::ofstream textFile( path, std::ios_base::out | std::ios_base::binary );
    for ( uint64_t i = 0; i < size; ++i ) {
        const auto c = i % 80 == 0 ? '\n' : 'A' + ( rand() % ( 'Z' - 'A' ) );
        textFile << static_cast<char>( c );
    }
}


inline void
fillWithRandomData( void* const  data,
                    const size_t size )
{
    std::mt19937_64 randomEngine;
    std::array<uint64_t, 1_Ki> buffer{};  // 8 KiB of buffer
    for ( size_t nBytesWritten = 0; nBytesWritten < size; ) {
        for ( auto& x : buffer ) {
            x = randomEngine();
        }
        const auto nBytesToWrite = std::min<uint64_t>( buffer.size() * sizeof( buffer[0] ), size - nBytesWritten );
        std::memcpy( reinterpret_cast<char*>( data ) + nBytesWritten, buffer.data(), nBytesToWrite);
        nBytesWritten += nBytesToWrite;
    }
}


inline void
createRandomFile( const std::string& path,
                  const uint64_t     size )
{
    std::ofstream textFile( path, std::ios_base::out | std::ios_base::binary );

    std::mt19937_64 randomEngine;
    std::array<uint64_t, 4_Ki> buffer{};  // 32 KiB of buffer
    for ( size_t nBytesWritten = 0; nBytesWritten < size; ) {
        for ( auto& x : buffer ) {
            x = randomEngine();
        }
        const auto nBytesToWrite = std::min<uint64_t>( buffer.size() * sizeof( buffer[0] ), size - nBytesWritten );
        textFile.write( reinterpret_cast<const char*>( buffer.data() ), static_cast<std::streamsize>( nBytesToWrite ) );
        nBytesWritten += nBytesToWrite;
    }
}


template<typename Container>
void
fillWithRandomBase64( Container& container )
{
    for ( size_t i = 0; i < container.size(); ++i ) {
        container[i] = static_cast<typename Container::value_type>(
            ( i + 1 == container.size() ) || ( ( i + 1 ) % 77 == 0 )
            ? '\n' : BASE64_SYMBOLS[static_cast<size_t>( rand() ) % BASE64_SYMBOLS.size()] );
    }
}


[[nodiscard]] inline std::vector<char>
createRandomBase64( const size_t size )
{
    std::vector<char> result( size );
    fillWithRandomBase64( result );
    return result;
}


inline void
createRandomBase64( const std::string& filePath,
                    const size_t       fileSize )
{
    std::ofstream file{ filePath, std::ios_base::out | std::ios_base::binary };
    for ( size_t i = 0; i < fileSize; ++i ) {
        file << ( ( i + 1 == fileSize ) || ( ( i + 1 ) % 77 == 0 )
                  ? '\n' : BASE64_SYMBOLS[static_cast<size_t>( rand() ) % BASE64_SYMBOLS.size()] );
    }
}


template<typename Container>
void
fillWithRandomNumbers( Container& container )
{
    constexpr std::string_view DIGITS = "0123456789";
    for ( size_t i = 0; i < container.size(); ++i ) {
        container[i] = static_cast<typename Container::value_type>(
            ( i + 1 == container.size() ) || ( ( i + 1 ) % 77 == 0 )
            ? '\n' : DIGITS[static_cast<size_t>( rand() ) % DIGITS.size()] );
    }
}


[[nodiscard]] inline std::vector<char>
createRandomNumbers( const size_t size )
{
    std::vector<char> result( size );
    fillWithRandomNumbers( result );
    return result;
}


inline void
createRandomNumbers( const std::string& filePath,
                     const size_t       fileSize )
{
    constexpr std::string_view DIGITS = "0123456789";
    std::ofstream file{ filePath, std::ios_base::out | std::ios_base::binary };
    for ( size_t i = 0; i < fileSize; ++i ) {
        file << ( ( i + 1 == fileSize ) || ( ( i + 1 ) % 77 == 0 )
                  ? '\n' : DIGITS[static_cast<size_t>( rand() ) % DIGITS.size()] );
    }
}


inline void
createZeros( const std::string& filePath,
             const size_t       fileSize )
{
    std::ofstream file{ filePath, std::ios_base::out | std::ios_base::binary };
    static constexpr std::array<char, 4_Ki> BUFFER{};
    for ( size_t i = 0; i < fileSize; i += BUFFER.size() ) {
        const auto size = std::min( fileSize - i, BUFFER.size() );
        file.write( BUFFER.data(), size );
    }
}


inline void
createRandomWords( const std::string& filePath,
                   const size_t       fileSize )
{
    static constexpr size_t WORD_SIZE{ 16 };
    std::vector<std::array<char, WORD_SIZE> > words( 32 );
    for ( auto& word : words ) {
        for ( auto& c : word ) {
            c = static_cast<char>( rand() );
        }
    }

    std::ofstream file{ filePath, std::ios_base::out | std::ios_base::binary };
    for ( size_t i = 0; i < fileSize; ) {
        const auto iWord = static_cast<size_t>( rand() ) % words.size();
        file.write( words[iWord].data(), words[iWord].size() );
        i += words[iWord].size();
    }
}
}  // namespace rapidgzip
