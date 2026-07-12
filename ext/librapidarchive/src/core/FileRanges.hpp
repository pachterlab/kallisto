#pragma once

#include <charconv>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "common.hpp"

namespace rapidgzip
{
struct FileRange
{
    size_t offset{ 0 };
    size_t size{ 0 };
    bool offsetIsLine{ false };
    bool sizeIsLine{ false };

    [[nodiscard]] bool
    operator==( const FileRange& other ) const noexcept
    {
        return ( offset == other.offset ) && ( size == other.size )
               && ( offsetIsLine == other.offsetIsLine ) && ( sizeIsLine == other.sizeIsLine );
    }
};


inline std::ostream&
operator<<( std::ostream&    out,
            const FileRange& range )
{
    out << range.size << ( range.sizeIsLine ? "L" : "" ) << "@" << range.offset << ( range.offsetIsLine ? "L" : "" );
    return out;
}


[[nodiscard]] inline const char*
skipWhitespaces( const char* const first,
                 const char* const last )
{
    constexpr std::string_view WHITESPACES{ " \t" };
    const std::string_view view{ first, static_cast<size_t>( std::distance( first, last ) ) };
    const auto skippedCount = view.find_first_not_of( WHITESPACES );
    if ( skippedCount < view.size() ) {
        return first + skippedCount;
    }
    return last;
}


[[nodiscard]] inline const char*
readNumber( const char* const first,
            const char* const last,
            size_t&           value,
            bool&             valueIsLine )
{
    const auto result = std::from_chars( first, last, value );
    if ( ( result.ec != std::errc() ) || ( result.ptr == first ) ) {
        throw std::invalid_argument( "Failed to parse number at the start of the remaining expression: "
                                     + std::string( first, last ) );
    }

    static const std::vector<std::pair<std::string_view, size_t> > PREFIXES{
        { "ki", 1ULL << 10U },  /* Strily speaking neither valid SI or IEC, but I like it for consistency with SI. */
        { "Ki", 1ULL << 10U },
        { "Mi", 1ULL << 20U },
        { "Gi", 1ULL << 30U },
        { "Ti", 1ULL << 40U },
        { "Pi", 1ULL << 50U },
        { "Ei", 1ULL << 60U },
        { "k", 1000ULL },
        { "M", 1000'000ULL },
        { "G", 1000'000'000ULL },
        { "T", 1000'000'000'000ULL },
        { "P", 1000'000'000'000'000ULL },
        { "E", 1000'000'000'000'000'000ULL },
        { "", 1ULL },  /* No prefix special case. Required to detect "B" suffix. */
    };

    const auto* current = skipWhitespaces( result.ptr, last );
    const std::string_view unitString{ current, static_cast<size_t>( std::distance( current, last ) ) };

    size_t longestMatch{ 0 };
    size_t longestMatchFactor{ 1 };
    bool longestMatchIsLine{ false };
    for ( const auto suffix : { std::string_view( "B" ), std::string_view( "L" ), std::string_view( "" ) } ) {
        for ( const auto& [prefix, factor] : PREFIXES ) {
            if ( ( prefix.size() + suffix.size() > longestMatch )
                 && startsWith( unitString, prefix )
                 && startsWith( std::string_view( unitString.data() + prefix.size(),
                                                  unitString.size() - prefix.size() ),
                                suffix ) )
            {
                longestMatch = prefix.size() + suffix.size();
                longestMatchFactor = factor;
                longestMatchIsLine = suffix == "L";
            }
        }
    }

    valueIsLine = longestMatchIsLine;
    if ( longestMatch == 0 ) {
        return result.ptr;
    }
    value *= longestMatchFactor;
    return current + longestMatch;
}


[[nodiscard]] inline std::vector<FileRange>
parseFileRanges( const std::string_view expression )
{
    constexpr char OFFSET_PREFIX{ '@' };
    constexpr char SEPARATOR{ ',' };
    static constexpr std::string_view INFINITY_STRING{ "inf" };

    /**
     * Allowed state transitions:
     * @verbatim
     *           23           @                      10            ,
     * TUPLE_END -> SIZE_END -> OFFSET_SEPARATOR_END -> OFFSET_END -> TUPLE_END
     *   +---^ (simply ignore duplicated tuple separators: ",,,")
     * @endverbatim
     * Any whitespace between these states is ignored. Empty sizes or offsets are not allowed, i.e.,
     * SIZE_END may not skip to TUPLE_END and so on.
     * I introduced the state machine to avoid duplicate end-of-string checks. Thanks to the state machine,
     * we can simply enter the next loop iteration, which checks for end-of-string, after reading each part
     * or skipping over whitespace characters.
     */
    enum class State
    {
        TUPLE_END,
        SIZE_END,
        OFFSET_SEPARATOR_END,
        OFFSET_END,
    };

    std::vector<FileRange> ranges;

    const auto* const expressionEnd = expression.data() + expression.size();
    auto state = State::TUPLE_END;
    FileRange range;

    const auto* current = skipWhitespaces( expression.data(), expressionEnd );
    for ( ; current != expressionEnd; current = skipWhitespaces( current, expressionEnd ) ) {
        switch ( state )
        {
        case State::TUPLE_END:
            if ( *current == SEPARATOR ) {
                ++current;
                continue;
            }

            range.size = 0;
            if ( startsWith( std::string_view( current,
                                               static_cast<size_t>( std::distance( current, expressionEnd ) ) ),
                             INFINITY_STRING ) )
            {
                range.size = std::numeric_limits<size_t>::max();
                current += INFINITY_STRING.size();
            } else {
                current = readNumber( current, expressionEnd, range.size, range.sizeIsLine );
            }
            state = State::SIZE_END;
            break;

        case State::SIZE_END:
            if ( *current != OFFSET_PREFIX ) {
                std::stringstream message;
                message << "Expected " << OFFSET_PREFIX << " after a size at position "
                        << std::distance( expression.data(), current ) << " in expression: " << expression;
                throw std::invalid_argument( std::move( message ).str() );
            }
            state = State::OFFSET_SEPARATOR_END;
            ++current;
            break;

        case State::OFFSET_SEPARATOR_END:
            range.offset = 0;
            current = readNumber( current, expressionEnd, range.offset, range.offsetIsLine );
            ranges.emplace_back( range );
            state = State::OFFSET_END;
            break;

        case State::OFFSET_END:
            if ( *current != SEPARATOR ) {
                std::stringstream message;
                message << "Expected " << SEPARATOR << " after an size@offset tuple at position "
                        << std::distance( expression.data(), current ) << " in expression: " << expression;
                throw std::invalid_argument( std::move( message ).str() );
            }
            ++current;
            state = State::TUPLE_END;
            break;
        }
    }

    if ( ( state != State::TUPLE_END ) && ( state != State::OFFSET_END ) ) {
        throw std::invalid_argument( "Incomplete size@offset tuple at end of expression: "
                                     + std::string( expression ) );
    }

    return ranges;
}
}  // namespace rapidgzip
