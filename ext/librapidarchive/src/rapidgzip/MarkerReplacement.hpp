#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <core/VectorView.hpp>
#include <rapidgzip/gzip/definitions.hpp>         // MAX_WINDOW_SIZE


namespace rapidgzip::deflate
{
template<bool FULL_WINDOW>
struct MapMarkers
{
    MapMarkers( VectorView<uint8_t> const& window ) :
        m_window( window )
    {
        assert( ( m_window.size() >= MAX_WINDOW_SIZE ) == FULL_WINDOW );
    }

    [[nodiscard]] constexpr uint8_t
    operator()( uint16_t value ) const
    {
        if ( value <= std::numeric_limits<uint8_t>::max() ) {
            return static_cast<uint8_t>( value );
        }

        if ( value < MAX_WINDOW_SIZE ) {
            throw std::invalid_argument( "Cannot replace unknown 2 B code!" );
        }

        if constexpr ( !FULL_WINDOW ) {
            if ( value - MAX_WINDOW_SIZE >= m_window.size() ) {
                throw std::invalid_argument( "Window too small!" );
            }
        }

        return m_window[value - MAX_WINDOW_SIZE];
    }

private:
    const VectorView<uint8_t> m_window;
};


inline void
replaceMarkerBytes( WeakVector<std::uint16_t>       buffer,
                    VectorView<std::uint8_t> const& window )
{
    /* For maximum size windows, we can skip one check because even UINT16_MAX is valid. */
    if ( window.size() >= MAX_WINDOW_SIZE ) {
        std::transform( buffer.begin(), buffer.end(), buffer.begin(), MapMarkers<true>( window ) );
    } else {
        std::transform( buffer.begin(), buffer.end(), buffer.begin(), MapMarkers<false>( window ) );
    }
}
}  // namespace rapidgzip::deflate
