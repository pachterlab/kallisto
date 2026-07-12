#pragma once

#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "FasterVector.hpp"


namespace rapidgzip
{
/**
 * Views are by their name read-only. This represents a read-only non-owned memory chunk.
 * Probably could be removed in favor of C++20 std::span.
 */
template<typename T>
class VectorView
{
public:
    using value_type = T;

public:
    constexpr
    VectorView() noexcept = default;

    constexpr
    VectorView( const VectorView<T>& ) noexcept = default;

    constexpr
    VectorView( VectorView<T>&& ) noexcept = default;

    constexpr VectorView<T>&
    operator=( const VectorView<T>& ) noexcept = default;

    constexpr VectorView<T>&
    operator=( VectorView<T>&& ) noexcept = default;

    template<typename Container,
             std::enable_if_t<std::is_same_v<Container, std::vector<T> >
                              || std::is_same_v<Container, FasterVector<T> >
                              || ( ( std::is_same_v<T, uint8_t>
                                     || std::is_same_v<T, char>
                                     || std::is_same_v<T, std::byte> ) &&
                                   ( std::is_same_v<typename Container::value_type, uint8_t>
                                     || std::is_same_v<typename Container::value_type, char>
                                     || std::is_same_v<typename Container::value_type, std::byte> ) )
             >* = nullptr>
    constexpr
    VectorView( const Container& vector ) noexcept :  // NOLINT
        m_data( reinterpret_cast<const T*>( vector.data() ) ),
        m_size( vector.size() )
    {}

    constexpr
    VectorView( const T* data,
                size_t   size ) noexcept :
        m_data( data ),
        m_size( size )
    {}

    constexpr
    VectorView( const T* data,
                const T* dataEnd ) noexcept :
        m_data( data ),
        m_size( std::distance( data, dataEnd ) )
    {}

    [[nodiscard]] constexpr T
    front() const noexcept
    {
        return *m_data;
    }

    [[nodiscard]] constexpr const T*
    begin() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const T*
    end() const noexcept
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr const T*
    data() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr size_t
    size() const noexcept
    {
        return m_size;
    }

    [[nodiscard]] constexpr bool
    empty() const noexcept
    {
        return m_size == 0;
    }

    [[nodiscard]] constexpr T
    operator[]( size_t i ) const noexcept
    {
        return m_data[i];
    }

    [[nodiscard]] constexpr T
    at( size_t i ) const
    {
        if ( i >= m_size ) {
            throw std::out_of_range( "VectorView index larger than size!" );
        }
        return m_data[i];
    }

    [[nodiscard]] explicit
    operator std::vector<T>() const
    {
        return std::vector<T>( begin(), end() );
    }

private:
    const T* m_data{ nullptr };
    size_t   m_size{ 0 };
};


/**
 * This is in basically a modifiable VectorView. The vector cannot be resized but
 * it is a simple pointer and length tuple to a read-write memory chunk.
 */
template<typename T>
class WeakVector
{
public:
    using value_type = T;

public:
    constexpr
    WeakVector() = default;

    constexpr
    WeakVector( const WeakVector& ) = default;

    constexpr
    WeakVector( WeakVector&& ) noexcept = default;

    constexpr WeakVector&
    operator=( const WeakVector& ) = default;

    constexpr WeakVector&
    operator=( WeakVector&& ) noexcept = default;

    constexpr
    WeakVector( std::vector<T>* vector ) noexcept :  // NOLINT
        m_data( vector->data() ),
        m_size( vector->size() )
    {}

    constexpr
    WeakVector( T*     data,
                size_t size ) noexcept :
        m_data( data ),
        m_size( size )
    {}

    [[nodiscard]] constexpr T
    front() const noexcept
    {
        return *m_data;
    }

    [[nodiscard]] constexpr const T*
    begin() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr T*
    begin() noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const T*
    end() const noexcept
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr T*
    end() noexcept
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr T*
    data() noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const T*
    data() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr size_t
    size() const noexcept
    {
        return m_size;
    }

    [[nodiscard]] constexpr bool
    empty() const noexcept
    {
        return m_size == 0;
    }

    [[nodiscard]] constexpr const T&
    operator[]( size_t i ) const noexcept
    {
        return m_data[i];
    }

    [[nodiscard]] constexpr T&
    operator[]( size_t i ) noexcept
    {
        return m_data[i];
    }

    [[nodiscard]] constexpr T
    at( size_t i ) const
    {
        if ( i >= m_size ) {
            throw std::out_of_range( "VectorView index larger than size!" );
        }
        return m_data[i];
    }

private:
    T* m_data{ nullptr };
    size_t m_size{ 0 };
};


template<typename T,
         size_t   T_size>
class ArrayView
{
public:
    using value_type = T;

public:
    explicit constexpr
    ArrayView( const std::array<T, T_size>& array ) noexcept :
        m_data( array.data() )
    {}

    [[nodiscard]] constexpr T
    front() const noexcept
    {
        return *m_data;
    }

    [[nodiscard]] constexpr const T*
    begin() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const T*
    end() const noexcept
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr const T*
    data() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr size_t
    size() const noexcept
    {
        return m_size;
    }

    [[nodiscard]] constexpr bool
    empty() const noexcept
    {
        return m_size == 0;
    }

    [[nodiscard]] constexpr T
    operator[]( size_t i ) const noexcept
    {
        return m_data[i];
    }

    [[nodiscard]] T
    at( size_t i ) const
    {
        if ( i >= m_size ) {
            throw std::out_of_range( "VectorView index larger than size!" );
        }
        return m_data[i];
    }

private:
    const T* const m_data;
    constexpr static size_t m_size = T_size;
};


template<typename T,
         size_t   T_size>
class WeakArray
{
public:
    using value_type = T;

public:
    constexpr explicit
    WeakArray( T* data ) noexcept :
        m_data( data )
    {}

    [[nodiscard]] constexpr T
    front() const noexcept
    {
        return *m_data;
    }

    [[nodiscard]] constexpr const T*
    begin() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const T*
    end() const noexcept
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr T*
    data() const noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr size_t
    size() const noexcept
    {
        return m_size;
    }

    [[nodiscard]] constexpr bool
    empty() const noexcept
    {
        return m_size == 0;
    }

    [[nodiscard]] constexpr T&
    operator[]( size_t i ) const noexcept
    {
        return m_data[i];
    }

    [[nodiscard]] constexpr T
    at( size_t i ) const
    {
        if ( i >= m_size ) {
            throw std::out_of_range( "VectorView index larger than size!" );
        }
        return m_data[i];
    }

private:
    T* const m_data;
    constexpr static size_t m_size = T_size;
};
}  // namespace rapidgzip
