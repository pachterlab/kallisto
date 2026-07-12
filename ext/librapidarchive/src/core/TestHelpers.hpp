#pragma once

#include <iostream>
#include <filesystem>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "common.hpp"  // duration

#ifdef __APPLE_CC__
    #include <AvailabilityMacros.h>
#endif


namespace rapidgzip
{
int gnTests = 0;  // NOLINT
int gnTestErrors = 0;  // NOLINT


template<typename A,
         typename B>
void
requireEqual( const A&  a,
              const B&  b,
              const int line )
{
    ++gnTests;
    if ( a != b ) {
        ++gnTestErrors;
        std::cerr << "[FAIL on line " << line << "] " << a << " != " << b << "\n";
    }
}


void
require( bool               condition,
         std::string const& conditionString,
         int                line )
{
    ++gnTests;
    if ( !condition ) {
        ++gnTestErrors;
        std::cerr << "[FAIL on line " << line << "] " << conditionString << "\n";
    }
}


#define REQUIRE_EQUAL( a, b ) requireEqual( a, b, __LINE__ )  // NOLINT
#define REQUIRE( condition ) require( condition, #condition, __LINE__ )  // NOLINT
#define REQUIRE_THROWS( condition ) require( [&] () { \
    try { \
        (void)condition; \
    } catch ( const std::exception& ) { \
        return true; \
    } \
    return false; \
} (), #condition, __LINE__ )  // NOLINT


template<
    size_t   REPETITIONS,
    typename Functor,
    typename FunctorResult = std::invoke_result_t<Functor>
>
[[nodiscard]] std::pair<FunctorResult, std::vector<double> >
benchmarkFunction( Functor functor )
{
    std::optional<FunctorResult> result;
    std::vector<double> durations;
    for ( size_t i = 0; i < REPETITIONS; ++i ) {
        const auto t0 = now();
        const auto currentResult = functor();
        const auto t1 = now();
        durations.push_back( duration( t0, t1 ) );

        if ( !result ) {
            result = std::move( currentResult );
        } else if ( *result != currentResult ) {
            throw std::logic_error( "Function to benchmark returns indeterministic results!" );
        }
    }

    return { *result, durations };
}


template<
    size_t   REPETITIONS,
    typename Functor,
    typename SetupFunctor,
    typename FunctorResult = std::invoke_result_t<Functor, std::invoke_result_t<SetupFunctor> >
>
[[nodiscard]] std::pair<FunctorResult, std::vector<double> >
benchmarkFunction( SetupFunctor setup,
                   Functor      functor )
{
    decltype( setup() ) setupResult;
    try {
        setupResult = setup();
    } catch ( const std::exception& e ) {
        std::cerr << "Failed to run setup with exception: " << e.what() << "\n";
        return {};
    }

    std::optional<FunctorResult> result;
    std::vector<double> durations;
    for ( size_t i = 0; i < REPETITIONS; ++i ) {
        const auto t0 = now();
        auto currentResult = functor( setupResult );
        const auto t1 = now();
        durations.push_back( duration( t0, t1 ) );

        if ( !result ) {
            result = std::move( currentResult );
        } else if ( *result != currentResult ) {
            throw std::logic_error( "Function to benchmark returns indeterministic results!" );
        }
    }

    return { *result, durations };
}


class ThreadSafeStreamBuffer:
    public std::stringbuf
{
protected:
    /* Put area */

    std::streamsize
    xsputn( char_type const* s,
            std::streamsize  count ) override
    {
        const std::scoped_lock lock( m_mutex );
        return std::stringbuf::xsputn( s, count );
    }

    int
    overflow( int c ) override
    {
        const std::scoped_lock lock( m_mutex );
        return std::stringbuf::overflow( c );
    }

    /* Get area */

    std::streamsize
    showmanyc() override
    {
        throw std::logic_error( "Not supported!" );
    }

    int_type
    underflow() override
    {
        throw std::logic_error( "Not supported!" );
    }

    int_type
    uflow() override
    {
        throw std::logic_error( "Not supported!" );
    }

    std::streamsize
    xsgetn( char_type*      /* s */,
            std::streamsize /* count */ ) override
    {
        throw std::logic_error( "Not supported!" );
    }

    /* Locales */

    void
    imbue( const std::locale& /* loc */ ) override
    {
        throw std::logic_error( "Not supported!" );
    }

    /* Positioning */

    std::stringbuf*
    setbuf( char_type*      /* s */,
            std::streamsize /* n */ ) override
    {
        throw std::logic_error( "Not supported!" );
    }

    pos_type
    seekoff( off_type                /* off */,
             std::ios_base::seekdir  /* dir */,
             std::ios_base::openmode /* which */ = std::ios_base::in | std::ios_base::out ) override
    {
        throw std::logic_error( "Not supported!" );
    }

    pos_type
    seekpos( pos_type                /* pos */,
             std::ios_base::openmode /* which */ = std::ios_base::in | std::ios_base::out ) override
    {
        throw std::logic_error( "Not supported!" );
    }

    int
    sync() override
    {
        const std::scoped_lock lock( m_mutex );
        return std::stringbuf::sync();
    }

    /* Putback */

    int_type
    pbackfail( int_type /* c */ ) override
    {
        throw std::logic_error( "Not supported!" );
    }

private:
    std::recursive_mutex m_mutex;
};


class StreamInterceptor :
    public ThreadSafeStreamBuffer
{
public:
    explicit
    StreamInterceptor( std::ostream& out ) :
        m_out( out ),
        m_rdbuf( m_out.rdbuf( this ) )
    {}

    ~StreamInterceptor()
    {
        close();
    }

    void
    close()
    {
        if ( m_rdbuf.has_value() ) {
            /* Expected to return this and therefore can be ignored. */
            m_out.rdbuf( *m_rdbuf );
            m_rdbuf.reset();
        }
    }

    StreamInterceptor( const StreamInterceptor& ) = delete;
    StreamInterceptor( StreamInterceptor&& ) = delete;
    StreamInterceptor& operator=( const StreamInterceptor& ) = delete;
    StreamInterceptor& operator=( StreamInterceptor&& ) = delete;

private:
    std::ostream& m_out;
    std::optional<std::basic_streambuf<char>*> m_rdbuf;
};


/* error: 'std::filesystem::path' is unavailable: introduced in macOS 10.15.
 * Fortunately, this is only needed for the tests, so the incomplete std::filesystem support
 * is not a problem for building the manylinux wheels on the pre 10.15 macOS kernel.
 * https://opensource.apple.com/source/xnu/xnu-2050.7.9/EXTERNAL_HEADERS/AvailabilityMacros.h.auto.html */
#if !defined(__APPLE_CC__ ) || ( defined(MAC_OS_X_VERSION_MIN_REQUIRED) \
    && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_15 )
class TemporaryDirectory
{
public:
    explicit
    TemporaryDirectory( std::filesystem::path path ) :
        m_path( std::move( path ) )
    {}

    TemporaryDirectory( TemporaryDirectory&& ) = default;

    TemporaryDirectory( const TemporaryDirectory& ) = delete;

    TemporaryDirectory&
    operator=( TemporaryDirectory&& ) = default;

    TemporaryDirectory&
    operator=( const TemporaryDirectory& ) = delete;

    ~TemporaryDirectory()
    {
        if ( !m_path.empty() ) {
            std::filesystem::remove_all( m_path );
        }
    }

    [[nodiscard]] operator std::filesystem::path() const
    {
        return m_path;
    }

    [[nodiscard]] const std::filesystem::path&
    path() const
    {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};


[[nodiscard]] inline TemporaryDirectory
createTemporaryDirectory( const std::string& title = "tmpTest" )
{
    const std::filesystem::path tmpFolderName = title + "." + std::to_string( unixTimeInNanoseconds() );
    std::filesystem::create_directory( tmpFolderName );
    return TemporaryDirectory( tmpFolderName );
}
#endif
}  // namespace rapidgzip
