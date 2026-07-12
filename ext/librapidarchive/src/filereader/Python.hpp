#pragma once

#include <algorithm>
#include <cstring>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <Python.h>

#include <core/ScopedGIL.hpp>

#include "FileReader.hpp"


namespace rapidgzip
{
template<typename Value,
         typename std::enable_if_t<std::is_integral_v<Value> && std::is_signed_v<Value>, void*> = nullptr>
[[nodiscard]] PyObject*
toPyObject( Value value )
{
    auto* const result = PyLong_FromLongLong( value );
    if ( result == nullptr ) {
        throw std::runtime_error( "PyLong_FromLongLong returned null for: " + std::to_string( value ) + "!" );
    }
    return result;
}


template<typename Value,
         typename std::enable_if_t<std::is_integral_v<Value> && !std::is_signed_v<Value>, void*> = nullptr>
[[nodiscard]] PyObject*
toPyObject( Value value )
{
    auto* const result = PyLong_FromUnsignedLongLong( value );
    if ( result == nullptr ) {
        throw std::runtime_error( "PyLong_FromUnsignedLongLong returned null for: " + std::to_string( value ) + "!" );
    }
    return result;
}


[[nodiscard]] PyObject*
toPyObject( PyObject* value )
{
    if ( value == nullptr ) {
        throw std::runtime_error( "Got null PyObject as argument to toPyObject!" );
    }
    return value;
}


template<typename Result>
[[nodiscard]] Result
fromPyObject( PyObject* pythonObject );

/* For "long" objects, see https://docs.python.org/3/c-api/long.html */

template<>
[[nodiscard]] long long int
fromPyObject<long long int>( PyObject* pythonObject )
{
    /** @todo How to do better error reporting if the object can't be converted to the specified type? */
    return PyLong_AsLongLong( pythonObject );
}


template<>
[[nodiscard]] size_t
fromPyObject<size_t>( PyObject* pythonObject )
{
    return PyLong_AsSize_t( pythonObject );
}

/* For bool objects, see https://docs.python.org/3/c-api/long.html */

template<>
[[nodiscard]] bool
fromPyObject<bool>( PyObject* pythonObject )
{
    return pythonObject == Py_True;
}


template<>
[[nodiscard]] PyObject*
fromPyObject<PyObject*>( PyObject* pythonObject )
{
    return pythonObject;
}


template<typename Result,
         typename... Args>
Result
callPyObject( PyObject* pythonObject,
              Args...   args )
{
    constexpr auto nArgs = sizeof...( Args );

    if ( pythonObject == nullptr ) {
        throw std::invalid_argument( "[callPyObject] Got null PyObject!" );
    }

    const ScopedGILLock gilLock;

    if constexpr ( std::is_same_v<Result, void> ) {
        PyObject_Call( pythonObject, PyTuple_Pack( nArgs, toPyObject( args )... ), nullptr );
    } else {
        const auto result = PyObject_Call( pythonObject, PyTuple_Pack( nArgs, toPyObject( args )... ), nullptr );
        if ( result == nullptr ) {
            std::stringstream message;
            message << "Cannot convert nullptr Python object to the requested result type ("
                    << typeid( Result ).name() << ")!";
            if ( ( pythonObject != nullptr ) && ( pythonObject->ob_type != nullptr ) ) {
                message << " Got no result when calling: " << pythonObject->ob_type->tp_name;
            }
            throw std::invalid_argument( std::move( message ).str() );
        }
        return fromPyObject<Result>( result );
    }
}


class PythonFileReader :
    public FileReader
{
public:
    explicit
    PythonFileReader( PyObject* pythonObject ) :
        m_pythonObject( checkNullptr( pythonObject ) ),

        mpo_tell( getAttribute( m_pythonObject, "tell" ) ),
        mpo_seek( getAttribute( m_pythonObject, "seek" ) ),
        mpo_read( getAttribute( m_pythonObject, "read" ) ),
        mpo_write( getAttribute( m_pythonObject, "write" ) ),
        mpo_seekable( getAttribute( m_pythonObject, "seekable" ) ),
        mpo_close( getAttribute( m_pythonObject, "close" ) ),

        m_initialPosition( callPyObject<long long int>( mpo_tell ) ),
        m_seekable( callPyObject<bool>( mpo_seekable ) )
    {
        if ( m_seekable ) {
            m_fileSizeBytes = seek( 0, SEEK_END );
            seek( 0, SEEK_SET );
        }

        /* No throws should come after incrementing the reference count because an exception thrown inside the
         * constructor means the destructor will not be called! */
        Py_INCREF( m_pythonObject );
    }

    ~PythonFileReader()
    {
        close();
    }

    [[nodiscard]] UniqueFileReader
    cloneRaw() const override
    {
        throw std::invalid_argument( "Cloning file path reader not allowed because the internal file position "
                                     "should not be modified by multiple owners!" );
    }

    /* Copying is simply not allowed because that might interfere with the file position state, use SharedFileReader! */

    void
    close() override
    {
        if ( m_pythonObject == nullptr ) {
            return;
        }

        /* Try to restore the file position the file had before it was given to us. */
        if ( m_seekable ) {
            seek( m_initialPosition );
        }

        const ScopedGILLock gilLock;
        if ( Py_REFCNT( m_pythonObject ) == 1 ) {
            callPyObject<void>( mpo_close );
        }
        Py_DECREF( m_pythonObject );
        m_pythonObject = nullptr;
    }

    [[nodiscard]] bool
    closed() const override
    {
        return m_pythonObject == nullptr;
    }

    [[nodiscard]] bool
    eof() const override
    {
        return m_seekable ? tell() >= size() : !m_lastReadSuccessful;
    }

    [[nodiscard]] bool
    fail() const override
    {
        return false;
    }

    [[nodiscard]] int
    fileno() const override
    {
        throw std::invalid_argument( "This Python file-like object has no valid fileno!" );
    }

    [[nodiscard]] bool
    seekable() const override
    {
        return m_seekable;
    }

    [[nodiscard]] size_t
    read( char*  buffer,
          size_t nMaxBytesToRead ) override
    {
        if ( m_pythonObject == nullptr ) {
            throw std::invalid_argument( "Invalid or file can't be read from!" );
        }

        if ( nMaxBytesToRead == 0 ) {
            return 0;
        }

        const ScopedGILLock gilLock;

        /** @todo better to use readinto because read might return less than requested even before the EOF! */
        auto* const bytes = callPyObject<PyObject*>( mpo_read, nMaxBytesToRead );
        if ( !PyBytes_Check( bytes ) ) {
            Py_XDECREF( bytes );
            throw std::runtime_error( "Expected a bytes object to be returned by read!" );
        }

        const auto nBytesRead = PyBytes_Size( bytes );
        if ( buffer != nullptr ) {
            std::memset( buffer, '\0', nBytesRead );
            std::memcpy( buffer, PyBytes_AsString( bytes ), nBytesRead );
        }
        Py_XDECREF( bytes );

        if ( nBytesRead < 0 ) {
            std::stringstream message;
            message
                << "[PythonFileReader] Read call failed (" << nBytesRead << " B read)!\n"
                << "  Buffer: " << (void*)buffer << "\n"
                << "  nMaxBytesToRead: " << nMaxBytesToRead << " B\n"
                << "  File size: " << m_fileSizeBytes << " B\n"
                << "  m_currentPosition: " << m_currentPosition << "\n"
                << "  tell: " << tell() << "\n"
                << "\n";
            std::cerr << message.str();
            throw std::domain_error( std::move( message ).str() );
        }

        m_currentPosition += nBytesRead;
        m_lastReadSuccessful = static_cast<size_t>( nBytesRead ) == nMaxBytesToRead;

        return nBytesRead;
    }

    /**
     * Should not be mixed with read calls on the same PythonFileReader object!
     */
    [[nodiscard]] size_t
    write( const char* buffer,
           size_t      nBytesToWrite )
    {
        if ( m_pythonObject == nullptr ) {
            throw std::invalid_argument( "Invalid or file can't be written to!" );
        }

        if ( nBytesToWrite == 0 ) {
            return 0;
        }

        const ScopedGILLock gilLock;

        auto* const bytes = PyBytes_FromStringAndSize( buffer, nBytesToWrite );
        const auto nBytesWritten = callPyObject<long long int>( mpo_write, bytes );

        if ( ( nBytesWritten < 0 ) || ( static_cast<size_t>( nBytesWritten ) < nBytesToWrite ) ) {
            std::stringstream message;
            message
                << "[PythonFileReader] Write call failed (" << nBytesWritten << " B written)!\n"
                << "  Buffer: " << (void*)buffer << "\n"
                << "  tell: " << tell() << "\n"
                << "\n";
            std::cerr << message.str();
            throw std::domain_error( std::move( message ).str() );
        }

        return static_cast<size_t>( nBytesWritten );
    }

    size_t
    seek( long long int offset,
          int           origin = SEEK_SET ) override
    {
        if ( ( m_pythonObject == nullptr ) || !m_seekable ) {
            throw std::invalid_argument( "Invalid or unseekable file can't be seeked!" );
        }

        /* see https://docs.python.org/3/library/io.html#io.IOBase.seek */
        int pythonWhence = 0;
        switch ( origin )
        {
        case SEEK_SET:
            pythonWhence = 0;
            break;
        case SEEK_CUR:
            pythonWhence = 1;
            break;
        case SEEK_END:
            pythonWhence = 2;
            break;
        }

        m_currentPosition = callPyObject<size_t>( mpo_seek, offset, pythonWhence );

        return m_currentPosition;
    }

    [[nodiscard]] std::optional<size_t>
    size() const override
    {
        return m_fileSizeBytes;
    }

    [[nodiscard]] size_t
    tell() const override
    {
        if ( m_pythonObject == nullptr ) {
            throw std::invalid_argument( "Cannot call 'tell' on invalid file!" );
        }
        return callPyObject<size_t>( mpo_tell );
    }

    void
    clearerr() override
    {
        /* There exists nothing like this for Python's RawIOBase */
    }

private:
    [[nodiscard]] static PyObject*
    checkNullptr( PyObject* pythonObject )
    {
        if ( pythonObject == nullptr ) {
            throw std::invalid_argument( "PythonFileReader may not be constructed with a nullptr PyObject!" );
        }
        return pythonObject;
    }

    [[nodiscard]] static PyObject*
    getAttribute( PyObject*   pythonObject,
                  const char* name )
    {
        auto* const attribute = PyObject_GetAttrString( pythonObject, name );

        if ( attribute == nullptr ) {
            std::stringstream message;
            message << "The given Python file-like object must have a '" << name << "' method!";
            throw std::invalid_argument( std::move( message ).str() );
        }

        return attribute;
    }

protected:
    PyObject* m_pythonObject{ nullptr };
    PyObject* const mpo_tell;
    PyObject* const mpo_seek;
    PyObject* const mpo_read;
    PyObject* const mpo_write;
    PyObject* const mpo_seekable;
    PyObject* const mpo_close;

    const long long int m_initialPosition;
    const bool m_seekable;
    size_t m_fileSizeBytes;

    size_t m_currentPosition{ 0 };  /**< Only necessary for unseekable files. */
    bool m_lastReadSuccessful{ true };
};
}  // namespace rapidgzip
