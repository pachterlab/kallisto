#pragma once

#include <stdexcept>
#include <string>


namespace rapidgzip
{
class DecompressionError :
    public std::runtime_error
{
public:
    DecompressionError( const std::string& message ) :
        std::runtime_error( message )
    {}
};


class NoBlockInRange :
    public DecompressionError
{
public:
    NoBlockInRange( const std::string& message ) :
        DecompressionError( message )
    {}
};
}  // namespace rapidgzip
