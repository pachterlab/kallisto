
#include <iostream>
#include <fstream>
#include <string_view>
#include <vector>

#include <igzip_lib.h>


std::vector<char>
readFile( const std::string& fileName )
{
    std::ifstream file( fileName );

    file.seekg( 0, std::ios_base::end );
    const auto fileSize = file.tellg();
    if ( fileSize < 0 ) {
        throw std::invalid_argument( "Could not get size of specified file!" );
    }

    std::vector<char> buffer( fileSize );
    file.seekg( 0 );
    file.read( buffer.data(), buffer.size() );
    return buffer;
}


std::vector<char>
compressWithISAL( const std::vector<char>& toCompress )
{
    /* > Parameter avail_out must be large enough to fit the entire compressed output.
     * > Max expansion is limited to the input size plus the header size of a stored/raw block. */
    std::vector<char> compressed( toCompress.size() + 1000 );

    isal_zstream stream;
    isal_deflate_stateless_init( &stream );
    stream.next_in = reinterpret_cast<uint8_t*>( const_cast<char*>( toCompress.data() ) );
    stream.avail_in = toCompress.size();
    stream.next_out = reinterpret_cast<uint8_t*>( const_cast<char*>( compressed.data() ) );
    stream.avail_out = compressed.size();
    stream.gzip_flag = IGZIP_GZIP;
    const auto result = isal_deflate_stateless( &stream );
    if ( result != ISAL_DECOMP_OK ) {
        throw std::runtime_error( "Compression failed with error code: " + std::to_string( result ) );
    }
    if ( stream.avail_out >= compressed.size() ) {
        throw std::logic_error( "Something went wrong. Avail_out should be smaller than it was before." );
    }
    compressed.resize( compressed.size() - stream.avail_out );

    return compressed;
}


int main( int    argc,
          char** argv )
{
    if ( argc <= 1 ) {
        std::cerr << "Please specify a file to decompress.\n";
        return 1;
    }

    const auto toCompress = readFile( argv[1] );
    const auto compressed = compressWithISAL( toCompress );
    std::cerr << "Compressed " << toCompress.size() << " B down to " << compressed.size() << " B\n";
    std::fwrite( compressed.data(), compressed.size(), 1, stdout );
    return 0;
}
