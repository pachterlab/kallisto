#pragma once

#include <string>

#include <cxxopts.hpp>


[[nodiscard]] inline std::string
getFilePath( cxxopts::ParseResult const& parsedArgs,
             std::string          const& argument )
{
    if ( parsedArgs.count( argument ) > 1 ) {
        if ( parsedArgs.count( "quiet" ) == 0 ) {
            std::cerr << "[Warning] Multiple output files specified. Will only use the last one: "
                      << parsedArgs[argument].as<std::string>() << "!\n";
        }
    }
    if ( parsedArgs.count( argument ) > 0 ) {
        auto path = parsedArgs[argument].as<std::string>();
        if ( path != "-" ) {
            return path;
        }
    }
    return {};
}
