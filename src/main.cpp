#include <cstdio>

#include <string>
#include <iostream>

#include "example.h"

void usage()
{
    // TODO: write me !
    std::cerr << "Usage: kallisto TODO" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        usage();
        return 1;
    }
    std::cout << "Hello world!" << std::endl;
    std::cout << sq(3.0) << std::endl;

    return 0;
}
