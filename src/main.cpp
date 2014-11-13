#include <cstdio>

#include <string>
#include <iostream>

#include "example.h"

#include "gff.h"

void usage()
{
    std::cerr << "Usage: xprs_wiggle output.bam" << std::endl;
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
    const std::string fname = "/Users/hjp/lmcb/intron-retention/src/tests/inputs/refGene_07.23.2014_CHL1.gtf";
    FILE* fhandle = fopen(fname.c_str(), "r");

    GffReader gff_reader;
    gff_reader.init(fhandle, true);

    std::cout << "here!" << std::endl;

    gff_reader.readAll();

    std::cout << "there!" << std::endl;

    fclose(fhandle);

    return 0;
}
