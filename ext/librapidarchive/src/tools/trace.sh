#!/usr/bin/env bash

set -e

export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=1G
export SCOREP_EXPERIMENT_DIRECTORY=scorep-bzcat
export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=48
#export SCOREP_FILTERING_FILE=custom.scorep.filter # better to use --instrument-filter but might be required when using --nocompiler, however only the compiler plugin can isntrument inline functions
# These require Score-P being built with libunwind-dev!
# Does not work for me ... When I activate it, I only get I/O functions traced none of my own ...
#export SCOREP_ENABLE_UNWINDING=true # add calling context like source code location!
#export SCOREP_SAMPLING_EVENTS= # deactivate sampling (only unwind events)
#export SCOREP_SAMPLING_EVENTS=perf_cycles@100000
#export SCOREP_SAMPLING_EVENTS=PAPI_TOT_CYC@100000

#scorep --io=posix --instrument-filter=custom.scorep.filter g++ -I ../indexed_bzip2 -std=c++11 -Wall -Wextra -Wshadow -c bzcat.cpp -o bzcat.o
#scorep --io=posix --instrument-filter=custom.scorep.filter g++ bzcat.o -o bzcat

#scorep --io=posix g++ -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# only main an I/O calls

#cat <<EOF > custom.scorep.filter
#SCOREP_REGION_NAMES_BEGIN
#  INCLUDE *
#SCOREP_REGION_NAMES_END
#EOF
#scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# only main and I/O calls

#cat <<EOF > custom.scorep.filter
#SCOREP_REGION_NAMES_BEGIN
#  INCLUDE BZ2Reader*
#SCOREP_REGION_NAMES_END
#EOF
#scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# only BZ2Reader::BlockHeader BZ2Reader::readBlockHeader and BZ2Reader::BZ2Reader and I/O calls
#  -> looks like return typ is indeed matched in these filter rules

#cat <<EOF > custom.scorep.filter
#SCOREP_REGION_NAMES_BEGIN
#  INCLUDE *BZ2Reader*
#SCOREP_REGION_NAMES_END
#EOF
#scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# all BZ2Reader functions and I/O calls! Considerably slower though
# 40% of time is in BZ2Reader::getBits

#cat <<EOF > custom.scorep.filter
#SCOREP_REGION_NAMES_BEGIN
#  INCLUDE *BZ2Reader*
#  INCLUDE *BitReader*
#SCOREP_REGION_NAMES_END
#EOF
#scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# 20% in BitReader::read, 42% in getBits, 22.5% readBlockData

cat <<EOF > custom.scorep.filter
SCOREP_REGION_NAMES_BEGIN
  INCLUDE *BZ2Reader*
  INCLUDE *BitReader*
  EXCLUDE *getBits*
  EXCLUDE *BitReader::read*
SCOREP_REGION_NAMES_END
EOF
#scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# 38% write (to stdout), 26% decodeStream, 25% readBlockData (for 4k read buffers)

#scorep --nocompiler --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat
# only I/O calls are traced, not even main. Custom functions are not traced because they are inlined but I can't explain the missing main

scorep --instrument-filter=custom.scorep.filter --io=posix g++ -O3 -I ../indexed_bzip2 -std=c++11 bzcat.cpp -o bzcat

#time ./bzcat test.bz2 4096 | wc -c
#  => many small calls to write and so on are very bad for performance!
#  => what kind of calls is FUSE doing?
#  => 0.46s
time ./bzcat test-100k.bz2 $(( 1024 * 1024 )) | wc -c
#  => 0.18s
# all measurements are including tracing overhead
#  - after increasing internal buffers from 4kiB to 1MiB => 0.13s, ~42% readBlockData, ~42% decodeStream
#

scorep-score -r scorep-bzcat/profile.cubex
