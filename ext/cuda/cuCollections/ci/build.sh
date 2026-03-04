#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2023 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -eo pipefail

ORIGINAL_DIR=$(pwd)

resolve_path() {
    local input_path=$1
    # Check if the input is an absolute path
    if [[ "$input_path" = /* ]]; then
        echo "$input_path"
    else
        # Treat as a relative path or executable name
        # Check if it's in the PATH
        if command -v "$input_path" >/dev/null 2>&1; then
            echo "$input_path"
        else
            echo "$ORIGINAL_DIR/$input_path"
        fi
    fi
}

# Ensure the script is being executed in its containing directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";

# Determine repo root as the parent of the `ci` directory
REPO_ROOT="$(cd .. && pwd)"

# Script defaults
BUILD_TESTS=${BUILD_TESTS:-OFF}
BUILD_EXAMPLES=${BUILD_EXAMPLES:-OFF}
BUILD_BENCHMARKS=${BUILD_BENCHMARKS:-OFF}
CLEAN_BUILD=0 # Re-use existing artifacts by-default
BUILD_PREFIX=../build # <repo_root>/build
BUILD_INFIX=local # <repo_root>/build/local
DEBUG_BUILD=0 # Default build type is Release
PARALLEL_LEVEL=${PARALLEL_LEVEL:-$(nproc)} # defaults to number of cores in the system
CUDA_COMPILER=${CUDACXX:-nvcc} # $CUDACXX if set, otherwise `nvcc`
HOST_COMPILER=${CXX:-g++} # $CXX if set, otherwise `g++`
CUDA_ARCHS=native # detect system's GPU architectures
CXX_STANDARD=17

# Initialize CMAKE_ARGS from environment variable if available
if [ -n "${CMAKE_ARGS:-}" ]; then
    read -ra CMAKE_ARGS <<< "$CMAKE_ARGS"
else
    CMAKE_ARGS=()
fi

function usage {
    echo "cuCollections build script"
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -t/--tests: Build tests"
    echo "  -e/--examples: Build examples"
    echo "  -b/--benchmarks: Build benchmarks"
    echo "  -c/--clean: Clean (re-)build"
    echo "  --prefix: Build directory prefix (Defaults to <repo_root>/build)"
    echo "  -i/--infix: Build directory infix (Defaults to local)"
    echo "  -d/--debug: Debug build"
    echo "  -p/--parallel: Build parallelism (Defaults to \$PARALLEL_LEVEL if set, otherwise the system's number of CPU cores)"
    echo "  --cuda: CUDA compiler (Defaults to \$CUDACXX if set, otherwise nvcc)"
    echo "  --cxx: Host compiler (Defaults to \$CXX if set, otherwise g++)"
    echo "  --arch: Target CUDA arches, e.g. \"60-real;70;80-virtual\" (Defaults to the system's native GPU archs)"
    echo "  --std: CUDA/C++ standard (Defaults to 17)"
    echo "  -v/-verbose/--verbose: Enable shell echo for debugging"
    echo "  -h/-help/--help: Show this usage message"
    echo
    echo "Examples:"
    echo "  Basic Build:"
    echo "    $ $0"
    echo "    Runs a basic build with default settings, i.e., builds tests, examples, and benchmarks."
    echo "    Build files will be written to <repo_root>/build/local and symlinked to <repo_root>/build/latest."
    echo
    echo "  Custom Build Infix Directory:"
    echo "    $ $0 -i my_build"
    echo "    Build files will be written to the <repo_root>/build/my_build directory and symlinked to <repo_root>/build/latest."
    echo
    echo "  Parallel Build with Specific CUDA Architecture and CUDA Compiler:"
    echo "    $ PARALLEL_LEVEL=8 $0 --cuda /my_cuda_compiler/nvcc --arch 70;80"
    echo "    $ $0 -p 8 --cuda /my_cuda_compiler/nvcc --arch 70;80"
    echo "    Specifies parallel build level of 8 and CUDA architecture 70 and 80 with the specified CUDA compiler."
    echo "    Build files will be written to <repo_root>/build/local and symlinked to <repo_root>/build/latest."
    echo
    echo "  Debug Build with Tests and Examples:"
    echo "    $ CXX=g++-9 $0 -t -e -d"
    echo "    $ $0 --cxx g++-9 -t -e -d"
    echo "    Sets the host compiler to g++-9, builds tests and examples, and enables debug mode."
    echo "    Build files will be written to <repo_root>/build/local and symlinked to <repo_root>/build/latest."
    echo
    echo "  Custom Build Directory with Benchmarks:"
    echo "    $ BUILD_BENCHMARKS=ON $0 --prefix /custom/build --infix my_build"
    echo "    $ $0 --prefix /custom/build --infix my_build -b"
    echo "    Builds benchmarks only."
    echo "    Build files will be written to /custom/build/my_build and symlinked to /custom/build/latest."
    echo
    echo "  Verbose Mode for Debugging:"
    echo "    $ $0 -v --std 17"
    echo "    Enables verbose mode for detailed output and builds with C++17 standard."
    echo "    Build files will be written to <repo_root>/build/local and symlinked to <repo_root>/build/latest."
    echo
    echo "  Using CMAKE_ARGS Environment Variable:"
    echo "    $ CMAKE_ARGS=\"-DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_FEATURE=ON\" $0 -t"
    echo "    $ export CMAKE_ARGS=\"-DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_FEATURE=ON\""
    echo "    $ $0 -t"
    echo "    Uses CMAKE_ARGS environment variable to pass additional CMake options."
    echo "    Can be overridden by using -- followed by specific arguments."
    echo
    echo "  Pass-through to CMake:"
    echo "    -- [CMake args...]  Anything after -- is forwarded to CMake (overrides CMAKE_ARGS env var)"
    echo
    exit 1
}

# Parse options

# Copy the args into a temporary array, since we will modify them and
# the parent script may still need them.
args=("$@")
while [ "${#args[@]}" -ne 0 ]; do
    case "${args[0]}" in
    -t | --tests) BUILD_TESTS=ON; args=("${args[@]:1}");;
    -e | --examples) BUILD_EXAMPLES=ON; args=("${args[@]:1}");;
    -b | --benchmarks) BUILD_BENCHMARKS=ON; args=("${args[@]:1}");;
    -c | --clean) CLEAN_BUILD=1; args=("${args[@]:1}");;
    --prefix) BUILD_PREFIX=$(resolve_path "${args[1]}"); args=("${args[@]:2}");;
    -i | --infix) BUILD_INFIX="${args[1]}"; args=("${args[@]:2}");;
    -d | --debug) DEBUG_BUILD=1; args=("${args[@]:1}");;
    -p | --parallel)  PARALLEL_LEVEL="${args[1]}";  args=("${args[@]:2}");;
    --cuda) CUDA_COMPILER=$(resolve_path "${args[1]}"); args=("${args[@]:2}");;
    --cxx)  HOST_COMPILER=$(resolve_path "${args[1]}"); args=("${args[@]:2}");;
    --arch) CUDA_ARCHS="${args[1]}";    args=("${args[@]:2}");;
    --std)  CXX_STANDARD="${args[1]}";  args=("${args[@]:2}");;
    -v | -verbose | --verbose) VERBOSE=1; args=("${args[@]:1}");;
    --) CMAKE_ARGS=("${args[@]:1}"); break;;
    -h | -help | --help) usage ;;
    *) echo "Unrecognized option: ${args[0]}"; usage ;;
    esac
done

if [ $VERBOSE ]; then
    set -x
fi

# Convert to full paths:
HOST_COMPILER=$(which ${HOST_COMPILER})
CUDA_COMPILER=$(which ${CUDA_COMPILER})
# Make CUDA arch list compatible with cmake
CUDA_ARCHS=$(echo "$CUDA_ARCHS" | tr ' ,' ';;')

# Begin processing unsets after option parsing
set -u

if [ "$BUILD_INFIX" = "latest" ] || [ -z "$BUILD_INFIX" ]; then
    echo "Error: BUILD_INFIX cannot be empty or 'latest'" >&2
    exit 1
fi

# If no build target is specified, build all targets
if [ "$BUILD_TESTS" == "OFF" ] && [ "$BUILD_EXAMPLES" == "OFF" ] && [ "$BUILD_BENCHMARKS" == "OFF" ]; then
    BUILD_TESTS=ON
    BUILD_EXAMPLES=ON
    BUILD_BENCHMARKS=ON
fi

BUILD_DIR="$BUILD_PREFIX/$BUILD_INFIX"
# Trigger clean (re-)build
if [ "$CLEAN_BUILD" -eq 1 ]; then
    rm -rf BUILD_DIR
fi
mkdir -p $BUILD_DIR

# The most recent build will be symlinked to cuCollections/build/latest
rm -f $BUILD_PREFIX/latest
ln -sf $BUILD_DIR $BUILD_PREFIX/latest

# Now that BUILD_DIR exists, use readlink to canonicalize the path:
BUILD_DIR=$(readlink -f "${BUILD_DIR}")

BUILD_TYPE=$( [ "$DEBUG_BUILD" -eq 1 ] && echo "Debug" || echo "Release" )

CMAKE_OPTIONS="
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_CXX_STANDARD=${CXX_STANDARD} \
    -DCMAKE_CUDA_STANDARD=${CXX_STANDARD} \
    -DCMAKE_CXX_COMPILER=${HOST_COMPILER} \
    -DCMAKE_CUDA_COMPILER=${CUDA_COMPILER} \
    -DCMAKE_CUDA_HOST_COMPILER=${HOST_COMPILER} \
    -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHS} \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DBUILD_TESTS=${BUILD_TESTS} \
    -DBUILD_EXAMPLES=${BUILD_EXAMPLES} \
    -DBUILD_BENCHMARKS=${BUILD_BENCHMARKS} \
    ${CMAKE_ARGS[*]}
"

echo "[INFO]=============================================="
echo "-- TIMESTAMP: $(date -u +"%Y-%m-%d %H:%M:%S UTC")"
echo "-- GIT_SHA: $(git rev-parse HEAD 2>/dev/null || echo 'N/A')"
echo "-- SRC_DIR: ${REPO_ROOT}"
echo "-- BUILD_DIR: ${BUILD_DIR}"
echo "-- BUILD_TYPE: ${BUILD_TYPE}"
echo "-- PARALLEL_LEVEL: ${PARALLEL_LEVEL}"
echo "-- CUDA_ARCHS: ${CUDA_ARCHS}"
echo "-- BUILD_TESTS: ${BUILD_TESTS}"
echo "-- BUILD_EXAMPLES: ${BUILD_EXAMPLES}"
echo "-- BUILD_BENCHMARKS: ${BUILD_BENCHMARKS}"

if [ ${#CMAKE_ARGS[@]} -gt 0 ]; then
    echo "-- CMAKE_ARGS: ${CMAKE_ARGS[*]}"
else
    echo "-- CMAKE_ARGS: (none)"
fi

# configure
echo "[CONFIGURE]========================================"
cmake -S .. -B $BUILD_DIR $CMAKE_OPTIONS

if command -v sccache >/dev/null; then
    source "./sccache_stats.sh" start
else
    echo "sccache stats: N/A"
fi

#build
echo "[BUILD]============================================"
cmake --build $BUILD_DIR --parallel $PARALLEL_LEVEL
echo "Build complete"

if command -v sccache >/dev/null; then
    source "./sccache_stats.sh" end
else
    echo "sccache stats: N/A"
fi