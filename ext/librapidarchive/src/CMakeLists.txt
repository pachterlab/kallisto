# Vendored for gpukallisto (offline builds). Do not fetch via git submodules here.
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts/include/cxxopts.hpp OR
   NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/zlib-ng/inflate.h OR
   NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/rpmalloc/rpmalloc/rpmalloc.h)
    message(FATAL_ERROR
        "librapidarchive external dependencies are incomplete. "
        "See ext/librapidarchive/README.gpukallisto for vendoring instructions.")
endif()

add_library(cxxopts INTERFACE)
target_include_directories(cxxopts SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts/include)
target_sources(cxxopts INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts/include/cxxopts.hpp)


# rpmalloc

if(LIBRAPIDARCHIVE_WITH_RPMALLOC)
    project(rpmalloc C)

    add_library(rpmalloc STATIC)
    set(RPMALLOC_HOME "${CMAKE_CURRENT_SOURCE_DIR}/external/rpmalloc/rpmalloc")
    target_include_directories(rpmalloc SYSTEM INTERFACE ${RPMALLOC_HOME})
    target_sources(rpmalloc PRIVATE
        ${RPMALLOC_HOME}/rpmalloc.c
        ${RPMALLOC_HOME}/rpmalloc.h
        # Do not use! Overriding operator new/delete is too bug-prone.
        #${RPMALLOC_HOME}/rpnew.h
    )
    set_target_properties(rpmalloc PROPERTIES C_STANDARD 11)
    set_target_properties(rpmalloc PROPERTIES LINKER_LANGUAGE C)
    # ENABLE_OVERRIDE=1 results in a segfault from _dl_start_user
    target_compile_definitions(rpmalloc PUBLIC ENABLE_OVERRIDE=0 ENABLE_ASSERTS=0)
    if(MSVC)
        # To fix errors with atomic support and the -D_WIN32_WINNT to fix this error:
        # rpmalloc.c(2110): error C2440: 'function': cannot convert from 'void (__cdecl *)(void *)' to 'PFLS_CALLBACK_FUNCTION'
        target_compile_options(rpmalloc PUBLIC /experimental:c11atomics -D_WIN32_WINNT=0xA00)
    endif()

    option(ENABLE_RPMALLOC_STATISTICS "Enables statistics for rpmalloc" OFF)
    if (ENABLE_RPMALLOC_STATISTICS)
        target_compile_definitions(rpmalloc PRIVATE -DENABLE_STATISTICS=1)
    endif()
endif()

#add_subdirectory(external/zlib)

if(NOT LIBRAPIDARCHIVE_USE_SYSTEM_ZLIB AND LIBRAPIDARCHIVE_USE_ZLIB_NG)
    project(zlibstatic C)
    add_library(zlibstatic STATIC)
    set(ZLIB_HOME "${CMAKE_CURRENT_SOURCE_DIR}/external/zlib-ng")
    target_include_directories(zlibstatic SYSTEM INTERFACE ${ZLIB_HOME} ${CMAKE_CURRENT_BINARY_DIR}/zlib-ng)
    target_include_directories(zlibstatic PUBLIC ${ZLIB_HOME} ${CMAKE_CURRENT_BINARY_DIR}/zlib-ng)

    # Not working, see comment for Z_PREFIX.
    #set(ZLIB_SYMBOL_PREFIX LIBRAPIDARCHIVE_)

    configure_file(${ZLIB_HOME}/zconf-ng.h.in ${ZLIB_HOME}/zconf-ng.h @ONLY)
    configure_file(${ZLIB_HOME}/zconf.h.in ${ZLIB_HOME}/zconf.h @ONLY)
    configure_file(${ZLIB_HOME}/zlib.h.in ${ZLIB_HOME}/zlib.h @ONLY)
    #configure_file(${ZLIB_HOME}/zlib_name_mangling.h.in ${ZLIB_HOME}/zlib_name_mangling.h @ONLY)
    configure_file(${ZLIB_HOME}/zlib_name_mangling.h.empty ${ZLIB_HOME}/zlib_name_mangling.h COPYONLY)

    target_sources(zlibstatic PRIVATE
        ${ZLIB_HOME}/cpu_features.c
        ${ZLIB_HOME}/inflate.c
        ${ZLIB_HOME}/inflate.h
        ${ZLIB_HOME}/adler32.c
        ${ZLIB_HOME}/crc32.c
        ${ZLIB_HOME}/crc32_braid_comb.c
        ${ZLIB_HOME}/functable.c
        ${ZLIB_HOME}/inftrees.c
        ${ZLIB_HOME}/inftrees.h
        ${ZLIB_HOME}/gzguts.h
        ${ZLIB_HOME}/zutil.c
        ${ZLIB_HOME}/zutil.h
        # Only required for compression, e.g., used for the tests and benchmarks
        ${ZLIB_HOME}/deflate.c
        ${ZLIB_HOME}/deflate.h
        ${ZLIB_HOME}/deflate_fast.c
        ${ZLIB_HOME}/deflate_huff.c
        ${ZLIB_HOME}/deflate_medium.c
        ${ZLIB_HOME}/deflate_quick.c
        ${ZLIB_HOME}/deflate_rle.c
        ${ZLIB_HOME}/deflate_slow.c
        ${ZLIB_HOME}/deflate_stored.c
        ${ZLIB_HOME}/insert_string.c
        ${ZLIB_HOME}/insert_string_roll.c
        ${ZLIB_HOME}/trees.c
        ${ZLIB_HOME}/trees.h
        # Generic implementations
        ${ZLIB_HOME}/arch/generic/adler32_c.c
        ${ZLIB_HOME}/arch/generic/adler32_fold_c.c
        ${ZLIB_HOME}/arch/generic/chunk_permute_table.h
        ${ZLIB_HOME}/arch/generic/chunkset_c.c
        ${ZLIB_HOME}/arch/generic/compare256_c.c
        ${ZLIB_HOME}/arch/generic/compare256_p.h
        ${ZLIB_HOME}/arch/generic/crc32_braid_c.c
        ${ZLIB_HOME}/arch/generic/crc32_fold_c.c
        ${ZLIB_HOME}/arch/generic/generic_functions.h
        ${ZLIB_HOME}/arch/generic/slide_hash_c.c
    )
    set_target_properties(zlibstatic PROPERTIES LINKER_LANGUAGE C)
    set_target_properties(zlibstatic PROPERTIES C_STANDARD 11)
    target_compile_definitions(zlibstatic PUBLIC ZLIB_COMPAT)
elseif(NOT LIBRAPIDARCHIVE_USE_SYSTEM_ZLIB)
    project(zlibstatic C)
    add_library(zlibstatic STATIC)
    set(ZLIB_HOME "${CMAKE_CURRENT_SOURCE_DIR}/external/zlib")
    target_include_directories(zlibstatic SYSTEM INTERFACE ${ZLIB_HOME})
    target_sources(zlibstatic PRIVATE
        ${ZLIB_HOME}/inflate.c
        ${ZLIB_HOME}/inflate.h
        ${ZLIB_HOME}/adler32.c  # zutil.h is the corresponding header
        ${ZLIB_HOME}/crc32.c
        ${ZLIB_HOME}/crc32.h
        ${ZLIB_HOME}/inffixed.h
        ${ZLIB_HOME}/inffast.c
        ${ZLIB_HOME}/inffast.h
        ${ZLIB_HOME}/inftrees.c
        ${ZLIB_HOME}/inftrees.h
        ${ZLIB_HOME}/gzguts.h
        ${ZLIB_HOME}/zutil.c
        ${ZLIB_HOME}/zutil.h
        ${ZLIB_HOME}/zlib.h
        # Only required for compression, e.g., used for the tests and benchmarks
        ${ZLIB_HOME}/deflate.c
        ${ZLIB_HOME}/deflate.h
        ${ZLIB_HOME}/trees.c
        ${ZLIB_HOME}/trees.h
    )
    set_target_properties(zlibstatic PROPERTIES LINKER_LANGUAGE C)
    # It might be better to add a prefix to avoid symbol name collision when this is built statically
    # into a shared library! However, if zlib.h is not the first include, then it leads to compilation errors
    # because member calls named crc32 are also searched-and-replaced by the preprocessor but not the definition.
    # src/gzip/zlib.hpp: In member function ‘rapidgzip::Footer rapidgzip::ZlibInflateWrapper::readGzipFooter()’:
    # src/gzip/zlib.hpp:470:16: error: ‘struct rapidgzip::gzip::Footer’ has no member named ‘LIBRAPIDARCHIVE_crc32’
    #   470 |         footer.crc32 += static_cast<uint32_t>( subbyte ) << ( i * BYTE_SIZE );
    # Furthermore, it accidentally changes the symbols of the public librapidarchive API! We can't have that.
    #target_compile_definitions(zlibstatic PUBLIC Z_PREFIX=LIBRAPIDARCHIVE_)
endif()


# ISA-l

if(LIBRAPIDARCHIVE_WITH_ISAL)
    include(FetchContent)

    FetchContent_Declare(isal_project SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/isa-l)
    # https://discourse.cmake.org/t/fetchcontent-declare-with-cmake-args-does-not-work/11227
    # https://discourse.cmake.org/t/what-is-the-correct-way-to-set-options-of-a-project-before-fetch-content/268
    set(ISAL_BUILD_TESTS OFF CACHE BOOL "Disable ISA-L tests" FORCE)
    set(ISAL_BUILD_SHARED_LIBS OFF CACHE BOOL "Build static ISA-L" FORCE)
    FetchContent_MakeAvailable(isal_project)

    # I am not using the guarded <isa-l/igzip_lib.h> includes because it would not work in the current setup.py :(.
    # The guarded includes are also not necessary or possible with FetchContent, only with ExternalProject_Add!
    add_library(ISAL::isal ALIAS isal)
endif()


add_subdirectory(core)
add_subdirectory(filereader)
add_subdirectory(huffman)
add_subdirectory(indexed_bzip2)
add_subdirectory(rapidgzip)
add_subdirectory(tools)
if (LIBRAPIDARCHIVE_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
if (LIBRAPIDARCHIVE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
