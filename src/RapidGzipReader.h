#ifndef RAPIDGZIP_READER_H
#define RAPIDGZIP_READER_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

// Thin wrapper around librapidarchive ParallelGzipReader (or zlib fallback).
// Keeps rapidgzip headers out of CUDA translation units.
class RapidGzipReader {
public:
    explicit RapidGzipReader(const std::string& path,
                             size_t num_threads = 0,
                             uint64_t chunk_size = 4ULL * 1024 * 1024);

    ~RapidGzipReader();

    RapidGzipReader(const RapidGzipReader&) = delete;
    RapidGzipReader& operator=(const RapidGzipReader&) = delete;

    // Returns bytes read, 0 at EOF, -1 on error.
    ssize_t read(char* buf, size_t size);

    bool eof() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

#endif
