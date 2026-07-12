#include "RapidGzipReader.h"

#include <cstdio>
#include <stdexcept>

#ifdef USE_RAPIDGZIP
#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>
#else
#include <zlib.h>
#endif

#ifdef USE_RAPIDGZIP

struct RapidGzipReader::Impl {
    rapidgzip::ParallelGzipReader<> reader;

    Impl(const std::string& path, size_t num_threads, uint64_t chunk_size)
        : reader(std::make_unique<rapidgzip::StandardFileReader>(path),
                 num_threads,
                 chunk_size)
    {
        reader.setCRC32Enabled(false);
    }
};

RapidGzipReader::RapidGzipReader(const std::string& path,
                                 size_t num_threads,
                                 uint64_t chunk_size)
    : impl_(std::make_unique<Impl>(path, num_threads, chunk_size))
{
}

RapidGzipReader::~RapidGzipReader() = default;

ssize_t RapidGzipReader::read(char* buf, size_t size)
{
    try {
        const size_t n = impl_->reader.read(buf, size);
        return static_cast<ssize_t>(n);
    } catch (const std::exception& e) {
        fprintf(stderr, "RapidGzipReader error: %s\n", e.what());
        return -1;
    }
}

bool RapidGzipReader::eof() const
{
    return impl_->reader.eof();
}

#else

struct RapidGzipReader::Impl {
    gzFile file;

    explicit Impl(const std::string& path)
        : file(gzopen(path.c_str(), "rb"))
    {
        if (!file) {
            throw std::runtime_error("Cannot open gzip file: " + path);
        }
        gzbuffer(file, 256 * 1024);
    }

    ~Impl()
    {
        if (file) {
            gzclose(file);
        }
    }
};

RapidGzipReader::RapidGzipReader(const std::string& path,
                                 size_t /*num_threads*/,
                                 uint64_t /*chunk_size*/)
    : impl_(std::make_unique<Impl>(path))
{
}

RapidGzipReader::~RapidGzipReader() = default;

ssize_t RapidGzipReader::read(char* buf, size_t size)
{
    const int n = gzread(impl_->file, buf, static_cast<unsigned>(size));
    if (n < 0) {
        int err = 0;
        const char* msg = gzerror(impl_->file, &err);
        fprintf(stderr, "Gzip error: %s\n", msg);
        return -1;
    }
    return n;
}

bool RapidGzipReader::eof() const
{
    return gzeof(impl_->file) != 0;
}

#endif
