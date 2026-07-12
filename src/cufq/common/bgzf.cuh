#pragma once

// Shared BGZF (Block GZIP Format) structures and utilities
// Used by both the standalone cufq tools and the gpukallisto GPU pipeline

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// BGZF block header structure (RFC 1952 + SAM spec)
#pragma pack(push, 1)
struct BgzfHeader {
    uint8_t id1;        // 31
    uint8_t id2;        // 139
    uint8_t cm;         // 8
    uint8_t flg;        // 4
    uint32_t mtime;
    uint8_t xfl;
    uint8_t os;
    uint16_t xlen;      // 6
    uint8_t si1;        // 66 'B'
    uint8_t si2;        // 67 'C'
    uint16_t slen;      // 2
    uint16_t bsize;     // Total block size - 1
};
#pragma pack(pop)

struct BgzfBlock {
    size_t file_offset;
    size_t compressed_size;
    size_t uncompressed_size;
};

// Check if memory buffer is BGZF format
inline bool is_bgzf_memory(const char* data, size_t size) {
    if (size < sizeof(BgzfHeader)) return false;
    
    const BgzfHeader* header = reinterpret_cast<const BgzfHeader*>(data);
    return header->id1 == 31 && header->id2 == 139 &&
           header->cm == 8 && (header->flg & 4) &&
           header->si1 == 66 && header->si2 == 67;
}

// Check BGZF format from file header only (no full-file read)
inline bool is_bgzf_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char buf[sizeof(BgzfHeader)];
    f.read(buf, sizeof(buf));
    if (static_cast<size_t>(f.gcount()) < sizeof(BgzfHeader)) return false;
    return is_bgzf_memory(buf, sizeof(BgzfHeader));
}

// Parse BGZF blocks from memory buffer (no file seeks!)
// Scans sequentially through the buffer to find all blocks
inline std::vector<BgzfBlock> parse_bgzf_blocks_from_memory(const char* data, size_t file_size) {
    std::vector<BgzfBlock> blocks;
    blocks.reserve(file_size / (32 * 1024)); // Estimate ~32KB average block
    
    size_t offset = 0;
    size_t invalid_blocks = 0;
    
    while (offset + 28 < file_size) { // Need at least header + trailer
        const BgzfHeader* header = reinterpret_cast<const BgzfHeader*>(data + offset);
        
        // Validate BGZF magic
        if (header->id1 != 31 || header->id2 != 139) break;
        
        // Validate extra field for BGZF
        if (header->si1 != 66 || header->si2 != 67) break;
        
        size_t block_size = header->bsize + 1;
        if (offset + block_size > file_size) break;
        if (block_size < 28) break; // Minimum valid BGZF block size
        
        // Read isize from end of block (last 4 bytes)
        uint32_t isize;
        memcpy(&isize, data + offset + block_size - 4, sizeof(uint32_t));
        
        // Validate ISIZE (BGZF max uncompressed size is 65536)
        if (isize > 65536) {
            if (invalid_blocks == 0) {
                fprintf(stderr, "Warning: Block at offset %zu has invalid ISIZE=%u (max 65536)\n",
                        offset, isize);
            }
            invalid_blocks++;
            isize = 65536;
        }
        
        BgzfBlock block;
        block.file_offset = offset;
        block.compressed_size = block_size;
        block.uncompressed_size = isize;
        blocks.push_back(block);
        
        offset += block_size;
    }
    
    if (invalid_blocks > 0) {
        fprintf(stderr, "Warning: %zu blocks had invalid ISIZE values\n", invalid_blocks);
    }
    
    return blocks;
}
