#ifndef GPU_INDEX_FORMAT_H
#define GPU_INDEX_FORMAT_H

#include "common.h"
#include "KmerIndex.h"

#include <string>
#include <vector>
#include <cstdint>

// GPU index format version
constexpr uint32_t GPU_INDEX_FORMAT_VERSION = 1;
constexpr uint32_t GPU_INDEX_MAGIC = 0x47494458;  // "GIDX" in little-endian

struct ECBlock {
  uint32_t start;  // 0-based k-mer position
  uint32_t end;    // exclusive end (k-mer position)
  int32_t ec_id;
};

struct GPUIndex {
  int k;
  size_t num_transcripts;
  size_t num_ecs;
  size_t num_contigs;

  std::vector<std::string> target_names_;
  std::vector<uint32_t> target_lens_;

  // EC map inverse: for each EC id, the set of transcript IDs (as Roaring)
  // Built from serialized format (count + ids per EC)
  EcMapInv ecmapinv;

  // Contig data: sequences and EC blocks per contig
  std::vector<std::string> contig_sequences;
  std::vector<std::vector<ECBlock>> contig_ec_blocks;

  bool load(const std::string& path);
  bool write(const std::string& path) const;
};

// Convert kallisto index to GPU index format
// Loads the full kallisto index and writes .gpuidx
bool convert_kallisto_to_gpu_index(const std::string& kallisto_index_path,
                                   const std::string& gpuidx_path,
                                   ProgramOptions& opt);

// Check if file is a GPU index (by magic)
bool is_gpu_index(const std::string& path);

#endif  // GPU_INDEX_FORMAT_H
