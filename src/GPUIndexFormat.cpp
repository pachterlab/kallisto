#include "GPUIndexFormat.h"

#include <fstream>
#include <algorithm>
#include <cstring>
#include <iostream>

#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#define HAVE_MMAP 1
#endif

#include "CompactedDBG.hpp"
#include "Node.hpp"

namespace {

uint8_t encode_base(char c) {
  c = ::toupper(static_cast<unsigned char>(c));
  if (c == 'A') return 0;
  if (c == 'C') return 1;
  if (c == 'G') return 2;
  if (c == 'T') return 3;
  return 0;  // Invalid: default to A
}

char decode_base(uint8_t b) {
  const char* table = "ACGT";
  return table[b & 3];
}

void pack_sequence(const std::string& seq, std::vector<uint8_t>& out) {
  size_t packed_size = (seq.size() + 3) / 4;
  out.resize(packed_size, 0);
  for (size_t i = 0; i < seq.size(); ++i) {
    uint8_t code = encode_base(seq[i]);
    size_t byte_idx = i / 4;
    size_t bit_idx = (3 - (i % 4)) * 2;
    out[byte_idx] |= (code << bit_idx);
  }
}

void unpack_sequence(const uint8_t* packed, size_t seq_len, std::string& out) {
  out.resize(seq_len);
  for (size_t i = 0; i < seq_len; ++i) {
    size_t byte_idx = i / 4;
    size_t bit_idx = (3 - (i % 4)) * 2;
    uint8_t code = (packed[byte_idx] >> bit_idx) & 3;
    out[i] = decode_base(code);
  }
}

}  // namespace

bool GPUIndex::load(const std::string& path) {
#ifdef HAVE_MMAP
  int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    std::cerr << "Error: could not open GPU index file: " << path << std::endl;
    return false;
  }
  struct stat st;
  if (fstat(fd, &st) != 0) {
    close(fd);
    return false;
  }
  size_t file_size = static_cast<size_t>(st.st_size);
  void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mapped == MAP_FAILED) {
    std::cerr << "Error: mmap failed for " << path << std::endl;
    return false;
  }
  const uint8_t* p = static_cast<const uint8_t*>(mapped);
  const uint8_t* end = p + file_size;

  auto read = [&](void* dst, size_t n) {
    if (p + n > end) return false;
    memcpy(dst, p, n);
    p += n;
    return true;
  };

  uint32_t magic = 0;
  if (!read(&magic, sizeof(magic)) || magic != GPU_INDEX_MAGIC) {
    munmap(mapped, file_size);
    std::cerr << "Error: invalid GPU index magic (not a .gpuidx file)" << std::endl;
    return false;
  }

  uint32_t version = 0;
  if (!read(&version, sizeof(version)) || version != GPU_INDEX_FORMAT_VERSION) {
    munmap(mapped, file_size);
    std::cerr << "Error: incompatible GPU index version " << version
              << ", expected " << GPU_INDEX_FORMAT_VERSION << std::endl;
    return false;
  }

  uint32_t nt = 0, ne = 0, nc = 0;
  if (!read(&k, sizeof(k)) || !read(&nt, sizeof(nt)) || !read(&ne, sizeof(ne)) ||
      !read(&nc, sizeof(nc))) {
    munmap(mapped, file_size);
    return false;
  }
  num_transcripts = nt;
  num_ecs = ne;
  num_contigs = nc;

  ecmapinv.clear();
  ecmapinv.reserve(ne);
  for (uint32_t ec_id = 0; ec_id < ne; ++ec_id) {
    uint32_t count = 0;
    if (!read(&count, sizeof(count))) {
      munmap(mapped, file_size);
      return false;
    }
    if (count > 0 && p + count * sizeof(uint32_t) > end) {
      munmap(mapped, file_size);
      return false;
    }
    Roaring r;
    const uint32_t* ids = reinterpret_cast<const uint32_t*>(p);
    for (uint32_t i = 0; i < count; ++i) {
      r.add(ids[i]);
    }
    p += count * sizeof(uint32_t);
    r.runOptimize();
    ecmapinv.insert({std::move(r), static_cast<int32_t>(ec_id)});
  }

  target_lens_.resize(num_transcripts);
  if (!read(target_lens_.data(), num_transcripts * sizeof(uint32_t))) {
    munmap(mapped, file_size);
    return false;
  }

  target_names_.resize(num_transcripts);
  for (size_t i = 0; i < num_transcripts; ++i) {
    uint32_t name_len = 0;
    if (!read(&name_len, sizeof(name_len))) {
      munmap(mapped, file_size);
      return false;
    }
    target_names_[i].resize(name_len);
    if (name_len > 0 && !read(&target_names_[i][0], name_len)) {
      munmap(mapped, file_size);
      return false;
    }
  }

  contig_sequences.clear();
  contig_ec_blocks.clear();
  contig_sequences.reserve(num_contigs);
  contig_ec_blocks.reserve(num_contigs);

  for (size_t c = 0; c < num_contigs; ++c) {
    uint32_t seq_len = 0;
    if (!read(&seq_len, sizeof(seq_len))) {
      munmap(mapped, file_size);
      return false;
    }
    size_t packed_size = (seq_len + 3) / 4;
    if (p + packed_size > end) {
      munmap(mapped, file_size);
      return false;
    }
    std::string seq;
    unpack_sequence(p, seq_len, seq);
    p += packed_size;
    contig_sequences.push_back(std::move(seq));

    uint32_t num_blocks = 0;
    if (!read(&num_blocks, sizeof(num_blocks))) {
      munmap(mapped, file_size);
      return false;
    }
    std::vector<ECBlock> blocks(num_blocks);
    for (uint32_t b = 0; b < num_blocks; ++b) {
      if (!read(&blocks[b].start, sizeof(blocks[b].start)) ||
          !read(&blocks[b].end, sizeof(blocks[b].end)) ||
          !read(&blocks[b].ec_id, sizeof(blocks[b].ec_id))) {
        munmap(mapped, file_size);
        return false;
      }
    }
    contig_ec_blocks.push_back(std::move(blocks));
  }

  munmap(mapped, file_size);
#else
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) {
    std::cerr << "Error: could not open GPU index file: " << path << std::endl;
    return false;
  }

  uint32_t magic = 0;
  in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
  if (magic != GPU_INDEX_MAGIC) {
    std::cerr << "Error: invalid GPU index magic (not a .gpuidx file)" << std::endl;
    return false;
  }

  uint32_t version = 0;
  in.read(reinterpret_cast<char*>(&version), sizeof(version));
  if (version != GPU_INDEX_FORMAT_VERSION) {
    std::cerr << "Error: incompatible GPU index version " << version
              << ", expected " << GPU_INDEX_FORMAT_VERSION << std::endl;
    return false;
  }

  in.read(reinterpret_cast<char*>(&k), sizeof(k));
  uint32_t nt = 0, ne = 0, nc = 0;
  in.read(reinterpret_cast<char*>(&nt), sizeof(nt));
  in.read(reinterpret_cast<char*>(&ne), sizeof(ne));
  in.read(reinterpret_cast<char*>(&nc), sizeof(nc));
  num_transcripts = nt;
  num_ecs = ne;
  num_contigs = nc;

  ecmapinv.clear();
  ecmapinv.reserve(ne);
  for (uint32_t ec_id = 0; ec_id < ne; ++ec_id) {
    uint32_t count = 0;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));
    std::vector<uint32_t> ids(count);
    if (count > 0) {
      in.read(reinterpret_cast<char*>(ids.data()), count * sizeof(uint32_t));
    }
    Roaring r;
    for (uint32_t id : ids) {
      r.add(id);
    }
    r.runOptimize();
    ecmapinv.insert({std::move(r), static_cast<int32_t>(ec_id)});
  }

  target_lens_.resize(num_transcripts);
  in.read(reinterpret_cast<char*>(target_lens_.data()),
          num_transcripts * sizeof(uint32_t));

  target_names_.resize(num_transcripts);
  for (size_t i = 0; i < num_transcripts; ++i) {
    uint32_t name_len = 0;
    in.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
    std::string name(name_len, '\0');
    if (name_len > 0) {
      in.read(&name[0], name_len);
    }
    target_names_[i] = std::move(name);
  }

  contig_sequences.clear();
  contig_ec_blocks.clear();
  contig_sequences.reserve(num_contigs);
  contig_ec_blocks.reserve(num_contigs);

  for (size_t c = 0; c < num_contigs; ++c) {
    uint32_t seq_len = 0;
    in.read(reinterpret_cast<char*>(&seq_len), sizeof(seq_len));
    size_t packed_size = (seq_len + 3) / 4;
    std::vector<uint8_t> packed(packed_size);
    if (packed_size > 0) {
      in.read(reinterpret_cast<char*>(packed.data()), packed_size);
    }
    std::string seq;
    unpack_sequence(packed.data(), seq_len, seq);
    contig_sequences.push_back(std::move(seq));

    uint32_t num_blocks = 0;
    in.read(reinterpret_cast<char*>(&num_blocks), sizeof(num_blocks));
    std::vector<ECBlock> blocks(num_blocks);
    for (uint32_t b = 0; b < num_blocks; ++b) {
      in.read(reinterpret_cast<char*>(&blocks[b].start), sizeof(blocks[b].start));
      in.read(reinterpret_cast<char*>(&blocks[b].end), sizeof(blocks[b].end));
      in.read(reinterpret_cast<char*>(&blocks[b].ec_id), sizeof(blocks[b].ec_id));
    }
    contig_ec_blocks.push_back(std::move(blocks));
  }
#endif

  std::cerr << "[index] loaded GPU index: " << num_contigs << " contigs, "
            << num_ecs << " ECs, " << num_transcripts << " transcripts" << std::endl;
  return true;
}

bool GPUIndex::write(const std::string& path) const {
  std::ofstream out(path, std::ios::binary);
  if (!out.is_open()) {
    std::cerr << "Error: could not open output file: " << path << std::endl;
    return false;
  }

  out.write(reinterpret_cast<const char*>(&GPU_INDEX_MAGIC), sizeof(GPU_INDEX_MAGIC));
  out.write(reinterpret_cast<const char*>(&GPU_INDEX_FORMAT_VERSION),
            sizeof(GPU_INDEX_FORMAT_VERSION));
  out.write(reinterpret_cast<const char*>(&k), sizeof(k));
  uint32_t nt = static_cast<uint32_t>(num_transcripts);
  uint32_t ne = static_cast<uint32_t>(num_ecs);
  uint32_t nc = static_cast<uint32_t>(num_contigs);
  out.write(reinterpret_cast<const char*>(&nt), sizeof(nt));
  out.write(reinterpret_cast<const char*>(&ne), sizeof(ne));
  out.write(reinterpret_cast<const char*>(&nc), sizeof(nc));

  // Write ecmapinv: need to serialize in EC-id order
  // ecmapinv maps Roaring -> ec_id, we need to output by ec_id
  std::vector<Roaring> by_ec(num_ecs);
  for (const auto& p : ecmapinv) {
    int32_t ec_id = p.second;
    if (ec_id >= 0 && static_cast<size_t>(ec_id) < num_ecs) {
      by_ec[ec_id] = p.first;
    }
  }
  for (size_t ec_id = 0; ec_id < num_ecs; ++ec_id) {
    const Roaring& r = by_ec[ec_id];
    uint32_t count = static_cast<uint32_t>(r.cardinality());
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    for (uint32_t tr : r) {
      out.write(reinterpret_cast<const char*>(&tr), sizeof(tr));
    }
  }

  // Write target metadata
  out.write(reinterpret_cast<const char*>(target_lens_.data()),
            num_transcripts * sizeof(uint32_t));
  for (const auto& name : target_names_) {
    uint32_t name_len = static_cast<uint32_t>(name.size());
    out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
    if (name_len > 0) {
      out.write(name.data(), name_len);
    }
  }

  // Write contigs
  for (size_t c = 0; c < num_contigs; ++c) {
    const std::string& seq = contig_sequences[c];
    uint32_t seq_len = static_cast<uint32_t>(seq.size());
    out.write(reinterpret_cast<const char*>(&seq_len), sizeof(seq_len));
    std::vector<uint8_t> packed;
    pack_sequence(seq, packed);
    out.write(reinterpret_cast<const char*>(packed.data()), packed.size());

    const auto& blocks = contig_ec_blocks[c];
    uint32_t num_blocks = static_cast<uint32_t>(blocks.size());
    out.write(reinterpret_cast<const char*>(&num_blocks), sizeof(num_blocks));
    for (const auto& blk : blocks) {
      out.write(reinterpret_cast<const char*>(&blk.start), sizeof(blk.start));
      out.write(reinterpret_cast<const char*>(&blk.end), sizeof(blk.end));
      out.write(reinterpret_cast<const char*>(&blk.ec_id), sizeof(blk.ec_id));
    }
  }

  out.close();
  std::cerr << "[index] wrote GPU index: " << path << std::endl;
  return true;
}

bool convert_kallisto_to_gpu_index(const std::string& kallisto_index_path,
                                  const std::string& gpuidx_path,
                                  ProgramOptions& opt) {
  opt.index = kallisto_index_path;
  KmerIndex index(opt);
  index.load(opt, true, true, /*gpuMode=*/true);

  GPUIndex gpu_idx;
  gpu_idx.k = index.k;
  gpu_idx.num_transcripts = index.num_trans;
  gpu_idx.num_ecs = index.ecmapinv.size();
  gpu_idx.num_contigs = index.dbg.size();
  gpu_idx.target_names_ = index.target_names_;
  gpu_idx.target_lens_ = index.target_lens_;
  gpu_idx.ecmapinv = index.ecmapinv;

  gpu_idx.contig_sequences.reserve(gpu_idx.num_contigs);
  gpu_idx.contig_ec_blocks.reserve(gpu_idx.num_contigs);

  std::vector<SparseVector<uint32_t>> vals;
  for (const auto& contig : index.dbg) {
    auto n = contig.getData();
    std::string seq = contig.referenceUnitigToString();
    gpu_idx.contig_sequences.push_back(seq);

    n->ec.get_vals(vals);
    std::vector<ECBlock> blocks;
    int j = 0;
    size_t contigpos = 0;
    while (contigpos < contig.len) {
      auto mc = n->ec.get_block_at(contigpos);
      const auto& val = vals[j];
      const Roaring& trs = val.getIndices();
      auto ec_it = index.ecmapinv.find(trs);
      if (ec_it != index.ecmapinv.end()) {
        ECBlock blk;
        blk.start = static_cast<uint32_t>(mc.first);
        blk.end = static_cast<uint32_t>(mc.second);
        blk.ec_id = ec_it->second;
        blocks.push_back(blk);
      }
      contigpos = mc.second;
      ++j;
    }
    gpu_idx.contig_ec_blocks.push_back(std::move(blocks));
  }

  return gpu_idx.write(gpuidx_path);
}

bool is_gpu_index(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) {
    return false;
  }
  uint32_t magic = 0;
  in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
  return (magic == GPU_INDEX_MAGIC);
}
