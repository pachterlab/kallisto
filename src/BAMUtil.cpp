#include <assert.h>
#include <vector>
#include <iostream>
#include <cstring>
#include "BAMUtil.hpp"

std::string seq_enc = "=ACMGRSVTWYHKDBN";

void parseBAMHeader(gzFile &f, BAMHeader &h) {
  int br, be;

  be = sizeof(h.magic) + sizeof(h.l_text);
  br = gzread(f, &h, be);
  assert(br == be);

  h.text = std::string(h.l_text, '\x00');
  be = h.l_text;
  br = gzread(f, (void *) h.text.data(), be);
  assert(br == be);

  be = sizeof(h.n_ref);
  br = gzread(f, &(h.n_ref), be);
  assert(br == be);
}

void parseBAMRefInfo(gzFile &f, std::vector<BAMRefInfo> &v, int32_t n_ref) {
  int br, be;
  v.reserve(n_ref);
  for (int i = 0; i < n_ref; ++i) {
    v.emplace_back();

    be = sizeof(v[i].l_name);
    br = gzread(f, &v[i], be);
    assert(br == be);

    v[i].name = std::string(v[i].l_name, '\x00');
    be = v[i].l_name;
    br = gzread(f, (void *) v[i].name.data(), be);
    assert(br == be);

    be = sizeof(v[i].l_ref);
    br = gzread(f, &(v[i].l_ref), be);
    assert (br == be);
  }
}

inline void parseSeq(uint8_t *buf, BAMSequence &v, int32_t len, bool odd) {
  if (odd) --len;
  int j = 0;
  for (int i = 0; i < len; ++i, ++buf) {
    v.seq[j++] = seq_enc[*buf >> 4];
    v.seq[j++] = seq_enc[*buf & 0x0F];
  }
  if (odd) {
    v.seq[j++] = seq_enc[*buf >> 4];
  }
  v.l_seq = j;
}

int getBAMSequence(gzFile &f, BAMSequence &s, char *buf) {
  int br, be;
  z_off_t end = gztell(f);

  BAMAlignment aln;
  int alnSize = sizeof(aln);
  br = gzread(f, &aln, alnSize);
  if (br == 0) {
    return 0;
  }
  assert(br == alnSize);
  end += aln.block_size + 4;

  // Skip read name and CIGAR string
  gzseek(f, aln.l_read_name + 4 * aln.n_cigar_op, SEEK_CUR);

  /* Get sequence. */
  assert(aln.l_seq <= BAMBUFSIZE);
  be = (aln.l_seq + 1) / 2;
  br = gzread(f, buf, be);
  assert(br == be);
  // Parse sequence (see BAM specification)
  bool odd = aln.l_seq % 2;
  if (odd) --be;
  s.l_seq = 0;
  for (int i = 0; i < be; ++i, ++buf) {
    s.seq[s.l_seq++] = seq_enc[*(uint8_t *) buf >> 4];
    s.seq[s.l_seq++] = seq_enc[*(uint8_t *) buf & 0x0F];
  }
  if (odd) {
    s.seq[s.l_seq++] = seq_enc[*(uint8_t *) buf >> 4];
  }

  // Skip qual string
  gzseek(f, aln.l_seq, SEEK_CUR);

  /* Parse tags. */
  BAMAuxData aux;
  be = sizeof(aux);
  bool bcumiRead = false;
  while (gztell(f) < end && !bcumiRead) {
    br = gzread(f, &aux, be);
    assert(br == be);

    int l_value;
    switch (aux.val_type) {
      case 'A':
      case 'c':
      case 'C':
        l_value = 1; break;
      case 's':
      case 'S':
        l_value = 2; break;
      case 'i':
      case 'I':
      case 'f':
        l_value = 4; break;
      case 'Z':
        l_value = -1; break;
      case 'H':
        l_value = -2; break;
      case 'B':
        l_value = -3; break;
      default: // Should never happen
        return -1;
    }
    if (l_value == -3) {
      char c;
      br = gzread(f, &c, 1);
      assert(br == 1);
      switch (c) {
        case 'c':
        case 'C':
          l_value = 1; break;
        case 's':
        case 'S':
          l_value = 2; break;
        case 'i':
        case 'I':
        case'f':
          l_value = 4; break;
        default: // Should never happen
          return -2;
      }
      int32_t count;
      br = gzread(f, &count, 4);
      assert(br == 4);

      gzseek(f, l_value * count, SEEK_CUR);
    } else if (l_value < 0) {
      l_value *= -1;
      char elt[2];
      elt[1] = '\xff';
      char *b = NULL;
      int i = 0;
      if (aux.tag[0] == 'C' && aux.tag[1] == 'R') {
        b = s.bc;
      } else if (aux.tag[0] == 'U' && aux.tag[1] == 'R') {
        b = s.umi;
      }
      while (true) {
        br = gzread(f, &elt, l_value);
        assert(br == l_value);
        if (!(elt[0] & elt[1])) {
          break;
        }
        if (b) {
          b[i++] = elt[0];
        }
      }
      assert(i <= 32);
      if (aux.tag[0] == 'C' && aux.tag[1] == 'R') {
        s.l_bc = i;
      } else if (aux.tag[0] == 'U' && aux.tag[1] == 'R') {
        s.l_umi = i;
      }
    } else if (aux.tag[0] == 'N' && aux.tag[1] == 'H') {
      br = gzread(f, &s.nh, l_value);
      assert(br == l_value);
    } else {
      gzseek(f, l_value, SEEK_CUR);
    }
  }

  gzseek(f, end, SEEK_SET);

  return s.l_seq;
}

