#include "common.h"

std::string pretty_num(int num) {
  return pretty_num(static_cast<size_t>(num));
}

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");

  if (s.size() <= 3) {
    return s;
  }

  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }

  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }

  ret += s.substr(start_pos, 3);

  return ret;
}
