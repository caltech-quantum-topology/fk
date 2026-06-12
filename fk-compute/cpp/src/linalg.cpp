#include "fk/linalg.hpp"

int computeDotProduct(const std::vector<int> &a, const std::vector<int> &b) {
  int accumulator = a[0];
  for (size_t z = 0; z < b.size(); z++) {
    accumulator += a[z + 1] * b[z];
  }
  return accumulator;
}
