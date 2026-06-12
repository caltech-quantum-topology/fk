#pragma once

#include <vector>

// Dot product with an affine constant term: a[0] + sum_i a[i+1] * b[i].
// Used to evaluate variable assignments at a lattice point.
int computeDotProduct(const std::vector<int> &a, const std::vector<int> &b);
