#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include "fk/bilvector.hpp"


// Multiply by q^power: shift exponents
template <typename T>
bilvector<T> multiplyByQPower(const bilvector<T> &poly, int power) {
  if (power == 0) {
    return poly;
  }

  int inMin = poly.getMaxNegativeIndex();
  int inMax = poly.getMaxPositiveIndex();
  int outMin = inMin + power;
  int outMax = inMax + power;

  int componentSize = poly.getComponentSize();

  int negativeCount = 0;
  if (outMin < 0) {
    negativeCount = -outMin;
  }

  int positiveCount = 0;
  if (outMax >= 0) {
    positiveCount = outMax + 1;
  }

  bilvector<T> result(negativeCount, positiveCount, componentSize, T{});

  for (int e = inMin; e <= inMax; ++e) {
    T c = poly[e];
    if (c == T{}) continue;
    result[e + power] += c;
  }

  return result;
}

// ============ Non-recursive positive q-binomial as bilvector<int> ===========

bilvector<int> computePositiveQBinomialHelper(int upperLimit, int lowerLimit) {
  int n = upperLimit;
  int k = lowerLimit;

  // Zero polynomial if out of range
  if (k < 0 || n < 0 || k > n) {
    return bilvector<int>(0, 1, 1, 0);  // all zero, only exponent 0 allocated
  }

  auto makeConstOne = []() {
    bilvector<int> p(0, 1, 1, 0);  // exponent 0 only
    p[0] = 1;
    return p;
  };

  // C[i][j] = [i choose j]_q
  std::vector<std::vector<bilvector<int>>> C(
      n + 1, std::vector<bilvector<int>>(k + 1, bilvector<int>(0, 1, 1, 0)));

  C[0][0] = makeConstOne();

  for (int i = 1; i <= n; ++i) {
    C[i][0] = makeConstOne();  // [i choose 0]_q = 1

    int jMax = std::min(i, k);
    for (int j = 1; j <= jMax; ++j) {
      if (j == i) {
        C[i][j] = makeConstOne();  // [i choose i]_q = 1
      } else {
        // [i choose j]_q = [i-1 choose j]_q + q^(i-j) [i-1 choose j-1]_q
        bilvector<int> shifted = multiplyByQPower(C[i - 1][j - 1], i - j);
        C[i][j] = C[i - 1][j] + shifted;
      }
    }
  }

  return C[n][k];
}

// ===================== Pretty-print q-polynomial ======================

void printQBinomial(int n, int k, const bilvector<int> &qb) {
  int minExp = qb.getMaxNegativeIndex();
  int maxExp = qb.getMaxPositiveIndex();

  std::cout << "[" << n << " choose " << k << "]_q = ";

  bool first = true;
  for (int e = minExp; e <= maxExp; ++e) {
    int c = qb[e];
    if (c == 0) continue;
    if (!first) {
      std::cout << " + ";
    }
    first = false;
    std::cout << c;
    if (e != 0) {
      std::cout << "*q^" << e;
    }
  }
  if (first) {
    std::cout << "0";
  }
  std::cout << "\n";
}

// =============================== main =================================

int main() {
  // Some sample q-binomials
  std::vector<std::pair<int,int>> tests = {
      {3, 1},
      {4, 2},
      {5, 2},
      {5, 3}
  };

  for (auto [n, k] : tests) {
    bilvector<int> qb = computePositiveQBinomialHelper(n, k);
    printQBinomial(n, k, qb);

    // Also dump coefficients explicitly
    int minExp = qb.getMaxNegativeIndex();
    int maxExp = qb.getMaxPositiveIndex();
    std::cout << "  Coefficients (q^e):\n";
    for (int e = minExp; e <= maxExp; ++e) {
      std::cout << "    e = " << e << " : " << qb[e] << "\n";
    }
    std::cout << "\n";
  }

  return 0;
}
