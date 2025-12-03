#include "fk/qalg_links.hpp"
#include "fk/linalg.hpp"
#include <map>

void computePositiveQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift) {
  if (upperLimit == lowerLimit) {
    binomialCoefficients[shift] += 1;
  } else if (lowerLimit == 0) {
    binomialCoefficients[shift] += 1;
  } else {
    computePositiveQBinomialHelper(binomialCoefficients, upperLimit - 1,
                                   lowerLimit, shift + lowerLimit);
    computePositiveQBinomialHelper(binomialCoefficients, upperLimit - 1,
                                   lowerLimit - 1, shift);
  }
}

void computePositiveQBinomial(std::vector<QPolynomialType> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg) {
  int maxQDegree = lowerLimit * (upperLimit - lowerLimit);
  std::vector<int> binomialCoefficients(maxQDegree + 1, 0);
  if (upperLimit == lowerLimit) {
    binomialCoefficients[0] = 1;
  } else if (lowerLimit == 0) {
    binomialCoefficients[0] = 1;
  } else {
    computePositiveQBinomialHelper(binomialCoefficients, upperLimit - 1,
                                   lowerLimit, lowerLimit);
    computePositiveQBinomialHelper(binomialCoefficients, upperLimit - 1,
                                   lowerLimit - 1, 0);
  }
  if (neg) {
    int componentSize = polynomialTerms[0].getComponentSize();
    if (componentSize <= 0) {
      throw std::runtime_error(
          "polynomialTerms[0] has invalid componentSize: " +
          std::to_string(componentSize));
    }
    QPolynomialType temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                  polynomialTerms[0].getPositiveSize(),
                                  componentSize, 0);
    for (int j = polynomialTerms[0].getMaxNegativeIndex();
         j <= polynomialTerms[0].getMaxPositiveIndex(); j++) {
      temporaryTerm[j] = polynomialTerms[0][j];
      polynomialTerms[0][j] = 0;
    }
    for (int j = temporaryTerm.getMaxNegativeIndex();
         j <= temporaryTerm.getMaxPositiveIndex(); j++) {
      if (temporaryTerm[j] != 0) {
        for (int k = 0; k < maxQDegree + 1; k++) {
          polynomialTerms[0][j - k] +=
              binomialCoefficients[k] * temporaryTerm[j];
        }
      }
    }
  } else {
    int componentSize = polynomialTerms[0].getComponentSize();
    if (componentSize <= 0) {
      throw std::runtime_error(
          "polynomialTerms[0] has invalid componentSize: " +
          std::to_string(componentSize));
    }
    QPolynomialType temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                  polynomialTerms[0].getPositiveSize(),
                                  componentSize, 0);
    for (int j = polynomialTerms[0].getMaxNegativeIndex();
         j <= polynomialTerms[0].getMaxPositiveIndex(); j++) {
      temporaryTerm[j] = polynomialTerms[0][j];
      polynomialTerms[0][j] = 0;
    }
    for (int j = temporaryTerm.getMaxNegativeIndex();
         j <= temporaryTerm.getMaxPositiveIndex(); j++) {
      if (temporaryTerm[j] != 0) {
        for (int k = 0; k < maxQDegree + 1; k++) {
          polynomialTerms[0][j + k] +=
              binomialCoefficients[k] * temporaryTerm[j];
        }
      }
    }
  }
}

void computeNegativeQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift,
                                    bool neg) {
  if (lowerLimit == 0) {
    if (neg) {
      binomialCoefficients[shift] += -1;
    } else {
      binomialCoefficients[shift] += 1;
    }
  } else if (lowerLimit < 0) {
    // Base case: if lowerLimit < 0, the binomial coefficient is 0
    // This prevents infinite recursion when lowerLimit keeps decreasing
    return;
  } else if (upperLimit == -1) {
    computeNegativeQBinomialHelper(binomialCoefficients, -1, lowerLimit - 1,
                                   shift - lowerLimit, !neg);
  } else {
    computeNegativeQBinomialHelper(binomialCoefficients, upperLimit,
                                   lowerLimit - 1,
                                   shift + 1 + upperLimit - lowerLimit, !neg);
    computeNegativeQBinomialHelper(binomialCoefficients, upperLimit + 1,
                                   lowerLimit, shift, neg);
  }
}

void computeNegativeQBinomial(std::vector<QPolynomialType> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg) {
  int qDegreeDelta = -(1 + upperLimit) * lowerLimit;
  int maxQDegree = -lowerLimit * (lowerLimit + 1) / 2;
  std::vector<int> binomialCoefficients(qDegreeDelta + 1, 0);
  if (lowerLimit == 0) {
    binomialCoefficients[0] = 1;
  } else if (upperLimit == -1) {
    computeNegativeQBinomialHelper(binomialCoefficients, -1, lowerLimit - 1,
                                   qDegreeDelta - maxQDegree - lowerLimit,
                                   true);
  } else {
    computeNegativeQBinomialHelper(
        binomialCoefficients, upperLimit, lowerLimit - 1,
        qDegreeDelta - maxQDegree + 1 + upperLimit - lowerLimit, true);
    computeNegativeQBinomialHelper(binomialCoefficients, upperLimit + 1,
                                   lowerLimit, qDegreeDelta - maxQDegree,
                                   false);
  }
  if (neg) {
    int componentSize = polynomialTerms[0].getComponentSize();
    if (componentSize <= 0) {
      throw std::runtime_error(
          "polynomialTerms[0] has invalid componentSize: " +
          std::to_string(componentSize));
    }
    QPolynomialType temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                  polynomialTerms[0].getPositiveSize(),
                                  componentSize, 0);
    for (int j = polynomialTerms[0].getMaxNegativeIndex();
         j <= polynomialTerms[0].getMaxPositiveIndex(); j++) {
      temporaryTerm[j] = polynomialTerms[0][j];
      polynomialTerms[0][j] = 0;
    }
    for (int j = temporaryTerm.getMaxNegativeIndex();
         j <= temporaryTerm.getMaxPositiveIndex(); j++) {
      if (temporaryTerm[j] != 0) {
        for (int k = 0; k < qDegreeDelta + 1; k++) {
          polynomialTerms[0][j - k + qDegreeDelta - maxQDegree] +=
              binomialCoefficients[k] * temporaryTerm[j];
        }
      }
    }
  } else {
    int componentSize = polynomialTerms[0].getComponentSize();
    if (componentSize <= 0) {
      throw std::runtime_error(
          "polynomialTerms[0] has invalid componentSize: " +
          std::to_string(componentSize));
    }
    QPolynomialType temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                  polynomialTerms[0].getPositiveSize(),
                                  componentSize, 0);
    for (int j = polynomialTerms[0].getMaxNegativeIndex();
         j <= polynomialTerms[0].getMaxPositiveIndex(); j++) {
      temporaryTerm[j] = polynomialTerms[0][j];
      polynomialTerms[0][j] = 0;
    }
    for (int j = temporaryTerm.getMaxNegativeIndex();
         j <= temporaryTerm.getMaxPositiveIndex(); j++) {
      if (temporaryTerm[j] != 0) {
        for (int k = 0; k < qDegreeDelta + 1; k++) {
          polynomialTerms[0][j + k - qDegreeDelta + maxQDegree] +=
              binomialCoefficients[k] * temporaryTerm[j];
        }
      }
    }
  }
}

void computeXQPochhammer(std::vector<QPolynomialType> &polynomialTerms,
                         int upperBound, int lowerBound, int componentIndex,
                         int totalComponents, std::vector<int> componentLengths,
                         std::vector<int> blockSizes) {
  for (int iterationVariable = lowerBound; iterationVariable <= upperBound;
       iterationVariable++) {
    for (int widthVariable = componentLengths[componentIndex];
         widthVariable > 0; widthVariable--) {
      matrixIndexColumn(totalComponents, componentLengths, componentIndex,
                        widthVariable - 1, polynomialTerms, 1,
                        iterationVariable, -1, blockSizes);
    }
  }
}

void computeXQInversePochhammer(std::vector<QPolynomialType> &polynomialTerms,
                                int upperBound, int lowerBound,
                                int componentIndex, int totalComponents,
                                std::vector<int> componentLengths,
                                std::vector<int> blockSizes) {
  for (int iterationVariable = lowerBound; iterationVariable <= upperBound;
       iterationVariable++) {
    for (int widthVariable = componentLengths[componentIndex];
         widthVariable > 0; widthVariable--) {
      for (int rankVariable = 1; rankVariable <= widthVariable;
           rankVariable++) {
        matrixIndexColumn(totalComponents, componentLengths, componentIndex,
                          widthVariable - rankVariable, polynomialTerms,
                          rankVariable, iterationVariable, 1, blockSizes);
      }
    }
  }
}

QPolynomialType QBinomialPositive(int upperLimit, int lowerLimit) {
  int n = upperLimit;
  int k = lowerLimit;

  // Zero polynomial if out of range
  if (k < 0 || n < 0 || k > n) {
    return QPolynomialType(0, 1, 1, 0); // all zero, only exponent 0 allocated
  }

  auto makeConstOne = []() {
    QPolynomialType p(0, 1, 1, 0); // exponent 0 only
    p[0] = 1;
    return p;
  };

  // C[i][j] = [i choose j]_q
  std::vector<std::vector<QPolynomialType>> C(
      n + 1, std::vector<QPolynomialType>(k + 1, QPolynomialType(0, 1, 1, 0)));

  C[0][0] = makeConstOne();

  for (int i = 1; i <= n; ++i) {
    C[i][0] = makeConstOne(); // [i choose 0]_q = 1

    int jMax = std::min(i, k);
    for (int j = 1; j <= jMax; ++j) {
      if (j == i) {
        C[i][j] = makeConstOne(); // [i choose i]_q = 1
      } else {
        // [i choose j]_q = [i-1 choose j]_q + q^(i-j) [i-1 choose j-1]_q
        QPolynomialType shifted = multiplyByQPower(C[i - 1][j - 1], i - j);
        C[i][j] = C[i - 1][j] + shifted;
      }
    }
  }

  return C[n][k];
}

/*
QPolynomialType QBinomialNegative(int upperLimit, int lowerLimit) {
  int k = lowerLimit;
  int u = upperLimit;

  // If k is out of range, return the zero polynomial
  if (k < 0) {
    return QPolynomialType(0, 1, 1, 0); // zero, only exponent 0 allocated
  }

  // If upperLimit is nonnegative, just defer to the positive version
  if (u >= 0) {
    return QBinomialPositive(u, k);
  }

  // Here upperLimit is negative: u = -n, with n > 0
  int n = -u;

  // Use the identity:
  //   [ -n choose k ]_q = (-1)^k q^(-n*k + k*(k-1)/2) [ n + k - 1 choose k ]_q
  //
  // First compute the positive q-binomial [ n + k - 1 choose k ]_q
  QPolynomialType base = QBinomialPositive(n + k - 1, k);

  // Compute the q-exponent shift: -n*k + k*(k-1)/2
  int shift = -n * k + (k * (k - 1)) / 2;

  // Apply the q^shift factor
  QPolynomialType result = multiplyByQPower(base, shift);

  // Apply the (-1)^k factor to all coefficients
  if (k % 2 != 0) { // k odd → multiply by -1
    int minExp = result.getMaxNegativeIndex();
    int maxExp = result.getMaxPositiveIndex();
    for (int e = minExp; e <= maxExp; ++e) {
      result[e] = -result[e];
    }
  }

  return result;
}*/

QPolynomialType QBinomialNegative(int upperLimit, int lowerLimit) {
  int k = lowerLimit;
  int u = upperLimit;

  // Handle nonsense k
  if (k < 0) {
    return QPolynomialType(0, 1, 1, 0); // zero polynomial
  }

  // If upperLimit >= 0, just reuse the positive version
  if (u >= 0) {
    return QBinomialPositive(u, k);
  }

  // u < 0: write u = -n with n > 0
  int n = -u;

  // Base positive q-binomial: [n + k - 1 choose k]_q
  QPolynomialType base = QBinomialPositive(n + k - 1, k);

  // Shift exponent so that we match computeNegativeQBinomial.
  // Desired shift: u * k - k*(k-1)/2
  // (this gives exponents [- (1+u)k - k(k+1)/2, -k(k+1)/2],
  //  which is exactly what your old code produces)
  int shift = u * k - (k * (k - 1)) / 2;

  QPolynomialType result = multiplyByQPower(base, shift);

  // Apply (-1)^k factor
  if (k % 2 != 0) { // k odd → multiply by -1
    int minExp = result.getMaxNegativeIndex();
    int maxExp = result.getMaxPositiveIndex();
    for (int e = minExp; e <= maxExp; ++e) {
      result[e] = -result[e];
    }
  }

  return result;
}

QPolynomialType QBinomial(int upperLimit, int lowerLimit) {
  return (upperLimit > 0) ? QBinomialPositive(upperLimit, lowerLimit)
                          : QBinomialNegative(upperLimit, lowerLimit);
}

// Compute (x q; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = ∏_{k=1}^n (1 - x q^{qpow + k})
// Optimized using direct coefficient computation
PolynomialType qpochhammer_xq_q(int n, int qpow) {
  const int numXVars = 1;
  const int maxXDegree = n;

  // Flattened coefficients map: coeffs[(x_degree, q_power)] = coefficient
  // Using pair keys instead of nested maps reduces allocations and improves cache locality
  std::map<std::pair<int, int>, int> coeffs;
  coeffs[{0, 0}] = 1; // Initialize with 1

  // For each factor (1 - x q^{qpow + k})
  for (int k = 0; k < n; ++k) {
    const int q_factor = qpow + k;
    std::map<std::pair<int, int>, int> new_coeffs;

    // Multiply current polynomial by (1 - x q^q_factor)
    for (const auto &[key, coeff] : coeffs) {
      const int x_deg = key.first;
      const int q_pow = key.second;

      if (coeff != 0) {
        // "1" term: contribute to same x_deg, q_pow
        new_coeffs[{x_deg, q_pow}] += coeff;

        // "-x q^q_factor" term: increment x_deg and add q_factor to q_pow
        if (x_deg + 1 <= maxXDegree) {
          new_coeffs[{x_deg + 1, q_pow + q_factor}] -= coeff;
        }
      }
    }
    coeffs = std::move(new_coeffs);
  }

  // Build result polynomial
  PolynomialType result(numXVars, maxXDegree);
  for (const auto &[key, coeff] : coeffs) {
    const int x_deg = key.first;
    const int q_pow = key.second;

    if (x_deg <= maxXDegree && coeff != 0) {
      result.addToCoefficient(q_pow, {x_deg}, coeff);
    }
  }

  return result;
}

// Compute 1/(x q^qpow; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = ∏_{l=0}^{n-1} ∑_{m=0}^{xMax+1} x^m q^{(l+qpow)*m}
// Optimized using direct coefficient computation


PolynomialType inverse_qpochhammer_xq_q(int n, int qpow, int xMax) {
  const int numXVars = 1;

  // Flattened coefficients map: coeffs[(x_degree, q_power)] = coefficient
  // Using pair keys instead of nested maps reduces allocations and improves cache locality
  std::map<std::pair<int, int>, int> coeffs;
  coeffs[{0, 0}] = 1; // Initialize with 1

  // For each factor (geometric series)
  for (int l = 0; l < n; ++l) {
    const int q_base = l + qpow;
    std::map<std::pair<int, int>, int> new_coeffs;

    // Multiply current polynomial by the l-th geometric series
    for (const auto &[key, coeff] : coeffs) {
      const int x_deg = key.first;
      const int q_pow = key.second;

      if (coeff != 0) {
        // Add terms from geometric series: 1 + x*q^q_base + x^2*q^(2*q_base) + ...
        for (int m = 0; x_deg + m <= xMax; ++m) {
          const int new_x_deg = x_deg + m;
          const int new_q_pow = q_pow + m * q_base;
          new_coeffs[{new_x_deg, new_q_pow}] += coeff;
        }
      }
    }
    coeffs = std::move(new_coeffs);
  }

  // Build result polynomial
  PolynomialType result(numXVars, 0);
  for (const auto &[key, coeff] : coeffs) {
    const int x_deg = key.first;
    const int q_pow = key.second;

    if (x_deg <= xMax && coeff != 0) {
      result.addToCoefficient(q_pow, {x_deg}, coeff);
    }
  }
  return result;
}
