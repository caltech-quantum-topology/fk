#include "fk/qalg_links.hpp"
#include "fk/linalg.hpp"
#include "fk/multivariable_polynomial.hpp"

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

void computePositiveQBinomial(std::vector<bilvector<int>> &polynomialTerms,
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
      throw std::runtime_error("polynomialTerms[0] has invalid componentSize: " + std::to_string(componentSize));
    }
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
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
      throw std::runtime_error("polynomialTerms[0] has invalid componentSize: " + std::to_string(componentSize));
    }
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
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

void computeNegativeQBinomial(std::vector<bilvector<int>> &polynomialTerms,
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
      throw std::runtime_error("polynomialTerms[0] has invalid componentSize: " + std::to_string(componentSize));
    }
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
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
      throw std::runtime_error("polynomialTerms[0] has invalid componentSize: " + std::to_string(componentSize));
    }
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
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

void computeXQPochhammer(std::vector<bilvector<int>> &polynomialTerms,
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

void computeXQInversePochhammer(std::vector<bilvector<int>> &polynomialTerms,
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



bilvector<int> QBinomialPositive(int upperLimit, int lowerLimit) {
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


/*
bilvector<int> QBinomialNegative(int upperLimit, int lowerLimit) {
  int k = lowerLimit;
  int u = upperLimit;

  // If k is out of range, return the zero polynomial
  if (k < 0) {
    return bilvector<int>(0, 1, 1, 0); // zero, only exponent 0 allocated
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
  bilvector<int> base = QBinomialPositive(n + k - 1, k);

  // Compute the q-exponent shift: -n*k + k*(k-1)/2
  int shift = -n * k + (k * (k - 1)) / 2;

  // Apply the q^shift factor
  bilvector<int> result = multiplyByQPower(base, shift);

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

bilvector<int> QBinomialNegative(int upperLimit, int lowerLimit) {
  int k = lowerLimit;
  int u = upperLimit;

  // Handle nonsense k
  if (k < 0) {
    return bilvector<int>(0, 1, 1, 0); // zero polynomial
  }

  // If upperLimit >= 0, just reuse the positive version
  if (u >= 0) {
    return QBinomialPositive(u, k);
  }

  // u < 0: write u = -n with n > 0
  int n = -u;

  // Base positive q-binomial: [n + k - 1 choose k]_q
  bilvector<int> base = QBinomialPositive(n + k - 1, k);

  // Shift exponent so that we match computeNegativeQBinomial.
  // Desired shift: u * k - k*(k-1)/2
  // (this gives exponents [- (1+u)k - k(k+1)/2, -k(k+1)/2],
  //  which is exactly what your old code produces)
  int shift = u * k - (k * (k - 1)) / 2;

  bilvector<int> result = multiplyByQPower(base, shift);

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

// Compute (x q; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = ∏_{k=1}^n (1 - x q^k)
MultivariablePolynomial qpochhammer_xq_q(int n, int qpow) {
    const int numXVars = 1;              // just x
    const int maxXDegree = static_cast<int>(n); // deg_x ≤ n

    // Start with P(q,x) = 1
    MultivariablePolynomial result(numXVars, maxXDegree);
    std::vector<int> zeroXPowers(numXVars, 0); // x^0
    result.addToCoefficient(0, zeroXPowers, 1); // q^0 * x^0 with coeff 1

    // Vector for x^1 (only x₁)
    std::vector<int> xXPowers(numXVars, 0);
    xXPowers[0] = 1; // x^1

    // Multiply factors (1 - x q^k) for k = 1..n
    for (int k = 1; k <= n; ++k) {
        MultivariablePolynomial factor(numXVars, 1);
        // 1 term: q^0 * x^0
        factor.addToCoefficient(0, zeroXPowers, 1);
        // - x q^k term: coefficient -1, q^k, x^1
        factor.addToCoefficient(qpow + k, xXPowers, -1);

        result *= factor;
    }

    return result;
}

// Compute 1/(x q^qpow; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = 1/∏_{k=1}^n (1 - x q^(k+qpow))
MultivariablePolynomial inverse_qpochhammer_xq_q(int n, int qpow, int xMax) {
  const int numXVars = 1;              // just x
  MultivariablePolynomial result(numXVars,0);
  result.setCoefficient(0,{0},1);
  for (int l = 0; l<n; ++l) {
     MultivariablePolynomial temp(numXVars,0);
     for (int m = 0; m <= std::min(n,xMax); ++m) {
        temp.setCoefficient((l+qpow)*m,{m},1);
    }
    result *= temp;
  }
  return result.truncate({xMax});
}
