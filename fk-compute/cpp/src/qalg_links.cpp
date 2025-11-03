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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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

MultivariablePolynomial computePositiveQBinomial(int upperLimit, int lowerLimit,
                                                 bool neg) {
  int maxQDegree = lowerLimit * (upperLimit - lowerLimit);
  MultivariablePolynomial binomialCoefficients();
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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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
    bilvector<int> temporaryTerm(polynomialTerms[0].getNegativeSize(),
                                 polynomialTerms[0].getPositiveSize(),
                                 polynomialTerms[0].getComponentSize(), 0);
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
