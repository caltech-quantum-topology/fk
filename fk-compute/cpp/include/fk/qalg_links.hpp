#pragma once

#include <functional>
#include <vector>

#include "fk/multivariable_polynomial.hpp"
#include "fk/bilvector.hpp"

// Forward declarations
void matrixIndexColumn(int &dimensions, std::vector<int> arrayLengths,
                       int &sliceIndex, int sliceValue,
                       std::vector<bilvector<int>> &polynomialTerms,
                       int rankOffset, int &zVariable, int signMultiplier,
                       std::vector<int> blockSizes);

// consider saving the multiplicands as separate variables before multiplying by
// term

void computePositiveQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift);

void computePositiveQBinomial(std::vector<bilvector<int>> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg);

void computeNegativeQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift,
                                    bool neg);

void computeNegativeQBinomial(std::vector<bilvector<int>> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg);

void computeXQPochhammer(std::vector<bilvector<int>> &polynomialTerms,
                         int upperBound, int lowerBound, int componentIndex,
                         int totalComponents, std::vector<int> componentLengths,
                         std::vector<int> blockSizes);

void computeXQInversePochhammer(std::vector<bilvector<int>> &polynomialTerms,
                                int upperBound, int lowerBound,
                                int componentIndex, int totalComponents,
                                std::vector<int> componentLengths,
                                std::vector<int> blockSizes);

bilvector<int> QBinomialPositive(int upperLimit, int lowerLimit);
bilvector<int> QBinomialNegative(int upperLimit, int lowerLimit);

MultivariablePolynomial qpochhammer_xq_q(int n, int qpow);
MultivariablePolynomial inverse_qpochhammer_xq_q(int n, int qpow, int xMax);
