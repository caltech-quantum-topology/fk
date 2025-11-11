#pragma once

#include <map>
#include <vector>

#include "fk/polynomial_config.hpp"

// Forward declarations
void matrixIndexColumn(int &dimensions, std::vector<int> arrayLengths,
                       int &sliceIndex, int sliceValue,
                       std::vector<QPolynomialType> &polynomialTerms,
                       int rankOffset, int &zVariable, int signMultiplier,
                       std::vector<int> blockSizes);

// consider saving the multiplicands as separate variables before multiplying by
// term

void computePositiveQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift);

void computePositiveQBinomial(std::vector<QPolynomialType> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg);

void computeNegativeQBinomialHelper(std::vector<int> &binomialCoefficients,
                                    int upperLimit, int lowerLimit, int shift,
                                    bool neg);

void computeNegativeQBinomial(std::vector<QPolynomialType> &polynomialTerms,
                              int upperLimit, int lowerLimit, bool neg);

void computeXQPochhammer(std::vector<QPolynomialType> &polynomialTerms,
                         int upperBound, int lowerBound, int componentIndex,
                         int totalComponents, std::vector<int> componentLengths,
                         std::vector<int> blockSizes);

void computeXQInversePochhammer(std::vector<QPolynomialType> &polynomialTerms,
                                int upperBound, int lowerBound,
                                int componentIndex, int totalComponents,
                                std::vector<int> componentLengths,
                                std::vector<int> blockSizes);

QPolynomialType QBinomialPositive(int upperLimit, int lowerLimit);
QPolynomialType QBinomialNegative(int upperLimit, int lowerLimit);

PolynomialType qpochhammer_xq_q(int n, int qpow, int lsign);
PolynomialType inverse_qpochhammer_xq_q(int n, int qpow, int xMax, int lsign,
                                        bool xInverse);
std::map<int, std::map<int, int>> inv_qpoch_xq_q_(int n, int qpow, int xMax,
                                                  int lsign);
std::map<int, std::map<int, int>> inv_qpoch_xinvq_q_(int n, int qpow, int xMax,
                                                     int lsign);
