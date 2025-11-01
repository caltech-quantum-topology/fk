#pragma once

#include "fk/bilvector.hpp"
#include <vector>

// x @ p <= xdeg <------- original inequality in angle variables; use a
// change-of-basis matrix M to examine the segment side (M^-1)^T p @ M(angle
// x-criterion)  <= xdeg <------- same inequality in segment variables; solve
// for (M^-1)^T using the transformed criterion p = M^T (M^-1)^T p <-------
// relation between respective points of the integer hulls; interestingly, M^-1
// isn't ever invoked explicitly

// when both ILP's are put in standard form, the nontrivial inequalities of a
// given side just becomes the trivial identity matrix expressing the sign of
// the variables in the other side

std::vector<int>
multiplyMatrixVector(std::vector<std::vector<int>> &inputMatrix,
                     std::vector<int> &inputVector);

// with transpose
std::vector<int>
multiplyMatrixVectorTranspose(std::vector<std::vector<int>> &inputMatrix,
                              const std::vector<int> &inputVector);

std::vector<double> mult(std::vector<std::vector<int>> &matrix,
                         std::vector<double> &vector);

int computeDotProduct(const std::vector<int> &a, const std::vector<int> &b);

void matrixIndexColumnRecursive(int &dimensions, std::vector<int> arrayLengths,
                                int &sliceIndex, int sliceValue,
                                int accumulator, int currentIndex,
                                std::vector<bilvector<int>> &polynomialTerms,
                                int rankOffset, int &zVariable,
                                int signMultiplier,
                                std::vector<int> blockSizes);

void matrixIndexColumn(int &dimensions, std::vector<int> arrayLengths,
                       int &sliceIndex, int sliceValue,
                       std::vector<bilvector<int>> &polynomialTerms,
                       int rankOffset, int &zVariable, int signMultiplier,
                       std::vector<int> blockSizes);

void offsetAdditionRecursive(std::vector<bilvector<int>> &targetArray,
                             std::vector<bilvector<int>> &sourceArray,
                             std::vector<int> &offsetVector,
                             int bilvectorOffset, int &dimensions,
                             std::vector<int> arrayLengths, int currentIndex,
                             int accumulator, int accumulator2,
                             int signMultiplier, std::vector<int> targetBlocks,
                             std::vector<int> sourceBlocks);

void performOffsetAddition(std::vector<bilvector<int>> &targetArray,
                           std::vector<bilvector<int>> sourceArray,
                           std::vector<int> &offsetVector, int bilvectorOffset,
                           int &dimensions, std::vector<int> arrayLengths,
                           int signMultiplier, std::vector<int> targetBlocks,
                           std::vector<int> sourceBlocks);