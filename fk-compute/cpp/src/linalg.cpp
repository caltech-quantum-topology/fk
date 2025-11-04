#include "fk/linalg.hpp"
#include <iostream>

std::vector<int>
multiplyMatrixVector(std::vector<std::vector<int>> &inputMatrix,
                     std::vector<int> &inputVector) {
  std::vector<int> resultVector(inputVector.size());
  for (int i = 0; i < inputVector.size(); i++) {
    int accumulator = 0;
    for (int j = 0; j < inputVector.size(); j++) {
      accumulator += inputMatrix[i][j] * inputVector[j];
    }
    resultVector[i] = accumulator;
  }
  return resultVector;
}

std::vector<int>
multiplyMatrixVectorTranspose(std::vector<std::vector<int>> &inputMatrix,
                              const std::vector<int> &inputVector) {
  std::vector<int> resultVector(inputVector.size());
  for (int i = 0; i < inputVector.size(); i++) {
    int accumulator = 0;
    for (int j = 0; j < inputVector.size(); j++) {
      accumulator += inputMatrix[1 + j][1 + i] * inputVector[j];
    }
    resultVector[i] = accumulator;
  }
  return resultVector;
}

std::vector<double> mult(std::vector<std::vector<int>> &matrix,
                         std::vector<double> &vector) {
  std::vector<double> out(vector.size());
  for (int i = 0; i < vector.size(); i++) {
    int acc = 0;
    for (int j = 0; j < vector.size(); j++) {
      acc += matrix[i][j] * vector[j];
    }
    out[i] = acc;
  }
  return out;
}

int computeDotProduct(const std::vector<int> &a, const std::vector<int> &b) {
  int accumulator = a[0];
  for (int z = 0; z < b.size(); z++) {
    accumulator += a[z + 1] * b[z];
  }
  return accumulator;
}

void matrixIndexColumnRecursive(int &dimensions, std::vector<int> arrayLengths,
                                int &sliceIndex, int sliceValue,
                                int accumulator, int currentIndex,
                                std::vector<bilvector<int>> &polynomialTerms,
                                int rankOffset, int &zVariable,
                                int signMultiplier,
                                std::vector<int> blockSizes) {
  currentIndex++;
  int oldAccumulator = accumulator;
  if (sliceIndex == currentIndex) {
    accumulator = oldAccumulator + sliceValue * blockSizes[currentIndex];
    if (dimensions > currentIndex + 1) {
      matrixIndexColumnRecursive(dimensions, arrayLengths, sliceIndex,
                                 sliceValue, accumulator, currentIndex,
                                 polynomialTerms, rankOffset, zVariable,
                                 signMultiplier, blockSizes);
    } else {
      int d = accumulator + rankOffset * blockSizes[sliceIndex];
      for (int k = polynomialTerms[accumulator].getMaxNegativeIndex();
           k <= polynomialTerms[accumulator].getMaxPositiveIndex(); k++) {
        polynomialTerms[d][k + rankOffset * zVariable] +=
            signMultiplier * polynomialTerms[accumulator][k];
      }
    }
  } else {
    for (int i = 0; i < arrayLengths[currentIndex] + 1; i++) {
      accumulator = oldAccumulator + i * blockSizes[currentIndex];
      if (dimensions > currentIndex + 1) {
        matrixIndexColumnRecursive(dimensions, arrayLengths, sliceIndex,
                                   sliceValue, accumulator, currentIndex,
                                   polynomialTerms, rankOffset, zVariable,
                                   signMultiplier, blockSizes);
      } else {
        int d = accumulator + rankOffset * blockSizes[sliceIndex];
        for (int k = polynomialTerms[accumulator].getMaxNegativeIndex();
             k <= polynomialTerms[accumulator].getMaxPositiveIndex(); k++) {
          polynomialTerms[d][k + rankOffset * zVariable] +=
              signMultiplier * polynomialTerms[accumulator][k];
        }
      }
    }
  }
}

void matrixIndexColumn(int &dimensions, std::vector<int> arrayLengths,
                       int &sliceIndex, int sliceValue,
                       std::vector<bilvector<int>> &polynomialTerms,
                       int rankOffset, int &zVariable, int signMultiplier,
                       std::vector<int> blockSizes) {
  if (sliceIndex == 0) {
    if (dimensions > 1) {
      matrixIndexColumnRecursive(
          dimensions, arrayLengths, sliceIndex, sliceValue, sliceValue, 0,
          polynomialTerms, rankOffset, zVariable, signMultiplier, blockSizes);
    } else {
      int d = sliceValue + rankOffset;
      for (int k = polynomialTerms[sliceValue].getMaxNegativeIndex();
           k <= polynomialTerms[sliceValue].getMaxPositiveIndex(); k++) {
        polynomialTerms[d][k + rankOffset * zVariable] +=
            signMultiplier * polynomialTerms[sliceValue][k];
      }
    }
  } else {
    for (int i = 0; i < arrayLengths[0] + 1; i++) {
      if (dimensions > 1) {
        matrixIndexColumnRecursive(
            dimensions, arrayLengths, sliceIndex, sliceValue, i, 0,
            polynomialTerms, rankOffset, zVariable, signMultiplier, blockSizes);
      } else {
        int d = i + rankOffset;
        for (int k = polynomialTerms[i].getMaxNegativeIndex();
             k <= polynomialTerms[i].getMaxPositiveIndex(); k++) {
          polynomialTerms[d][k + rankOffset * zVariable] +=
              signMultiplier * polynomialTerms[i][k];
        }
      }
    }
  }
}

void performOffsetAddition(std::vector<bilvector<int>> &targetArray,
                           std::vector<bilvector<int>> sourceArray,
                           std::vector<int> &offsetVector, int bilvectorOffset,
                           int &dimensions, std::vector<int> arrayLengths,
                           int signMultiplier, std::vector<int> targetBlocks,
                           std::vector<int> sourceBlocks) {
  if (dimensions == 1) {
    for (int i = std::max(0, -offsetVector[0]); i < arrayLengths[0] + 1; i++) {
      for (int q = sourceArray[i].getMaxNegativeIndex();
           q <= sourceArray[i].getMaxPositiveIndex(); q++) {
        targetArray[i + offsetVector[0]][q + bilvectorOffset] +=
            signMultiplier * sourceArray[i][q];
      }
    }
    return;
  }
  std::vector<int> indices(dimensions);
  std::vector<int> limits(dimensions);

  for (int d = 0; d < dimensions; d++) {
    indices[d] = std::max(0, -offsetVector[d]);
    limits[d] = arrayLengths[d] + 1;
  }

  bool done = false;
  while (!done) {
    int sourceAccumulator = 0;
    int targetAccumulator = 0;

    for (int d = 0; d < dimensions; d++) {
      sourceAccumulator += indices[d] * sourceBlocks[d];
      targetAccumulator += (indices[d] + offsetVector[d]) * targetBlocks[d];
    }

    for (int q = sourceArray[sourceAccumulator].getMaxNegativeIndex();
         q <= sourceArray[sourceAccumulator].getMaxPositiveIndex(); q++) {
      targetArray[targetAccumulator][q + bilvectorOffset] +=
          signMultiplier * sourceArray[sourceAccumulator][q];
    }

    int d = dimensions - 1;
    while (d >= 0) {
      indices[d]++;
      if (indices[d] < limits[d]) {
        break;
      }
      indices[d] = std::max(0, -offsetVector[d]);
      d--;
    }

    if (d < 0) {
      done = true;
    }
  }
}
