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

void offsetAdditionRecursive(std::vector<bilvector<int>> &targetArray,
                             std::vector<bilvector<int>> &sourceArray,
                             std::vector<int> &offsetVector,
                             int bilvectorOffset, int &dimensions,
                             std::vector<int> arrayLengths, int currentIndex,
                             int accumulator, int accumulator2,
                             int signMultiplier, std::vector<int> targetBlocks,
                             std::vector<int> sourceBlocks) {
  currentIndex++;
  int oldAccumulator = accumulator;
  int oldAccumulator2 =
      accumulator2 + offsetVector[currentIndex] * targetBlocks[currentIndex];
  for (int i = std::max(0, -offsetVector[currentIndex]);
       i < arrayLengths[currentIndex] + 1; i++) {
    accumulator = oldAccumulator + i * sourceBlocks[currentIndex];
    accumulator2 = oldAccumulator2 + i * targetBlocks[currentIndex];
    if (dimensions > currentIndex + 1) {
      offsetAdditionRecursive(targetArray, sourceArray, offsetVector,
                              bilvectorOffset, dimensions, arrayLengths,
                              currentIndex, accumulator, accumulator2,
                              signMultiplier, targetBlocks, sourceBlocks);
    } else {
      for (int q = sourceArray[accumulator].getMaxNegativeIndex();
           q <= sourceArray[accumulator].getMaxPositiveIndex(); q++) {
        targetArray[accumulator2][q + bilvectorOffset] +=
            signMultiplier * sourceArray[accumulator][q];
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
  for (int i = std::max(0, -offsetVector[0]); i < arrayLengths[0] + 1; i++) {
    if (dimensions > 1) {
      offsetAdditionRecursive(targetArray, sourceArray, offsetVector,
                              bilvectorOffset, dimensions, arrayLengths, 0, i,
                              i + offsetVector[0], signMultiplier, targetBlocks,
                              sourceBlocks);
    } else {
      for (int q = sourceArray[i].getMaxNegativeIndex();
           q <= sourceArray[i].getMaxPositiveIndex(); q++) {
        targetArray[i + offsetVector[0]][q + bilvectorOffset] +=
            signMultiplier * sourceArray[i][q];
      }
    }
  }
}
