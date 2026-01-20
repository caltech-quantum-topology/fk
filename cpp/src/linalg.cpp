#include "fk/linalg.hpp"
#include <iostream>
#include <algorithm>


int computeDotProduct(const std::vector<int> &a, const std::vector<int> &b) {
  int accumulator = a[0];
  for (size_t z = 0; z < b.size(); z++) {
    accumulator += a[z + 1] * b[z];
  }
  return accumulator;
}

void matrixIndexColumnRecursive(int &dimensions, std::vector<int> arrayLengths,
                                int &sliceIndex, int sliceValue,
                                int accumulator, int currentIndex,
                                std::vector<QPolynomialType> &polynomialTerms,
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
        int coeff = polynomialTerms[accumulator].getCoefficient(k);
        polynomialTerms[d].addToCoefficient(k + rankOffset * zVariable,
                                            signMultiplier * coeff);
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
          int coeff = polynomialTerms[accumulator].getCoefficient(k);
          polynomialTerms[d].addToCoefficient(k + rankOffset * zVariable,
                                              signMultiplier * coeff);
        }
      }
    }
  }
}

void matrixIndexColumn(int &dimensions, std::vector<int> arrayLengths,
                       int &sliceIndex, int sliceValue,
                       std::vector<QPolynomialType> &polynomialTerms,
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
        int coeff = polynomialTerms[sliceValue].getCoefficient(k);
        polynomialTerms[d].addToCoefficient(k + rankOffset * zVariable,
                                            signMultiplier * coeff);
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
          int coeff = polynomialTerms[i].getCoefficient(k);
          polynomialTerms[d].addToCoefficient(k + rankOffset * zVariable,
                                              signMultiplier * coeff);
        }
      }
    }
  }
}

void performOffsetAddition(std::vector<QPolynomialType> &targetArray,
                           std::vector<QPolynomialType> sourceArray,
                           std::vector<int> &offsetVector, int bilvectorOffset,
                           int &dimensions, std::vector<int> arrayLengths,
                           int signMultiplier, std::vector<int> targetBlocks,
                           std::vector<int> sourceBlocks) {
  if (dimensions == 1) {
    for (int i = std::max(0, -offsetVector[0]); i < arrayLengths[0] + 1; i++) {
      for (int q = sourceArray[i].getMaxNegativeIndex();
           q <= sourceArray[i].getMaxPositiveIndex(); q++) {
        int coeff = sourceArray[i].getCoefficient(q);
        targetArray[i + offsetVector[0]].addToCoefficient(
            q + bilvectorOffset, signMultiplier * coeff);
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
      int coeff = sourceArray[sourceAccumulator].getCoefficient(q);
      targetArray[targetAccumulator].addToCoefficient(
          q + bilvectorOffset, signMultiplier * coeff);
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

using Term = std::pair<std::vector<int>, QPolynomialType>;

void performOffsetAddition(
    std::vector<Term> &targetArray,
    const std::vector<Term> &sourceArray,
    const std::vector<int> &offsetVector,
    int bilvectorOffset,
    int signMultiplier,
    int dimensions,
    const std::vector<int> &arrayLengths) {
  // We assume:
  // - dimensions == arrayLengths.size()
  // - for each term, term.first.size() >= dimensions
  for (const auto &term : sourceArray) {
    const std::vector<int> &indices = term.first;
    const QPolynomialType   &srcBil  = term.second;

    if (static_cast<int>(indices.size()) < dimensions) {
      continue; // malformed, skip
    }

    // --- This is the key part: emulate dense iteration domain ---
    bool inDomain = true;
    for (int d = 0; d < dimensions; ++d) {
      int i    = indices[d];
      int low  = std::max(0, -offsetVector[d]);
      int high = arrayLengths[d];  // inclusive (because loop went to < +1)

      if (i < low || i > high) {
        inDomain = false;
        break;
      }
    }
    if (!inDomain) {
      continue; // dense code would never touch this index
    }

    // Compute shifted x-degrees (same as before)
    std::vector<int> shifted(dimensions);
    for (int d = 0; d < dimensions; ++d) {
      shifted[d] = indices[d] + offsetVector[d];
    }

    // Find (or create) matching term in target
    auto it = std::find_if(
        targetArray.begin(), targetArray.end(),
        [&](const Term &t) { return t.first == shifted; });

    if (it == targetArray.end()) {
      targetArray.emplace_back(shifted, QPolynomialType());
      it = std::prev(targetArray.end());
    }

    QPolynomialType &destBil = it->second;

    // Add shifted q-coefficients
    for (int q = srcBil.getMaxNegativeIndex();
         q <= srcBil.getMaxPositiveIndex(); ++q) {

      int coeff = srcBil.getCoefficient(q);
      if (coeff != 0) {
        destBil.addToCoefficient(q + bilvectorOffset, signMultiplier * coeff);
      }
    }
  }
}
