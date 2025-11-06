#include "fk/bmpoly.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <cmath>

// BMPoly implementation

BMPoly::BMPoly(int numVariables, int degree, const std::vector<int> &maxDegrees)
    : numXVariables(numVariables) {

  if (numVariables < 0) {
    throw std::invalid_argument("Number of variables must be non-negative");
  }

  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(numVariables, degree);
  } else {
    if (maxDegrees.size() != static_cast<size_t>(numVariables)) {
      throw std::invalid_argument(
          "Max degrees vector size must match number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  // Initialize ground degrees to 0 (no negative exponents initially)
  groundXDegrees = std::vector<int>(numVariables, 0);

  // Calculate initial block sizes and setup storage
  recalculateBlockSizes();
  resizeStorage();
}

BMPoly::BMPoly(const BMPoly &source, int newNumVariables,
               int targetVariableIndex, int degree,
               const std::vector<int> &maxDegrees)
    : numXVariables(newNumVariables) {

  if (newNumVariables < source.numXVariables) {
    throw std::invalid_argument(
        "New number of variables must be >= source's number of variables");
  }

  if (targetVariableIndex < 0 || targetVariableIndex >= newNumVariables) {
    throw std::invalid_argument(
        "Target variable index must be in range [0, newNumVariables)");
  }

  if (source.numXVariables != 1) {
    throw std::invalid_argument(
        "Source polynomial must have exactly 1 x variable");
  }

  // Set up max degrees
  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(newNumVariables, degree);
  } else {
    if (maxDegrees.size() != static_cast<size_t>(newNumVariables)) {
      throw std::invalid_argument(
          "Max degrees vector size must match new number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  // Initialize ground degrees
  groundXDegrees = std::vector<int>(newNumVariables, 0);
  groundXDegrees[targetVariableIndex] = source.groundXDegrees[0];

  // Calculate block sizes and setup storage
  recalculateBlockSizes();
  resizeStorage();

  // Copy coefficients from source, mapping the single variable to targetVariableIndex
  const auto &sourceCoeffs = source.coeffs_;
  for (size_t i = 0; i < sourceCoeffs.size(); ++i) {
    if (!sourceCoeffs[i].isZero()) {
      // Convert linear index back to multi-index for source
      std::vector<int> sourceXPowers = source.linearToMultiIndex(i);

      // Create new x-powers vector with zeros except at targetVariableIndex
      std::vector<int> newXPowers(newNumVariables, 0);
      newXPowers[targetVariableIndex] = sourceXPowers[0];

      // Copy the bilvector to the new location
      int newLinearIndex = multiIndexToLinear(newXPowers);
      if (newLinearIndex >= 0 && newLinearIndex < static_cast<int>(coeffs_.size())) {
        coeffs_[newLinearIndex] = sourceCoeffs[i];
      }
    }
  }
}

int BMPoly::multiIndexToLinear(const std::vector<int> &xPowers) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("xPowers size must match number of variables");
  }

  int linearIndex = 0;
  for (int i = 0; i < numXVariables; ++i) {
    int adjustedPower = xPowers[i] - groundXDegrees[i];
    if (adjustedPower < 0 || adjustedPower > maxXDegrees[i]) {
      return -1; // Out of bounds
    }
    linearIndex += adjustedPower * blockSizes[i];
  }
  return linearIndex;
}

std::vector<int> BMPoly::linearToMultiIndex(int linearIndex) const {
  std::vector<int> xPowers(numXVariables);

  for (int i = numXVariables - 1; i >= 0; --i) {
    int adjustedPower = linearIndex / blockSizes[i];
    xPowers[i] = adjustedPower + groundXDegrees[i];
    linearIndex %= blockSizes[i];
  }

  return xPowers;
}

void BMPoly::expandStorageIfNeeded(const std::vector<int> &xPowers) {
  bool needsExpansion = false;
  bool needsGroundUpdate = false;

  // Check if any power is outside current bounds
  for (int i = 0; i < numXVariables; ++i) {
    if (xPowers[i] < groundXDegrees[i]) {
      needsGroundUpdate = true;
    }
    if (xPowers[i] > groundXDegrees[i] + maxXDegrees[i]) {
      needsExpansion = true;
    }
  }

  if (needsGroundUpdate) {
    updateGroundDegrees(xPowers);
    needsExpansion = true; // Ground update requires storage resize
  }

  if (needsExpansion) {
    // Expand maxXDegrees to accommodate new powers
    for (int i = 0; i < numXVariables; ++i) {
      int requiredRange = xPowers[i] - groundXDegrees[i];
      if (requiredRange > maxXDegrees[i]) {
        maxXDegrees[i] = requiredRange;
      }
    }

    recalculateBlockSizes();
    resizeStorage();
  }
}

void BMPoly::updateGroundDegrees(const std::vector<int> &xPowers) {
  std::vector<bilvector<int>> oldCoeffs = std::move(coeffs_);
  std::vector<int> oldGroundDegrees = groundXDegrees;

  // Update ground degrees
  for (int i = 0; i < numXVariables; ++i) {
    groundXDegrees[i] = std::min(groundXDegrees[i], xPowers[i]);
  }

  // Recalculate block sizes and resize storage
  recalculateBlockSizes();
  resizeStorage();

  // Migrate old coefficients to new positions
  for (size_t oldIndex = 0; oldIndex < oldCoeffs.size(); ++oldIndex) {
    if (!oldCoeffs[oldIndex].isZero()) {
      // Convert old linear index to multi-index using old ground degrees
      std::vector<int> multiIndex(numXVariables);
      int tempIndex = oldIndex;

      for (int i = numXVariables - 1; i >= 0; --i) {
        int oldBlockSize = (i == 0) ? 1 : oldBlockSize * (maxXDegrees[i-1] + 1);
        int adjustedPower = tempIndex / oldBlockSize;
        multiIndex[i] = adjustedPower + oldGroundDegrees[i];
        tempIndex %= oldBlockSize;
      }

      // Find new linear index for this multi-index
      int newIndex = multiIndexToLinear(multiIndex);
      if (newIndex >= 0 && newIndex < static_cast<int>(coeffs_.size())) {
        coeffs_[newIndex] = std::move(oldCoeffs[oldIndex]);
      }
    }
  }
}

void BMPoly::recalculateBlockSizes() {
  blockSizes.resize(numXVariables);
  if (numXVariables > 0) {
    blockSizes[0] = 1;
    for (int i = 1; i < numXVariables; i++) {
      blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
    }
  }
}

void BMPoly::resizeStorage() {
  int totalSize = 1;
  for (int i = 0; i < numXVariables; ++i) {
    totalSize *= (maxXDegrees[i] + 1);
  }

  // Resize vector, preserving existing data
  size_t oldSize = coeffs_.size();
  coeffs_.resize(totalSize, bilvector<int>(0, 1, 20, 0));

  // Initialize new bilvectors if needed
  for (size_t i = oldSize; i < coeffs_.size(); ++i) {
    coeffs_[i] = bilvector<int>(0, 1, 20, 0);
  }
}

void BMPoly::checkIndexBounds(const std::vector<int> &xPowers) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("xPowers size must match number of variables");
  }
}

void BMPoly::pruneZeros() {
  // For vector-based storage, we can't really "prune" but we can detect empty bilvectors
  // This is mainly for consistency with the interface
}

int BMPoly::getCoefficient(int qPower, const std::vector<int> &xPowers) const {
  checkIndexBounds(xPowers);

  int linearIndex = multiIndexToLinear(xPowers);
  if (linearIndex < 0 || linearIndex >= static_cast<int>(coeffs_.size())) {
    return 0; // Out of bounds
  }

  return coeffs_[linearIndex][qPower];
}

void BMPoly::setCoefficient(int qPower, const std::vector<int> &xPowers,
                            int coefficient) {
  checkIndexBounds(xPowers);

  // Expand storage if needed
  const_cast<BMPoly*>(this)->expandStorageIfNeeded(xPowers);

  int linearIndex = multiIndexToLinear(xPowers);
  if (linearIndex >= 0 && linearIndex < static_cast<int>(coeffs_.size())) {
    coeffs_[linearIndex][qPower] = coefficient;
  }
}

void BMPoly::addToCoefficient(int qPower, const std::vector<int> &xPowers,
                              int coefficient) {
  checkIndexBounds(xPowers);

  // Expand storage if needed
  expandStorageIfNeeded(xPowers);

  int linearIndex = multiIndexToLinear(xPowers);
  if (linearIndex >= 0 && linearIndex < static_cast<int>(coeffs_.size())) {
    coeffs_[linearIndex][qPower] += coefficient;
  }
}

bilvector<int> &BMPoly::getQPolynomial(const std::vector<int> &xPowers) {
  checkIndexBounds(xPowers);

  // Expand storage if needed
  expandStorageIfNeeded(xPowers);

  int linearIndex = multiIndexToLinear(xPowers);
  if (linearIndex < 0 || linearIndex >= static_cast<int>(coeffs_.size())) {
    throw std::out_of_range("Invalid x-powers for getQPolynomial");
  }

  return coeffs_[linearIndex];
}

const bilvector<int> &BMPoly::getQPolynomial(const std::vector<int> &xPowers) const {
  checkIndexBounds(xPowers);

  int linearIndex = multiIndexToLinear(xPowers);
  if (linearIndex < 0 || linearIndex >= static_cast<int>(coeffs_.size())) {
    throw std::out_of_range("Invalid x-powers for getQPolynomial");
  }

  return coeffs_[linearIndex];
}

BMPoly BMPoly::invertVariable(const int target_index) {
  if (target_index < 0 || target_index >= numXVariables) {
    throw std::out_of_range("Variable index out of range");
  }

  BMPoly result(numXVariables, 0, maxXDegrees); // Create with same structure

  // Iterate through all stored coefficients
  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero()) {
      std::vector<int> xPowers = linearToMultiIndex(i);

      // Invert the target variable: x_i -> x_i^(-1), so power becomes -power
      xPowers[target_index] = -xPowers[target_index];

      // Copy the bilvector to the result
      result.expandStorageIfNeeded(xPowers);
      int resultIndex = result.multiIndexToLinear(xPowers);
      if (resultIndex >= 0 && resultIndex < static_cast<int>(result.coeffs_.size())) {
        result.coeffs_[resultIndex] = coeffs_[i];
      }
    }
  }

  return result;
}

BMPoly BMPoly::truncate(const std::vector<int> &maxXdegrees) {
  if (maxXdegrees.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("maxXdegrees size must match number of variables");
  }

  BMPoly result(numXVariables, 0, maxXdegrees);

  // Copy coefficients that fall within the truncation bounds
  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero()) {
      std::vector<int> xPowers = linearToMultiIndex(i);

      // Check if this term should be included
      bool include = true;
      for (int j = 0; j < numXVariables; ++j) {
        if (xPowers[j] > maxXdegrees[j]) {
          include = false;
          break;
        }
      }

      if (include) {
        result.expandStorageIfNeeded(xPowers);
        int resultIndex = result.multiIndexToLinear(xPowers);
        if (resultIndex >= 0 && resultIndex < static_cast<int>(result.coeffs_.size())) {
          result.coeffs_[resultIndex] = coeffs_[i];
        }
      }
    }
  }

  return result;
}

std::vector<bilvector<int>> BMPoly::getCoefficients() const {
  return coeffs_;
}

std::vector<bilvector<int>> &BMPoly::getCoefficients() {
  return coeffs_;
}

void BMPoly::syncFromDenseVector(const std::vector<bilvector<int>> &denseVector) {
  if (denseVector.size() != coeffs_.size()) {
    throw std::invalid_argument("Dense vector size must match internal storage size");
  }
  coeffs_ = denseVector;
}

int BMPoly::getNumXVariables() const {
  return numXVariables;
}

const std::vector<int> &BMPoly::getMaxXDegrees() const {
  return maxXDegrees;
}

const std::vector<int> &BMPoly::getBlockSizes() const {
  return blockSizes;
}

void BMPoly::clear() {
  for (auto &bilvec : coeffs_) {
    bilvec = bilvector<int>(0, 1, 20, 0);
  }
}

bool BMPoly::isZero() const {
  for (const auto &bilvec : coeffs_) {
    if (!bilvec.isZero()) {
      return false;
    }
  }
  return true;
}

void BMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream file(fileName);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file for writing: " + fileName);
  }

  file << "{\n";
  file << "  \"type\": \"BMPoly\",\n";
  file << "  \"numXVariables\": " << numXVariables << ",\n";
  file << "  \"maxXDegrees\": [";
  for (size_t i = 0; i < maxXDegrees.size(); ++i) {
    if (i > 0) file << ", ";
    file << maxXDegrees[i];
  }
  file << "],\n";
  file << "  \"groundXDegrees\": [";
  for (size_t i = 0; i < groundXDegrees.size(); ++i) {
    if (i > 0) file << ", ";
    file << groundXDegrees[i];
  }
  file << "],\n";
  file << "  \"coefficients\": {\n";

  bool first = true;
  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero()) {
      if (!first) file << ",\n";
      first = false;

      std::vector<int> xPowers = linearToMultiIndex(i);
      file << "    \"[";
      for (size_t j = 0; j < xPowers.size(); ++j) {
        if (j > 0) file << ",";
        file << xPowers[j];
      }
      file << "]\": {";

      // Export non-zero q-coefficients
      bool firstQ = true;
      int minQ = coeffs_[i].getMaxNegativeIndex();
      int maxQ = coeffs_[i].getMaxPositiveIndex();

      for (int q = minQ; q <= maxQ; ++q) {
        int coeff = coeffs_[i][q];
        if (coeff != 0) {
          if (!firstQ) file << ", ";
          firstQ = false;
          file << "\"" << q << "\": " << coeff;
        }
      }
      file << "}";
    }
  }

  file << "\n  }\n";
  file << "}\n";
  file.close();
}

void BMPoly::print(int maxTerms) const {
  std::cout << "BMPoly(" << numXVariables << " variables):\n";

  int termCount = 0;
  for (size_t i = 0; i < coeffs_.size() && termCount < maxTerms; ++i) {
    if (!coeffs_[i].isZero()) {
      std::vector<int> xPowers = linearToMultiIndex(i);

      // Print x-monomial
      std::cout << "  x^[";
      for (size_t j = 0; j < xPowers.size(); ++j) {
        if (j > 0) std::cout << ",";
        std::cout << xPowers[j];
      }
      std::cout << "]: ";

      // Print q-polynomial (simplified)
      int minQ = coeffs_[i].getMaxNegativeIndex();
      int maxQ = coeffs_[i].getMaxPositiveIndex();
      bool firstTerm = true;

      for (int q = minQ; q <= maxQ; ++q) {
        int coeff = coeffs_[i][q];
        if (coeff != 0) {
          if (!firstTerm) std::cout << " + ";
          firstTerm = false;
          std::cout << coeff << "*q^" << q;
        }
      }
      std::cout << "\n";
      termCount++;
    }
  }

  if (termCount == maxTerms && termCount < static_cast<int>(coeffs_.size())) {
    std::cout << "  ... (showing first " << maxTerms << " terms)\n";
  }
}

bilvector<int> BMPoly::evaluate(const std::vector<int> &point) const {
  if (point.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("Point dimension must match number of variables");
  }

  bilvector<int> result(0, 1, 20, 0);

  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero()) {
      std::vector<int> xPowers = linearToMultiIndex(i);

      // Calculate x-monomial value
      int monomialValue = 1;
      for (int j = 0; j < numXVariables; ++j) {
        if (xPowers[j] < 0 && point[j] == 0) {
          throw std::domain_error("Division by zero: negative exponent with zero value");
        }
        monomialValue *= static_cast<int>(std::pow(point[j], xPowers[j]));
      }

      // Add contribution to result
      int minQ = coeffs_[i].getMaxNegativeIndex();
      int maxQ = coeffs_[i].getMaxPositiveIndex();

      for (int q = minQ; q <= maxQ; ++q) {
        int coeff = coeffs_[i][q];
        if (coeff != 0) {
          result[q] += coeff * monomialValue;
        }
      }
    }
  }

  return result;
}

void BMPoly::checkCompatibility(const BMPoly &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument("Polynomials must have the same number of variables");
  }
}

BMPoly &BMPoly::operator+=(const BMPoly &other) {
  checkCompatibility(other);

  // Ensure our storage can accommodate all terms from other
  for (size_t i = 0; i < other.coeffs_.size(); ++i) {
    if (!other.coeffs_[i].isZero()) {
      std::vector<int> xPowers = other.linearToMultiIndex(i);
      expandStorageIfNeeded(xPowers);

      int ourIndex = multiIndexToLinear(xPowers);
      if (ourIndex >= 0 && ourIndex < static_cast<int>(coeffs_.size())) {
        // Add corresponding q-polynomials
        int minQ = std::min(coeffs_[ourIndex].getMaxNegativeIndex(),
                           other.coeffs_[i].getMaxNegativeIndex());
        int maxQ = std::max(coeffs_[ourIndex].getMaxPositiveIndex(),
                           other.coeffs_[i].getMaxPositiveIndex());

        for (int q = minQ; q <= maxQ; ++q) {
          coeffs_[ourIndex][q] += other.coeffs_[i][q];
        }
      }
    }
  }

  return *this;
}

BMPoly &BMPoly::operator-=(const BMPoly &other) {
  checkCompatibility(other);

  // Ensure our storage can accommodate all terms from other
  for (size_t i = 0; i < other.coeffs_.size(); ++i) {
    if (!other.coeffs_[i].isZero()) {
      std::vector<int> xPowers = other.linearToMultiIndex(i);
      expandStorageIfNeeded(xPowers);

      int ourIndex = multiIndexToLinear(xPowers);
      if (ourIndex >= 0 && ourIndex < static_cast<int>(coeffs_.size())) {
        // Subtract corresponding q-polynomials
        int minQ = std::min(coeffs_[ourIndex].getMaxNegativeIndex(),
                           other.coeffs_[i].getMaxNegativeIndex());
        int maxQ = std::max(coeffs_[ourIndex].getMaxPositiveIndex(),
                           other.coeffs_[i].getMaxPositiveIndex());

        for (int q = minQ; q <= maxQ; ++q) {
          coeffs_[ourIndex][q] -= other.coeffs_[i][q];
        }
      }
    }
  }

  return *this;
}

BMPoly &BMPoly::operator*=(const BMPoly &other) {
  checkCompatibility(other);

  BMPoly result(numXVariables, 0, maxXDegrees);

  // Multiply each term in this polynomial with each term in other
  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero()) {
      std::vector<int> xPowers1 = linearToMultiIndex(i);

      for (size_t j = 0; j < other.coeffs_.size(); ++j) {
        if (!other.coeffs_[j].isZero()) {
          std::vector<int> xPowers2 = other.linearToMultiIndex(j);

          // Calculate result x-powers (add exponents)
          std::vector<int> resultXPowers(numXVariables);
          for (int k = 0; k < numXVariables; ++k) {
            resultXPowers[k] = xPowers1[k] + xPowers2[k];
          }

          result.expandStorageIfNeeded(resultXPowers);
          int resultIndex = result.multiIndexToLinear(resultXPowers);

          if (resultIndex >= 0 && resultIndex < static_cast<int>(result.coeffs_.size())) {
            // Multiply the q-polynomials
            int minQ1 = coeffs_[i].getMaxNegativeIndex();
            int maxQ1 = coeffs_[i].getMaxPositiveIndex();
            int minQ2 = other.coeffs_[j].getMaxNegativeIndex();
            int maxQ2 = other.coeffs_[j].getMaxPositiveIndex();

            for (int q1 = minQ1; q1 <= maxQ1; ++q1) {
              int coeff1 = coeffs_[i][q1];
              if (coeff1 != 0) {
                for (int q2 = minQ2; q2 <= maxQ2; ++q2) {
                  int coeff2 = other.coeffs_[j][q2];
                  if (coeff2 != 0) {
                    result.coeffs_[resultIndex][q1 + q2] += coeff1 * coeff2;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  *this = std::move(result);
  return *this;
}

BMPoly &BMPoly::operator+=(const bilvector<int> &qPoly) {
  // Add qPoly to the x^0 term
  std::vector<int> zeroXPowers(numXVariables, 0);
  expandStorageIfNeeded(zeroXPowers);

  int index = multiIndexToLinear(zeroXPowers);
  if (index >= 0 && index < static_cast<int>(coeffs_.size())) {
    int minQ = std::min(coeffs_[index].getMaxNegativeIndex(), qPoly.getMaxNegativeIndex());
    int maxQ = std::max(coeffs_[index].getMaxPositiveIndex(), qPoly.getMaxPositiveIndex());

    for (int q = minQ; q <= maxQ; ++q) {
      coeffs_[index][q] += qPoly[q];
    }
  }

  return *this;
}

BMPoly &BMPoly::operator-=(const bilvector<int> &qPoly) {
  // Subtract qPoly from the x^0 term
  std::vector<int> zeroXPowers(numXVariables, 0);
  expandStorageIfNeeded(zeroXPowers);

  int index = multiIndexToLinear(zeroXPowers);
  if (index >= 0 && index < static_cast<int>(coeffs_.size())) {
    int minQ = std::min(coeffs_[index].getMaxNegativeIndex(), qPoly.getMaxNegativeIndex());
    int maxQ = std::max(coeffs_[index].getMaxPositiveIndex(), qPoly.getMaxPositiveIndex());

    for (int q = minQ; q <= maxQ; ++q) {
      coeffs_[index][q] -= qPoly[q];
    }
  }

  return *this;
}

BMPoly &BMPoly::operator*=(const bilvector<int> &qPoly) {
  // Multiply each x-monomial's q-polynomial by qPoly
  for (size_t i = 0; i < coeffs_.size(); ++i) {
    if (!coeffs_[i].isZero() && !qPoly.isZero()) {
      bilvector<int> result(0, 1, 20, 0);

      int minQ1 = coeffs_[i].getMaxNegativeIndex();
      int maxQ1 = coeffs_[i].getMaxPositiveIndex();
      int minQ2 = qPoly.getMaxNegativeIndex();
      int maxQ2 = qPoly.getMaxPositiveIndex();

      for (int q1 = minQ1; q1 <= maxQ1; ++q1) {
        int coeff1 = coeffs_[i][q1];
        if (coeff1 != 0) {
          for (int q2 = minQ2; q2 <= maxQ2; ++q2) {
            int coeff2 = qPoly[q2];
            if (coeff2 != 0) {
              result[q1 + q2] += coeff1 * coeff2;
            }
          }
        }
      }

      coeffs_[i] = std::move(result);
    }
  }

  return *this;
}

// Friend function implementations
BMPoly operator+(const BMPoly &lhs, const BMPoly &rhs) {
  BMPoly result = lhs;
  result += rhs;
  return result;
}

BMPoly operator-(const BMPoly &lhs, const BMPoly &rhs) {
  BMPoly result = lhs;
  result -= rhs;
  return result;
}

BMPoly operator*(const BMPoly &lhs, const BMPoly &rhs) {
  BMPoly result = lhs;
  result *= rhs;
  return result;
}

BMPoly operator+(const BMPoly &lhs, const bilvector<int> &rhs) {
  BMPoly result = lhs;
  result += rhs;
  return result;
}

BMPoly operator+(const bilvector<int> &lhs, const BMPoly &rhs) {
  BMPoly result = rhs;
  result += lhs;
  return result;
}

BMPoly operator-(const BMPoly &lhs, const bilvector<int> &rhs) {
  BMPoly result = lhs;
  result -= rhs;
  return result;
}

BMPoly operator-(const bilvector<int> &lhs, const BMPoly &rhs) {
  BMPoly result(rhs.getNumXVariables(), 0, rhs.getMaxXDegrees());
  result += lhs;
  result -= rhs;
  return result;
}

BMPoly operator*(const BMPoly &lhs, const bilvector<int> &rhs) {
  BMPoly result = lhs;
  result *= rhs;
  return result;
}

BMPoly operator*(const bilvector<int> &lhs, const BMPoly &rhs) {
  BMPoly result = rhs;
  result *= lhs;
  return result;
}