#pragma once

#include "fk/bilvector.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

/**
 * Custom hash function for vector<int> keys
 */
struct VectorHash {
  std::size_t operator()(const std::vector<int>& v) const {
    std::size_t seed = v.size();
    for (auto& i : v) {
      // Mix the hash with signed int support
      seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

/**
 * MultivariablePolynomial: A polynomial in variables q, x₁, x₂, ..., xₙ
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - Arbitrary positive/negative q powers
 * - Configurable number of x variables
 * - Sparse storage using unordered_map (only non-zero coefficients)
 * - Support for negative x exponents
 * - Efficient indexing and arithmetic operations
 * - Header-only implementation
 */
class MultivariablePolynomial {
private:
  int numXVariables;            // Number of x variables (components)
  std::vector<int> maxXDegrees; // Max degrees (advisory only, not enforced)
  std::vector<int> blockSizes;  // Block sizes (advisory only, for compatibility)

  // Sparse storage: map from x-exponent vector to q-polynomial
  std::unordered_map<std::vector<int>, bilvector<int>, VectorHash> coeffs_;

  // Prune zero coefficients from the map
  void pruneZeros() {
    auto it = coeffs_.begin();
    while (it != coeffs_.end()) {
      bool isZero = true;
      const auto& bilvec = it->second;

      // Check if this bilvector is all zeros
      for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
        if (bilvec[j] != 0) {
          isZero = false;
          break;
        }
      }

      if (isZero) {
        it = coeffs_.erase(it);
      } else {
        ++it;
      }
    }
  }

public:
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable (advisory only)
   */
  MultivariablePolynomial(int numVariables, int degree = 10,
                          const std::vector<int> &maxDegrees = {})
      : numXVariables(numVariables) {

    if (maxDegrees.empty()) {
      maxXDegrees = std::vector<int>(numVariables, degree);
    } else {
      if (maxDegrees.size() != static_cast<size_t>(numVariables)) {
        throw std::invalid_argument(
            "Max degrees vector size must match number of variables");
      }
      maxXDegrees = maxDegrees;
    }

    // Calculate block sizes for compatibility (not used for indexing)
    blockSizes.resize(numVariables);
    if (numVariables > 0) {
      blockSizes[0] = 1;
      for (int i = 1; i < numVariables; i++) {
        blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
      }
    }
  }

  /**
   * Get coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @return Coefficient value
   */
  int getCoefficient(int qPower, const std::vector<int> &xPowers) const {
    if (xPowers.size() != static_cast<size_t>(numXVariables)) {
      throw std::invalid_argument(
          "X powers vector size must match number of variables");
    }

    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      return 0;  // Not found means coefficient is 0
    }
    return it->second[qPower];
  }

  /**
   * Set coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient New coefficient value
   */
  void setCoefficient(int qPower, const std::vector<int> &xPowers,
                      int coefficient) {
    if (xPowers.size() != static_cast<size_t>(numXVariables)) {
      throw std::invalid_argument(
          "X powers vector size must match number of variables");
    }

    if (coefficient == 0) {
      // If setting to zero, just remove or don't add
      auto it = coeffs_.find(xPowers);
      if (it != coeffs_.end()) {
        it->second[qPower] = 0;
        // Check if entire bilvector became zero and remove if so
        bool allZero = true;
        for (int j = it->second.getMaxNegativeIndex();
             j <= it->second.getMaxPositiveIndex(); j++) {
          if (it->second[j] != 0) {
            allZero = false;
            break;
          }
        }
        if (allZero) {
          coeffs_.erase(it);
        }
      }
    } else {
      // Create entry if it doesn't exist, then set coefficient
      auto it = coeffs_.find(xPowers);
      if (it == coeffs_.end()) {
        // Create new bilvector for this x-monomial
        auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
        it = result.first;
      }
      it->second[qPower] = coefficient;
    }
  }

  /**
   * Add to coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient Value to add
   */
  void addToCoefficient(int qPower, const std::vector<int> &xPowers,
                        int coefficient) {
    if (xPowers.size() != static_cast<size_t>(numXVariables)) {
      throw std::invalid_argument(
          "X powers vector size must match number of variables");
    }

    if (coefficient == 0) {
      return;  // Adding zero does nothing
    }

    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      // Create new bilvector for this x-monomial
      auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
      it = result.first;
    }

    it->second[qPower] += coefficient;

    // If the result became zero, we might want to clean up
    if (it->second[qPower] == 0) {
      // Check if entire bilvector became zero
      bool allZero = true;
      for (int j = it->second.getMaxNegativeIndex();
           j <= it->second.getMaxPositiveIndex(); j++) {
        if (it->second[j] != 0) {
          allZero = false;
          break;
        }
      }
      if (allZero) {
        coeffs_.erase(it);
      }
    }
  }

  /**
   * Get access to the bilvector for a specific x-multi-index
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @return Reference to bilvector containing q-coefficients
   */
  bilvector<int> &getQPolynomial(const std::vector<int> &xPowers) {
    if (xPowers.size() != static_cast<size_t>(numXVariables)) {
      throw std::invalid_argument(
          "X powers vector size must match number of variables");
    }

    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      // Create new bilvector for this x-monomial
      auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
      it = result.first;
    }
    return it->second;
  }

  const bilvector<int> &getQPolynomial(const std::vector<int> &xPowers) const {
    if (xPowers.size() != static_cast<size_t>(numXVariables)) {
      throw std::invalid_argument(
          "X powers vector size must match number of variables");
    }

    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      // Return a static zero bilvector
      static bilvector<int> zeroBilvector(0, 1, 20, 0);
      return zeroBilvector;
    }
    return it->second;
  }

  /**
   * Get read-only access to coefficient map (replaces getCoefficients)
   */
  const std::unordered_map<std::vector<int>, bilvector<int>, VectorHash>&
  getCoefficientMap() const {
    return coeffs_;
  }

  /**
   * Backward compatibility: Convert sparse representation to dense vector
   * This is needed for existing code that expects vector<bilvector<int>>
   */
  std::vector<bilvector<int>> getCoefficients() const {
    // Calculate total size needed for dense representation
    int totalSize = 1;
    for (int i = 0; i < numXVariables; i++) {
      totalSize *= (maxXDegrees[i] + 1);
    }

    // Create dense vector initialized with zero bilvectors
    std::vector<bilvector<int>> result(totalSize, bilvector<int>(0, 1, 20, 0));

    // Fill in the non-zero entries from sparse map
    for (const auto& [xPowers, bilvec] : coeffs_) {
      // Convert x-powers to linear index using old formula
      int linearIndex = 0;
      bool withinBounds = true;

      for (int i = 0; i < numXVariables; i++) {
        if (xPowers[i] < 0 || xPowers[i] > maxXDegrees[i]) {
          withinBounds = false;
          break;
        }
        linearIndex += xPowers[i] * blockSizes[i];
      }

      if (withinBounds && linearIndex >= 0 && linearIndex < totalSize) {
        // Copy the bilvector - note: this is expensive but needed for compatibility
        for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
          if (bilvec[j] != 0) {
            result[linearIndex][j] = bilvec[j];
          }
        }
      }
    }

    return result;
  }

  /**
   * Backward compatibility: Non-const version
   * WARNING: This creates a mutable copy that must be manually synced back
   */
  std::vector<bilvector<int>>& getCoefficients() {
    // Create a thread-local copy for modification
    static thread_local std::vector<bilvector<int>> temp_result;
    temp_result = static_cast<const MultivariablePolynomial*>(this)->getCoefficients();
    return temp_result;
  }

  /**
   * Sync dense vector back to sparse representation
   * This must be called after any operations that modify the dense vector
   */
  void syncFromDenseVector(const std::vector<bilvector<int>>& denseVector) {
    // Clear current sparse representation
    coeffs_.clear();

    // Convert dense back to sparse
    for (size_t linearIndex = 0; linearIndex < denseVector.size(); linearIndex++) {
      // Convert linear index back to x-powers
      std::vector<int> xPowers(numXVariables);
      size_t remainingIndex = linearIndex;

      for (int i = numXVariables - 1; i >= 0; i--) {
        if (i < static_cast<int>(blockSizes.size()) && blockSizes[i] > 0) {
          xPowers[i] = static_cast<int>(remainingIndex / blockSizes[i]) % (maxXDegrees[i] + 1);
          remainingIndex %= blockSizes[i];
        }
      }

      // Check if this bilvector has any non-zero coefficients
      bool hasNonZero = false;
      for (int j = denseVector[linearIndex].getMaxNegativeIndex();
           j <= denseVector[linearIndex].getMaxPositiveIndex(); j++) {
        if (denseVector[linearIndex][j] != 0) {
          hasNonZero = true;
          break;
        }
      }

      if (hasNonZero) {
        // Copy the entire bilvector using emplace
        coeffs_.emplace(xPowers, denseVector[linearIndex]);
      }
    }
  }

  /**
   * Get number of x variables
   */
  int getNumXVariables() const {
    return numXVariables;
  }

  /**
   * Get maximum degrees for each x variable (advisory only)
   */
  const std::vector<int> &getMaxXDegrees() const {
    return maxXDegrees;
  }

  /**
   * Get block sizes (for compatibility, not used internally)
   */
  const std::vector<int> &getBlockSizes() const {
    return blockSizes;
  }

  /**
   * Clear all coefficients
   */
  void clear() {
    coeffs_.clear();
  }

  /**
   * Check if polynomial is zero
   */
  bool isZero() const {
    return coeffs_.empty();
  }

  /**
   * Export to JSON format with new sparse format
   * @param fileName Output file name
   */
  void exportToJson(const std::string &fileName) const {
    std::ofstream outputFile;
    outputFile.open(fileName + ".json");
    outputFile << "{\n\t\"terms\":[\n";

    bool firstTerm = true;
    for (const auto& [xPowers, bilvec] : coeffs_) {
      for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
        int coeff = bilvec[j];
        if (coeff != 0) {
          if (!firstTerm) {
            outputFile << ",\n";
          }
          firstTerm = false;

          outputFile << "\t\t{\"x\": [";
          for (size_t k = 0; k < xPowers.size(); k++) {
            outputFile << xPowers[k];
            if (k < xPowers.size() - 1) outputFile << ",";
          }
          outputFile << "], \"q\": " << j << ", \"c\": " << coeff << "}";
        }
      }
    }

    outputFile << "\n\t],\n";
    outputFile << "\t\"metadata\": {\n";
    outputFile << "\t\t\"num_x_variables\": " << numXVariables << ",\n";
    outputFile << "\t\t\"max_x_degrees\": [";
    for (int i = 0; i < numXVariables; i++) {
      outputFile << maxXDegrees[i];
      if (i < numXVariables - 1)
        outputFile << ",";
    }
    outputFile << "],\n";
    outputFile << "\t\t\"storage_type\": \"sparse\"\n";
    outputFile << "\t}\n}";
    outputFile.close();
  }

  /**
   * Print polynomial in human-readable format (for debugging)
   * @param maxTerms Maximum number of terms to print
   */
  void print(int maxTerms = 10) const {
    std::cout << "Multivariable Polynomial P(q";
    for (int i = 0; i < numXVariables; i++) {
      std::cout << ", x" << (i + 1);
    }
    std::cout << "):\n";

    if (coeffs_.empty()) {
      std::cout << "0\n";
      return;
    }

    // Collect all terms and sort for deterministic output
    std::vector<std::tuple<std::vector<int>, int, int>> terms;
    for (const auto& [xPowers, bilvec] : coeffs_) {
      for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
        int coeff = bilvec[j];
        if (coeff != 0) {
          terms.emplace_back(xPowers, j, coeff);
        }
      }
    }

    // Sort terms for consistent output
    std::sort(terms.begin(), terms.end());

    int termCount = 0;
    for (const auto& [xPowers, qPower, coeff] : terms) {
      if (termCount >= maxTerms) break;

      if (termCount > 0) std::cout << " + ";
      std::cout << coeff << "*q^" << qPower;
      for (size_t k = 0; k < xPowers.size(); k++) {
        if (xPowers[k] != 0) {
          std::cout << "*x" << (k + 1);
          if (xPowers[k] != 1) {
            std::cout << "^" << xPowers[k];
          }
        }
      }
      termCount++;
    }

    if (termCount == maxTerms && terms.size() > static_cast<size_t>(maxTerms)) {
      std::cout << " + ...";
    }
    std::cout << std::endl;
  }

  /**
   * Check if two polynomials are compatible for arithmetic operations
   */
  void checkCompatibility(const MultivariablePolynomial &other) const {
    if (numXVariables != other.numXVariables) {
      throw std::invalid_argument(
          "Polynomials must have the same number of x variables");
    }
  }

  /**
   * Add another polynomial to this one
   */
  MultivariablePolynomial &operator+=(const MultivariablePolynomial &other) {
    checkCompatibility(other);

    for (const auto& [xPowers, bilvec] : other.coeffs_) {
      for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
        int coeff = bilvec[j];
        if (coeff != 0) {
          addToCoefficient(j, xPowers, coeff);
        }
      }
    }

    return *this;
  }

  /**
   * Subtract another polynomial from this one
   */
  MultivariablePolynomial &operator-=(const MultivariablePolynomial &other) {
    checkCompatibility(other);

    for (const auto& [xPowers, bilvec] : other.coeffs_) {
      for (int j = bilvec.getMaxNegativeIndex(); j <= bilvec.getMaxPositiveIndex(); j++) {
        int coeff = bilvec[j];
        if (coeff != 0) {
          addToCoefficient(j, xPowers, -coeff);
        }
      }
    }

    return *this;
  }

  /**
   * Multiply this polynomial by another
   */
  MultivariablePolynomial &operator*=(const MultivariablePolynomial &other) {
    checkCompatibility(other);

    // Create a new coefficient map for the result
    std::unordered_map<std::vector<int>, bilvector<int>, VectorHash> result;

    // Multiply each term in this polynomial with each term in other polynomial
    for (const auto& [thisXPowers, thisBilvec] : coeffs_) {
      for (const auto& [otherXPowers, otherBilvec] : other.coeffs_) {

        // Calculate product x-powers
        std::vector<int> productXPowers(numXVariables);
        for (int i = 0; i < numXVariables; i++) {
          productXPowers[i] = thisXPowers[i] + otherXPowers[i];
        }

        // Multiply all q-coefficient combinations
        for (int thisQ = thisBilvec.getMaxNegativeIndex();
             thisQ <= thisBilvec.getMaxPositiveIndex(); thisQ++) {
          int thisCoeff = thisBilvec[thisQ];
          if (thisCoeff != 0) {
            for (int otherQ = otherBilvec.getMaxNegativeIndex();
                 otherQ <= otherBilvec.getMaxPositiveIndex(); otherQ++) {
              int otherCoeff = otherBilvec[otherQ];
              if (otherCoeff != 0) {
                int productQ = thisQ + otherQ;
                int productCoeff = thisCoeff * otherCoeff;

                // Add to result
                auto it = result.find(productXPowers);
                if (it == result.end()) {
                  auto insertResult = result.emplace(productXPowers, bilvector<int>(0, 1, 20, 0));
                  it = insertResult.first;
                }
                it->second[productQ] += productCoeff;
              }
            }
          }
        }
      }
    }

    // Replace coefficients with result
    coeffs_ = std::move(result);

    // Clean up any zeros that might have been created
    pruneZeros();

    return *this;
  }

  /**
   * Friend functions for binary operations
   */
  friend MultivariablePolynomial operator+(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs) {
    MultivariablePolynomial result = lhs;
    result += rhs;
    return result;
  }

  friend MultivariablePolynomial operator-(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs) {
    MultivariablePolynomial result = lhs;
    result -= rhs;
    return result;
  }

  friend MultivariablePolynomial operator*(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs) {
    MultivariablePolynomial result = lhs;
    result *= rhs;
    return result;
  }
};