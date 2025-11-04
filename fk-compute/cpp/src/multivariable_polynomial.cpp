#include "fk/multivariable_polynomial.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <tuple>

// VectorHash implementation
std::size_t VectorHash::operator()(const std::vector<int>& v) const {
  std::size_t seed = v.size();
  for (auto& i : v) {
    // Mix the hash with signed int support
    seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

// MultivariablePolynomial implementation

MultivariablePolynomial::MultivariablePolynomial(int numVariables, int degree,
                        const std::vector<int>&maxDegrees)
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

void MultivariablePolynomial::pruneZeros() {
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

int MultivariablePolynomial::getCoefficient(int qPower, const std::vector<int> &xPowers) const {
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

void MultivariablePolynomial::setCoefficient(int qPower, const std::vector<int> &xPowers,
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

void MultivariablePolynomial::addToCoefficient(int qPower, const std::vector<int> &xPowers,
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

bilvector<int> &MultivariablePolynomial::getQPolynomial(const std::vector<int> &xPowers) {
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

const bilvector<int> &MultivariablePolynomial::getQPolynomial(const std::vector<int> &xPowers) const {
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

const std::unordered_map<std::vector<int>, bilvector<int>, VectorHash>&
MultivariablePolynomial::getCoefficientMap() const {
  return coeffs_;
}

std::vector<bilvector<int>> MultivariablePolynomial::getCoefficients() const {
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

std::vector<bilvector<int>>& MultivariablePolynomial::getCoefficients() {
  // Create a thread-local copy for modification
  static thread_local std::vector<bilvector<int>> temp_result;
  temp_result = static_cast<const MultivariablePolynomial*>(this)->getCoefficients();
  return temp_result;
}

void MultivariablePolynomial::syncFromDenseVector(const std::vector<bilvector<int>>& denseVector) {
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

int MultivariablePolynomial::getNumXVariables() const {
  return numXVariables;
}

const std::vector<int> &MultivariablePolynomial::getMaxXDegrees() const {
  return maxXDegrees;
}

const std::vector<int> &MultivariablePolynomial::getBlockSizes() const {
  return blockSizes;
}

void MultivariablePolynomial::clear() {
  coeffs_.clear();
}

bool MultivariablePolynomial::isZero() const {
  return coeffs_.empty();
}

void MultivariablePolynomial::exportToJson(const std::string &fileName) const {
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

void MultivariablePolynomial::print(int maxTerms) const {
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

void MultivariablePolynomial::checkCompatibility(const MultivariablePolynomial &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }
}

MultivariablePolynomial &MultivariablePolynomial::operator+=(const MultivariablePolynomial &other) {
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

MultivariablePolynomial &MultivariablePolynomial::operator-=(const MultivariablePolynomial &other) {
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

MultivariablePolynomial &MultivariablePolynomial::operator*=(const MultivariablePolynomial &other) {
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

bilvector<int> MultivariablePolynomial::evaluate(const std::vector<int>& point) const {
  if (point.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("Point dimension must match number of variables");
  }

  bilvector<int> result(0, 0, 1, 0);

  for (const auto& [xPowers, qPoly] : coeffs_) {
    int termValue = 1;

    for (size_t i = 0; i < xPowers.size(); ++i) {
      if (xPowers[i] > 0) {
        for (int exp = 0; exp < xPowers[i]; ++exp) {
          termValue *= point[i];
        }
      } else if (xPowers[i] < 0) {
        if (point[i] == 0) {
          throw std::domain_error("Division by zero: negative exponent with zero value at variable " + std::to_string(i));
        }
        for (int exp = 0; exp > xPowers[i]; --exp) {
          if (termValue % point[i] != 0) {
            throw std::domain_error("Division would result in non-integer coefficient");
          }
          termValue /= point[i];
        }
      }
    }

    for (int qPower = qPoly.getMaxNegativeIndex(); qPower <= qPoly.getMaxPositiveIndex(); ++qPower) {
      int coefficient = qPoly[qPower];
      if (coefficient != 0) {
        result[qPower] += coefficient * termValue;
      }
    }
  }

  return result;
}

// Friend functions for binary operations
MultivariablePolynomial operator+(const MultivariablePolynomial &lhs,
                                 const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result += rhs;
  return result;
}

MultivariablePolynomial operator-(const MultivariablePolynomial &lhs,
                                 const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result -= rhs;
  return result;
}

MultivariablePolynomial operator*(const MultivariablePolynomial &lhs,
                                 const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result *= rhs;
  return result;
}

MultivariablePolynomial &
MultivariablePolynomial::operator+=(const bilvector<int> &qPoly) {
  // Monomial x_1^0 x_2^0 ... x_n^0
  std::vector<int> zeroXPowers(numXVariables, 0);

  for (int q = qPoly.getMaxNegativeIndex(); q <= qPoly.getMaxPositiveIndex(); ++q) {
    int coeff = qPoly[q];
    if (coeff != 0) {
      addToCoefficient(q, zeroXPowers, coeff);
    }
  }

  return *this;
}

MultivariablePolynomial &
MultivariablePolynomial::operator-=(const bilvector<int> &qPoly) {
  std::vector<int> zeroXPowers(numXVariables, 0);

  for (int q = qPoly.getMaxNegativeIndex(); q <= qPoly.getMaxPositiveIndex(); ++q) {
    int coeff = qPoly[q];
    if (coeff != 0) {
      addToCoefficient(q, zeroXPowers, -coeff);
    }
  }

  return *this;
}

MultivariablePolynomial &
MultivariablePolynomial::operator*=(const bilvector<int> &qPoly) {
  // If this polynomial is already zero, nothing to do
  if (coeffs_.empty()) {
    return *this;
  }

  // Check if qPoly is the zero polynomial
  bool qPolyIsZero = true;
  for (int q = qPoly.getMaxNegativeIndex(); q <= qPoly.getMaxPositiveIndex(); ++q) {
    if (qPoly[q] != 0) {
      qPolyIsZero = false;
      break;
    }
  }

  if (qPolyIsZero) {
    clear();  // P * 0 = 0
    return *this;
  }

  std::unordered_map<std::vector<int>, bilvector<int>, VectorHash> result;

  int factorMin = qPoly.getMaxNegativeIndex();
  int factorMax = qPoly.getMaxPositiveIndex();

  for (const auto &entry : coeffs_) {
    const std::vector<int> &xPowers = entry.first;
    const bilvector<int> &thisBilvec = entry.second;

    int thisMin = thisBilvec.getMaxNegativeIndex();
    int thisMax = thisBilvec.getMaxPositiveIndex();

    // Start with a minimal bilvector; it grows as needed
    bilvector<int> product(0, 1, 20, 0);

    for (int qi = thisMin; qi <= thisMax; ++qi) {
      int a = thisBilvec[qi];
      if (a == 0) continue;

      for (int qj = factorMin; qj <= factorMax; ++qj) {
        int b = qPoly[qj];
        if (b == 0) continue;

        product[qi + qj] += a * b;
      }
    }

    // Store only if product is not the zero polynomial
    bool allZero = true;
    for (int q = product.getMaxNegativeIndex(); q <= product.getMaxPositiveIndex(); ++q) {
      if (product[q] != 0) {
        allZero = false;
        break;
      }
    }

    if (!allZero) {
      result.emplace(xPowers, std::move(product));
    }
  }

  coeffs_ = std::move(result);
  return *this;
}

MultivariablePolynomial operator+(const MultivariablePolynomial &lhs,
                                 const bilvector<int> &rhs) {
  MultivariablePolynomial result = lhs;
  result += rhs;
  return result;
}

MultivariablePolynomial operator+(const bilvector<int> &lhs,
                                 const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = rhs;
  result += lhs;
  return result;
}

MultivariablePolynomial operator-(const MultivariablePolynomial &lhs,
                                 const bilvector<int> &rhs) {
  MultivariablePolynomial result = lhs;
  result -= rhs;
  return result;
}

MultivariablePolynomial operator-(const bilvector<int> &lhs,
                                 const MultivariablePolynomial &rhs) {
  // Construct a polynomial with the same x-structure as rhs, then do lhs - rhs
  MultivariablePolynomial result(rhs.getNumXVariables(),
                                 /*degree (unused when maxDegrees provided)*/ 0,
                                 rhs.getMaxXDegrees());

  result += lhs;   // attach lhs(q) to x^0...0
  result -= rhs;   // subtract rhs
  return result;
}

MultivariablePolynomial operator*(const MultivariablePolynomial &lhs,
                                 const bilvector<int> &rhs) {
  MultivariablePolynomial result = lhs;
  result *= rhs;
  return result;
}

MultivariablePolynomial operator*(const bilvector<int> &lhs,
                                 const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = rhs;
  result *= lhs;
  return result;
}
