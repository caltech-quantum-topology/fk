#include "fk/multivariable_polynomial.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <tuple>

// VectorHash implementation
std::size_t VectorHash::operator()(const std::vector<int> &v) const {
  std::size_t seed = v.size();
  for (auto &i : v) {
    // Mix the hash with signed int support
    seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

// MultivariablePolynomial implementation

MultivariablePolynomial::MultivariablePolynomial(
    int numVariables, int degree, const std::vector<int> &maxDegrees)
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

MultivariablePolynomial::MultivariablePolynomial(
    const MultivariablePolynomial &source, int newNumVariables,
    int targetVariableIndex, int degree, const std::vector<int> &maxDegrees)
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

  // Calculate block sizes for compatibility (not used for indexing)
  blockSizes.resize(newNumVariables);
  if (newNumVariables > 0) {
    blockSizes[0] = 1;
    for (int i = 1; i < newNumVariables; i++) {
      blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
    }
  }

  // Copy coefficients from source, mapping the single variable to
  // targetVariableIndex
  for (const auto &[sourceXPowers, bilvec] : source.coeffs_) {
    // Create new x-powers vector with zeros except at targetVariableIndex
    std::vector<int> newXPowers(newNumVariables, 0);
    newXPowers[targetVariableIndex] = sourceXPowers.e[0];

    // Copy the entire bilvector, using ExponentKey as map key
    coeffs_.emplace(makeKey(newXPowers), bilvec);
  }
}

void MultivariablePolynomial::pruneZeros() {
  auto it = coeffs_.begin();
  while (it != coeffs_.end()) {
    bool isZero = true;
    const auto &bilvec = it->second;

    // Check if this bilvector is all zeros
    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
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

void MultivariablePolynomial::setCoefficient(int qPower,
                                             const std::vector<int> &xPowers,
                                             int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  ExponentKey key = makeKey(xPowers);

  if (coefficient == 0) {
    // If setting to zero, just remove or don't add
    auto it = coeffs_.find(key);
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
    auto it = coeffs_.find(key);
    if (it == coeffs_.end()) {
      // Create new bilvector for this x-monomial
      auto result =
          coeffs_.emplace(makeKey(xPowers), bilvector<int>(0, 1, 20, 0));
      it = result.first;
    }
    it->second[qPower] = coefficient;
  }
}

void MultivariablePolynomial::addToCoefficient(int qPower,
                                               const std::vector<int> &xPowers,
                                               int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (coefficient == 0) {
    return; // Adding zero does nothing
  }
  ExponentKey key = makeKey(xPowers);

  auto it = coeffs_.find(key);
  if (it == coeffs_.end()) {
    // Create new bilvector for this x-monomial
    auto result =
        coeffs_.emplace(makeKey(xPowers), bilvector<int>(0, 1, 20, 0));
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

std::vector<std::pair<std::vector<int>, bilvector<int>>>
MultivariablePolynomial::getCoefficients() const {
  std::vector<std::pair<std::vector<int>, bilvector<int>>> result;
  result.reserve(coeffs_.size());

  for (const auto &entry : coeffs_) {
    const auto &key = entry.first;
    const auto &qPoly = entry.second;

    // Only include non-zero polynomials
    if (!qPoly.isZero()) {
      // Convert ExponentKey -> std::vector<int>
      std::vector<int> xPowers(numXVariables);
      for (int i = 0; i < numXVariables; ++i) {
        xPowers[i] = key.e[i];
      }

      result.emplace_back(std::move(xPowers), qPoly);
    }
  }

  return result;
}

MultivariablePolynomial
MultivariablePolynomial::truncate(const std::vector<int> &maxXdegrees) const {

  MultivariablePolynomial newPoly(this->numXVariables, 0);

  for (const auto &[key, qPoly] : this->coeffs_) {
    bool in_range = true;
    for (int j = 0; j < this->numXVariables; ++j) {
      if (key.e[j] > maxXdegrees[j]) {
        in_range = false;
        break;
      }
    }

    if (in_range) {
      // Copy the term to new polynomial (key is already an ExponentKey)
      newPoly.coeffs_.emplace(key, qPoly);
    }
  }

  return newPoly;
}

MultivariablePolynomial MultivariablePolynomial::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(this->numXVariables, maxDegree);
  return this->truncate(maxXdegrees);
}

void MultivariablePolynomial::clear() { coeffs_.clear(); }

void MultivariablePolynomial::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile;
  outputFile.open(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // Collect and sort x-power keys for deterministic output
  std::vector<ExponentKey> sortedXPowers;
  sortedXPowers.reserve(coeffs_.size());
  for (const auto &[xPowers, bilvec] : coeffs_) {
    sortedXPowers.push_back(xPowers);
  }
  std::sort(sortedXPowers.begin(), sortedXPowers.end());

  bool firstTerm = true;
  for (const auto &xPowers : sortedXPowers) {
    const auto &bilvec = coeffs_.at(xPowers);

    // Collect all non-zero q-terms for this x-power combination
    std::vector<std::pair<int, int>> qTerms; // (q_power, coefficient)
    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      int coeff = bilvec[j];
      if (coeff != 0) {
        qTerms.emplace_back(j, coeff);
      }
    }

    // Only output if there are non-zero terms
    if (!qTerms.empty()) {
      if (!firstTerm) {
        outputFile << ",\n";
      }
      firstTerm = false;

      outputFile << "\t\t{\"x\": [";
      for (size_t k = 0; k < xPowers.e.size(); k++) {
        outputFile << xPowers.e[k];
        if (k < xPowers.e.size() - 1)
          outputFile << ",";
      }
      outputFile << "], \"q_terms\": [";

      for (size_t i = 0; i < qTerms.size(); i++) {
        outputFile << "{\"q\": " << qTerms[i].first
                   << ", \"c\": " << qTerms[i].second << "}";
        if (i < qTerms.size() - 1)
          outputFile << ", ";
      }
      outputFile << "]}";
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
  for (const auto &[key, bilvec] : coeffs_) {
    // Convert ExponentKey -> std::vector<int> for output
    std::vector<int> xPowers(numXVariables);
    for (int i = 0; i < numXVariables; ++i) {
      xPowers[i] = key.e[i];
    }

    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      int coeff = bilvec[j];
      if (coeff != 0) {
        terms.emplace_back(xPowers, j, coeff);
      }
    }
  }
  // Sort terms for consistent output
  std::sort(terms.begin(), terms.end());

  int termCount = 0;
  for (const auto &[xPowers, qPower, coeff] : terms) {
    if (termCount >= maxTerms)
      break;

    if (termCount > 0)
      std::cout << " + ";
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

MultivariablePolynomial &
MultivariablePolynomial::operator+=(const MultivariablePolynomial &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  for (const auto &[key, bilvec] : other.coeffs_) {
    // Convert ExponentKey -> std::vector<int> for addToCoefficient
    std::vector<int> xPowers(numXVariables);
    for (int i = 0; i < numXVariables; ++i) {
      xPowers[i] = key.e[i];
    }

    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      int coeff = bilvec[j];
      if (coeff != 0) {
        addToCoefficient(j, xPowers, coeff);
      }
    }
  }

  return *this;
}

MultivariablePolynomial &
MultivariablePolynomial::operator*=(const MultivariablePolynomial &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  // Create a new coefficient map for the result
  std::unordered_map<ExponentKey, bilvector<int>, ExponentKeyHash> result;

  // Reserve space to avoid repeated rehashing/allocations during multiplication
  result.reserve(coeffs_.size() * other.coeffs_.size());

  // Multiply each term in this polynomial with each term in other polynomial
  for (const auto &[thisKey, thisBilvec] : coeffs_) {
    for (const auto &[otherKey, otherBilvec] : other.coeffs_) {

      // Calculate product x-powers
      std::vector<int> productXPowers(numXVariables);
      ExponentKey productKey;
      for (int i = 0; i < numXVariables; i++) {
        productKey.e[i] = thisKey.e[i] + otherKey.e[i];
      }

      // Multiply all q-coefficient combinations
      for (int thisQ = thisBilvec.getMaxNegativeIndex();
           thisQ <= thisBilvec.getMaxPositiveIndex(); thisQ++) {
        int thisCoeff = thisBilvec[thisQ];
        if (thisCoeff == 0) {
          continue;
        }
        for (int otherQ = otherBilvec.getMaxNegativeIndex();
             otherQ <= otherBilvec.getMaxPositiveIndex(); otherQ++) {
          int otherCoeff = otherBilvec[otherQ];
          if (otherCoeff == 0) {
            continue;
          }
          int productQ = thisQ + otherQ;
          int productCoeff = thisCoeff * otherCoeff;

          // Add to result (lookup + insert in one step)
          auto [it, inserted] = result.try_emplace(
              productKey, 0, 1, 20, 0); // bilvector<int>(0,1,20,0)
          it->second[productQ] += productCoeff;
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

MultivariablePolynomial
MultivariablePolynomial::operator*(const bilvector<int> &qPoly) const {
  MultivariablePolynomial result(numXVariables, 0, maxXDegrees);

  // If this polynomial is zero or qPoly is zero, return zero
  if (coeffs_.empty()) {
    return result;
  }

  // Check if qPoly is zero
  bool qPolyIsZero = true;
  for (int q = qPoly.getMaxNegativeIndex(); q <= qPoly.getMaxPositiveIndex();
       ++q) {
    if (qPoly[q] != 0) {
      qPolyIsZero = false;
      break;
    }
  }
  if (qPolyIsZero) {
    return result;
  }

  // Multiply each term
  for (const auto &[key, thisBilvec] : coeffs_) {
    bilvector<int> product = thisBilvec * qPoly;

    // Add to result
    result.coeffs_.emplace(key, std::move(product));
  }

  return result;
}

MultivariablePolynomial::ExponentKey
MultivariablePolynomial::makeKey(const std::vector<int> &xPowers) const {
  ExponentKey key{}; // all zeros

  // assume xPowers.size() == numXVariables
  for (int i = 0; i < numXVariables; ++i) {
    key.e[i] = xPowers[i];
  }

  // remaining entries stay 0
  return key;
}
