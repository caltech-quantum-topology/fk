#include "fk/zmpoly.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <tuple>

// VectorHashZM implementation
std::size_t VectorHashZM::operator()(const std::vector<int> &v) const {
  std::size_t seed = v.size();
  for (auto &i : v) {
    // Mix the hash with signed int support
    seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

// ZMPoly implementation

ZMPoly::ZMPoly(int numVariables, int degree,
               const std::vector<int> &maxDegrees)
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

ZMPoly::ZMPoly(const ZMPoly &source, int newNumVariables,
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
  for (const auto &[sourceXPowers, qPoly] : source.coeffs_) {
    // Create new x-powers vector with zeros except at targetVariableIndex
    std::vector<int> newXPowers(newNumVariables, 0);
    newXPowers[targetVariableIndex] = sourceXPowers[0];

    // Copy the QPolynomial
    coeffs_.emplace(newXPowers, qPoly);
  }
}

void ZMPoly::pruneZeros() {
  auto it = coeffs_.begin();
  while (it != coeffs_.end()) {
    if (it->second.isZero()) {
      it = coeffs_.erase(it);
    } else {
      ++it;
    }
  }
}

void ZMPoly::setCoefficient(int qPower, const std::vector<int> &xPowers,
                            int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (coefficient == 0) {
    // If setting to zero, just remove or don't add
    auto it = coeffs_.find(xPowers);
    if (it != coeffs_.end()) {
      it->second.setCoefficient(qPower, 0);
      // Check if entire QPolynomial became zero and remove if so
      if (it->second.isZero()) {
        coeffs_.erase(it);
      }
    }
  } else {
    // Create entry if it doesn't exist, then set coefficient
    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      // Create new QPolynomial for this x-monomial
      auto result = coeffs_.emplace(xPowers, QPolynomial());
      it = result.first;
    }
    it->second.setCoefficient(qPower, coefficient);
  }
}

void ZMPoly::addToCoefficient(int qPower, const std::vector<int> &xPowers,
                              int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (coefficient == 0) {
    return; // Adding zero does nothing
  }

  auto it = coeffs_.find(xPowers);
  if (it == coeffs_.end()) {
    // Create new QPolynomial for this x-monomial
    auto result = coeffs_.emplace(xPowers, QPolynomial());
    it = result.first;
  }

  it->second.addToCoefficient(qPower, coefficient);

  // If the result became zero, clean up
  if (it->second.isZero()) {
    coeffs_.erase(it);
  }
}

std::vector<std::pair<std::vector<int>, QPolynomial>>
ZMPoly::getCoefficients() const {
  std::vector<std::pair<std::vector<int>, QPolynomial>> result;
  result.reserve(coeffs_.size());

  for (const auto &entry : coeffs_) {
    const auto &xPowers = entry.first;
    const auto &qPoly = entry.second;

    // Only include non-zero polynomials
    if (!qPoly.isZero()) {
      result.emplace_back(xPowers, qPoly);
    }
  }

  return result;
}

ZMPoly ZMPoly::truncate(const std::vector<int> &maxXdegrees) const {

  ZMPoly newPoly(this->numXVariables, 0);
  for (const auto &[xPowers, qPoly] : this->coeffs_) {
    bool in_range = true;
    for (int j = 0; j < this->numXVariables; ++j) {
      if (xPowers[j] > maxXdegrees[j]) {
        in_range = false;
        break;
      }
    }

    if (in_range) {
      // Copy the term to new polynomial
      newPoly.coeffs_.emplace(xPowers, qPoly);
    }
  }

  return newPoly;
}

ZMPoly ZMPoly::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(this->numXVariables, maxDegree);
  return this->truncate(maxXdegrees);
}

void ZMPoly::clear() { coeffs_.clear(); }

void ZMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile;
  outputFile.open(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // Collect and sort x-power keys for deterministic output
  std::vector<std::vector<int>> sortedXPowers;
  sortedXPowers.reserve(coeffs_.size());
  for (const auto &[xPowers, qPoly] : coeffs_) {
    sortedXPowers.push_back(xPowers);
  }
  std::sort(sortedXPowers.begin(), sortedXPowers.end());

  bool firstTerm = true;
  for (const auto &xPowers : sortedXPowers) {
    const auto &qPoly = coeffs_.at(xPowers);

    // Collect all non-zero q-terms for this x-power combination
    std::vector<std::pair<int, int>> qTerms; // (q_power, coefficient)
    const int minQ = qPoly.getMinPower();
    const int maxQ = qPoly.getMaxPower();
    for (int j = minQ; j <= maxQ; j++) {
      int coeff = qPoly.getCoefficient(j);
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
      for (size_t k = 0; k < xPowers.size(); k++) {
        outputFile << xPowers[k];
        if (k < xPowers.size() - 1)
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
  outputFile << "\t\t\"storage_type\": \"sparse\",\n";
  outputFile << "\t\t\"coefficient_type\": \"arbitrary_precision\"\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void ZMPoly::print(int maxTerms) const {
  std::cout << "ZMPoly (Arbitrary Precision) P(q";
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
  for (const auto &[xPowers, qPoly] : coeffs_) {
    const int minQ = qPoly.getMinPower();
    const int maxQ = qPoly.getMaxPower();
    for (int j = minQ; j <= maxQ; j++) {
      int coeff = qPoly.getCoefficient(j);
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

ZMPoly &ZMPoly::operator+=(const ZMPoly &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  // Reserve space in coeffs_ if needed to reduce rehashing
  if (coeffs_.size() + other.coeffs_.size() > coeffs_.bucket_count()) {
    coeffs_.reserve(coeffs_.size() + other.coeffs_.size());
  }

  for (const auto &[xPowers, qPoly] : other.coeffs_) {
    // Cache range to avoid repeated function calls
    const int minQ = qPoly.getMinPower();
    const int maxQ = qPoly.getMaxPower();

    for (int j = minQ; j <= maxQ; j++) {
      int coeff = qPoly.getCoefficient(j);
      if (coeff != 0) {
        addToCoefficient(j, xPowers, coeff);
      }
    }
  }

  return *this;
}

ZMPoly &ZMPoly::operator*=(const ZMPoly &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  // Create a new coefficient map for the result
  // Reserve capacity to reduce rehashing (estimated size)
  std::unordered_map<std::vector<int>, QPolynomial, VectorHashZM> result;
  result.reserve(coeffs_.size() * other.coeffs_.size());

  // Multiply each term in this polynomial with each term in other polynomial
  for (const auto &[thisXPowers, thisQPoly] : coeffs_) {
    // Cache range for thisQPoly to avoid repeated function calls
    const int thisMinQ = thisQPoly.getMinPower();
    const int thisMaxQ = thisQPoly.getMaxPower();

    for (const auto &[otherXPowers, otherQPoly] : other.coeffs_) {
      // Cache range for otherQPoly to avoid repeated function calls
      const int otherMinQ = otherQPoly.getMinPower();
      const int otherMaxQ = otherQPoly.getMaxPower();

      // Calculate product x-powers
      std::vector<int> productXPowers(numXVariables);
      for (int i = 0; i < numXVariables; i++) {
        productXPowers[i] = thisXPowers[i] + otherXPowers[i];
      }

      // Find or create the result entry once for this x-power combination
      auto it = result.find(productXPowers);
      if (it == result.end()) {
        auto insertResult = result.emplace(productXPowers, QPolynomial());
        it = insertResult.first;
      }

      // Multiply all q-coefficient combinations and add to the same result entry
      for (int thisQ = thisMinQ; thisQ <= thisMaxQ; thisQ++) {
        int thisCoeff = thisQPoly.getCoefficient(thisQ);
        if (thisCoeff != 0) {
          for (int otherQ = otherMinQ; otherQ <= otherMaxQ; otherQ++) {
            int otherCoeff = otherQPoly.getCoefficient(otherQ);
            if (otherCoeff != 0) {
              int productQ = thisQ + otherQ;
              int productCoeff = thisCoeff * otherCoeff;
              it->second.addToCoefficient(productQ, productCoeff);
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

ZMPoly ZMPoly::operator*(const QPolynomial &qPoly) const {
  ZMPoly result(numXVariables, 0, maxXDegrees);

  // If this polynomial is zero or qPoly is zero, return zero
  if (coeffs_.empty() || qPoly.isZero()) {
    return result;
  }

  // Multiply each term
  for (const auto &[xPowers, thisQPoly] : coeffs_) {
    QPolynomial product = thisQPoly * qPoly;

    // Add to result (only if non-zero)
    if (!product.isZero()) {
      result.coeffs_.emplace(xPowers, std::move(product));
    }
  }

  return result;
}
