#include "fk/hmpoly.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

// HMPolyVectorHash implementation
std::size_t HMPolyVectorHash::operator()(const std::vector<int> &v) const {
  std::size_t seed = v.size();
  for (auto &i : v) {
    seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

// HMPoly implementation

HMPoly::HMPoly(int numVariables, int degree,
               const std::vector<int> &maxDegrees)
    : numXVariables(numVariables), ctx_initialized(false) {

  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(numVariables, degree);
  } else {
    if (maxDegrees.size() != static_cast<size_t>(numVariables)) {
      throw std::invalid_argument(
          "Max degrees vector size must match number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  // Calculate block sizes for compatibility
  blockSizes.resize(numVariables);
  if (numVariables > 0) {
    blockSizes[0] = 1;
    for (int i = 1; i < numVariables; i++) {
      blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
    }
  }
}

HMPoly::HMPoly(const HMPoly &source, int newNumVariables,
               int targetVariableIndex, int degree,
               const std::vector<int> &maxDegrees)
    : numXVariables(newNumVariables), ctx_initialized(false) {

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

  // Calculate block sizes
  blockSizes.resize(newNumVariables);
  if (newNumVariables > 0) {
    blockSizes[0] = 1;
    for (int i = 1; i < newNumVariables; i++) {
      blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
    }
  }

  // Copy coefficients from source
  for (const auto &[sourceXPowers, bilvec] : source.coeffs_) {
    std::vector<int> newXPowers(newNumVariables, 0);
    newXPowers[targetVariableIndex] = sourceXPowers[0];
    coeffs_.emplace(newXPowers, bilvec);
  }
}

HMPoly::HMPoly(const HMPoly &other)
    : numXVariables(other.numXVariables), maxXDegrees(other.maxXDegrees),
      blockSizes(other.blockSizes), coeffs_(other.coeffs_),
      ctx_initialized(false) {}

HMPoly &HMPoly::operator=(const HMPoly &other) {
  if (this != &other) {
    if (ctx_initialized) {
      clearFlintContext();
    }
    numXVariables = other.numXVariables;
    maxXDegrees = other.maxXDegrees;
    blockSizes = other.blockSizes;
    coeffs_ = other.coeffs_;
    ctx_initialized = false;
  }
  return *this;
}

HMPoly::~HMPoly() {
  if (ctx_initialized) {
    clearFlintContext();
  }
}

void HMPoly::initFlintContext() const {
  if (!ctx_initialized) {
    fmpz_mpoly_ctx_init(ctx, numXVariables + 1, ORD_LEX);
    ctx_initialized = true;
  }
}

void HMPoly::clearFlintContext() const {
  if (ctx_initialized) {
    fmpz_mpoly_ctx_clear(ctx);
    ctx_initialized = false;
  }
}

void HMPoly::toFlint(fmpz_mpoly_t poly) const {
  initFlintContext();
  fmpz_mpoly_init(poly, ctx);

  // Find minimum exponents for ground powers
  std::vector<int> minPowers(numXVariables + 1, 0);
  bool hasTerms = false;

  for (const auto &[xPowers, bilvec] : coeffs_) {
    for (int qPower = bilvec.getMaxNegativeIndex();
         qPower <= bilvec.getMaxPositiveIndex(); qPower++) {
      if (bilvec[qPower] != 0) {
        if (!hasTerms) {
          minPowers[0] = qPower;
          for (int i = 0; i < numXVariables; i++) {
            minPowers[i + 1] = xPowers[i];
          }
          hasTerms = true;
        } else {
          if (qPower < minPowers[0])
            minPowers[0] = qPower;
          for (int i = 0; i < numXVariables; i++) {
            if (xPowers[i] < minPowers[i + 1])
              minPowers[i + 1] = xPowers[i];
          }
        }
      }
    }
  }

  // Convert each term to FLINT
  for (const auto &[xPowers, bilvec] : coeffs_) {
    for (int qPower = bilvec.getMaxNegativeIndex();
         qPower <= bilvec.getMaxPositiveIndex(); qPower++) {
      int coeff = bilvec[qPower];
      if (coeff != 0) {
        // Create exponent array
        fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
        fmpz **exp_ptrs =
            (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));

        for (int i = 0; i <= numXVariables; i++) {
          fmpz_init(&exps[i]);
          exp_ptrs[i] = &exps[i];
        }

        // Set exponents (with ground power offset)
        fmpz_set_si(&exps[0], qPower - minPowers[0]);
        for (int i = 0; i < numXVariables; i++) {
          fmpz_set_si(&exps[i + 1], xPowers[i] - minPowers[i + 1]);
        }

        // Set coefficient
        fmpz_t fmpz_coeff;
        fmpz_init(fmpz_coeff);
        fmpz_set_si(fmpz_coeff, coeff);
        fmpz_mpoly_set_coeff_fmpz_fmpz(poly, fmpz_coeff, exp_ptrs, ctx);

        // Cleanup
        fmpz_clear(fmpz_coeff);
        for (int i = 0; i <= numXVariables; i++) {
          fmpz_clear(&exps[i]);
        }
        flint_free(exp_ptrs);
        flint_free(exps);
      }
    }
  }
}

void HMPoly::fromFlint(const fmpz_mpoly_t poly) {
  coeffs_.clear();
  initFlintContext();

  slong numTerms = fmpz_mpoly_length(poly, ctx);

  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    if (fmpz_is_zero(coeff)) {
      fmpz_clear(coeff);
      continue;
    }

    // Get exponents
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));

    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract powers
    int qPower = fmpz_get_si(&exps[0]);
    std::vector<int> xPowers(numXVariables);
    for (int j = 0; j < numXVariables; j++) {
      xPowers[j] = fmpz_get_si(&exps[j + 1]);
    }

    int coeffValue = fmpz_get_si(coeff);

    // Add to map
    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
      it = result.first;
    }
    it->second[qPower] = coeffValue;

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  pruneZeros();
}

void HMPoly::pruneZeros() {
  auto it = coeffs_.begin();
  while (it != coeffs_.end()) {
    bool isZero = true;
    const auto &bilvec = it->second;

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

void HMPoly::setCoefficient(int qPower, const std::vector<int> &xPowers,
                            int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (coefficient == 0) {
    auto it = coeffs_.find(xPowers);
    if (it != coeffs_.end()) {
      it->second[qPower] = 0;
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
    auto it = coeffs_.find(xPowers);
    if (it == coeffs_.end()) {
      auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
      it = result.first;
    }
    it->second[qPower] = coefficient;
  }
}

void HMPoly::addToCoefficient(int qPower, const std::vector<int> &xPowers,
                              int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (coefficient == 0) {
    return;
  }

  auto it = coeffs_.find(xPowers);
  if (it == coeffs_.end()) {
    auto result = coeffs_.emplace(xPowers, bilvector<int>(0, 1, 20, 0));
    it = result.first;
  }

  it->second[qPower] += coefficient;

  if (it->second[qPower] == 0) {
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
HMPoly::getCoefficients() const {
  std::vector<std::pair<std::vector<int>, bilvector<int>>> result;
  result.reserve(coeffs_.size());

  for (const auto &entry : coeffs_) {
    const auto &xPowers = entry.first;
    const auto &qPoly = entry.second;

    if (!qPoly.isZero()) {
      result.emplace_back(xPowers, qPoly);
    }
  }

  return result;
}

HMPoly HMPoly::truncate(const std::vector<int> &maxXdegrees) const {
  HMPoly newPoly(this->numXVariables, 0);

  for (const auto &[xPowers, qPoly] : this->coeffs_) {
    bool in_range = true;
    for (int j = 0; j < this->numXVariables; ++j) {
      if (xPowers[j] > maxXdegrees[j]) {
        in_range = false;
        break;
      }
    }

    if (in_range) {
      newPoly.coeffs_.emplace(xPowers, qPoly);
    }
  }

  return newPoly;
}

HMPoly HMPoly::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(this->numXVariables, maxDegree);
  return this->truncate(maxXdegrees);
}

void HMPoly::clear() { coeffs_.clear(); }

void HMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile;
  outputFile.open(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  std::vector<std::vector<int>> sortedXPowers;
  sortedXPowers.reserve(coeffs_.size());
  for (const auto &[xPowers, bilvec] : coeffs_) {
    sortedXPowers.push_back(xPowers);
  }
  std::sort(sortedXPowers.begin(), sortedXPowers.end());

  bool firstTerm = true;
  for (const auto &xPowers : sortedXPowers) {
    const auto &bilvec = coeffs_.at(xPowers);

    std::vector<std::pair<int, int>> qTerms;
    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      int coeff = bilvec[j];
      if (coeff != 0) {
        qTerms.emplace_back(j, coeff);
      }
    }

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
  outputFile << "\t\t\"storage_type\": \"hybrid\"\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void HMPoly::print(int maxTerms) const {
  std::cout << "HMPoly P(q";
  for (int i = 0; i < numXVariables; i++) {
    std::cout << ", x" << (i + 1);
  }
  std::cout << "):\n";

  if (coeffs_.empty()) {
    std::cout << "0\n";
    return;
  }

  std::vector<std::tuple<std::vector<int>, int, int>> terms;
  for (const auto &[xPowers, bilvec] : coeffs_) {
    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      int coeff = bilvec[j];
      if (coeff != 0) {
        terms.emplace_back(xPowers, j, coeff);
      }
    }
  }

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

HMPoly &HMPoly::operator+=(const HMPoly &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  // Use map directly - fast O(m) operation
  for (const auto &[xPowers, bilvec] : other.coeffs_) {
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

HMPoly &HMPoly::operator*=(const HMPoly &other) {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  // Use FLINT for multiplication - optimized for bulk operations
  fmpz_mpoly_t thisPoly, otherPoly, resultPoly;

  toFlint(thisPoly);
  other.toFlint(otherPoly);

  initFlintContext();
  fmpz_mpoly_init(resultPoly, ctx);

  // FLINT multiplication
  fmpz_mpoly_mul(resultPoly, thisPoly, otherPoly, ctx);

  // Convert back to map
  fromFlint(resultPoly);

  // Cleanup
  fmpz_mpoly_clear(thisPoly, ctx);
  fmpz_mpoly_clear(otherPoly, ctx);
  fmpz_mpoly_clear(resultPoly, ctx);

  return *this;
}

HMPoly HMPoly::operator*(const bilvector<int> &qPoly) const {
  HMPoly result(numXVariables, 0, maxXDegrees);

  if (coeffs_.empty()) {
    return result;
  }

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

  // Use map directly - element-wise multiplication is fast
  for (const auto &[xPowers, thisBilvec] : coeffs_) {
    bilvector<int> product = thisBilvec * qPoly;
    result.coeffs_.emplace(xPowers, std::move(product));
  }

  return result;
}
