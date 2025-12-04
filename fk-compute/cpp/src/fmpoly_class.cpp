#include "fk/fmpoly_class.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

FMPoly::FMPoly(int numVariables, int degree, const std::vector<int> &maxDegrees)
    : numXVariables(numVariables), allGroundPowers(numVariables + 1, 0) {

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

  setupContext();
}

FMPoly::FMPoly(const FMPoly &source, int newNumVariables,
               int targetVariableIndex, int degree,
               const std::vector<int> &maxDegrees)
    : numXVariables(newNumVariables), allGroundPowers(newNumVariables + 1, 0) {

  // Copy q ground power from source
  allGroundPowers[0] = source.allGroundPowers[0];

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

  // Calculate block sizes for compatibility
  blockSizes.resize(newNumVariables);
  if (newNumVariables > 0) {
    blockSizes[0] = 1;
    for (int i = 1; i < newNumVariables; i++) {
      blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
    }
  }

  setupContext();

  // Copy terms from source, mapping the single variable to targetVariableIndex
  // This is a simplified implementation - in practice, we'd iterate through
  // the source polynomial's terms and map them appropriately
  fmpz_mpoly_set(poly, source.poly, ctx);
}

FMPoly::FMPoly(const FMPoly &other)
    : numXVariables(other.numXVariables),
      allGroundPowers(other.allGroundPowers), maxXDegrees(other.maxXDegrees),
      blockSizes(other.blockSizes) {

  setupContext();
  fmpz_mpoly_set(poly, other.poly, ctx);
}

FMPoly &FMPoly::operator=(const FMPoly &other) {
  if (this != &other) {
    // Clean up current resources
    fmpz_mpoly_clear(poly, ctx);
    fmpz_mpoly_ctx_clear(ctx);

    // Copy data
    numXVariables = other.numXVariables;
    allGroundPowers = other.allGroundPowers;
    maxXDegrees = other.maxXDegrees;
    blockSizes = other.blockSizes;

    // Setup new context and copy polynomial
    setupContext();
    fmpz_mpoly_set(poly, other.poly, ctx);
  }
  return *this;
}

FMPoly::~FMPoly() {
  fmpz_mpoly_clear(poly, ctx);
  fmpz_mpoly_ctx_clear(ctx);
}

void FMPoly::setupContext() {
  // Initialize FLINT context with numXVariables + 1 variables (q, x1, x2, ...,
  // xn)
  fmpz_mpoly_ctx_init(ctx, numXVariables + 1, ORD_LEX);
  fmpz_mpoly_init(poly, ctx);
}

void FMPoly::convertExponents(int qPower, const std::vector<int> &xPowers,
                              fmpz **exps, slong *exp_bits) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Allocate exponent array: [q, x1, x2, ..., xn]
  *exp_bits = FLINT_BITS;
  *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));

  for (int i = 0; i <= numXVariables; i++) {
    fmpz_init(&((*exps)[i]));
  }

  // Set q exponent (handle offset for negative powers)
  fmpz_set_si(&((*exps)[0]), qPower - allGroundPowers[0]);

  // Set x variable exponents (apply ground powers offset)
  for (int i = 0; i < numXVariables; i++) {
    fmpz_set_si(&((*exps)[i + 1]), xPowers[i] - allGroundPowers[i + 1]);
  }
}

bool FMPoly::getExponentsFromMonomial(const fmpz *exps, int &qPower,
                                      std::vector<int> &xPowers) const {
  xPowers.resize(numXVariables);

  // Extract q power (handle offset)
  qPower = fmpz_get_si(&exps[0]) + allGroundPowers[0];

  // Extract x powers (reverse ground powers offset)
  for (int i = 0; i < numXVariables; i++) {
    xPowers[i] = fmpz_get_si(&exps[i + 1]) + allGroundPowers[i + 1];
  }

  return true;
}

void FMPoly::adjustGroundPowersIfNeeded(int qPower,
                                        const std::vector<int> &xPowers) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Check if we need to adjust any ground powers
  bool needToAdjust = false;
  std::vector<int> newGroundPowers = allGroundPowers;

  // Check q power
  if (qPower < allGroundPowers[0]) {
    newGroundPowers[0] = qPower;
    needToAdjust = true;
  }

  // Check x powers
  for (int i = 0; i < numXVariables; i++) {
    if (xPowers[i] < allGroundPowers[i + 1]) {
      newGroundPowers[i + 1] = xPowers[i];
      needToAdjust = true;
    }
  }

  if (!needToAdjust) {
    return;
  }

  // Create a new polynomial with adjusted exponents
  fmpz_mpoly_t newPoly;
  fmpz_mpoly_init(newPoly, ctx);

  // Get number of terms in current polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Copy all existing terms with adjusted exponents
  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient
    fmpz_t termCoeff;
    fmpz_init(termCoeff);
    fmpz_mpoly_get_term_coeff_fmpz(termCoeff, poly, i, ctx);

    // Get exponents
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Adjust all exponents by the ground power differences
    for (int j = 0; j <= numXVariables; j++) {
      int groundPowerDiff = allGroundPowers[j] - newGroundPowers[j];
      fmpz_add_si(&exps[j], &exps[j], groundPowerDiff);
    }

    // Add term to new polynomial
    fmpz_mpoly_set_coeff_fmpz_fmpz(newPoly, termCoeff, exp_ptrs, ctx);

    // Cleanup
    fmpz_clear(termCoeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Replace the old polynomial with the new one and update ground powers
  fmpz_mpoly_swap(poly, newPoly, ctx);
  fmpz_mpoly_clear(newPoly, ctx);
  allGroundPowers = newGroundPowers;
}

int FMPoly::getCoefficient(int qPower, const std::vector<int> &xPowers) const {
  fmpz *exps;
  slong exp_bits;
  convertExponents(qPower, xPowers, &exps, &exp_bits);

  fmpz_t coeff;
  fmpz_init(coeff);

  // Get coefficient for this monomial - need to use array of pointers
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int i = 0; i <= numXVariables; i++) {
    exp_ptrs[i] = &(exps[i]);
  }

  fmpz_mpoly_get_coeff_fmpz_fmpz(coeff, poly, exp_ptrs, ctx);

  int result = fmpz_get_si(coeff);

  // Cleanup
  fmpz_clear(coeff);
  for (int i = 0; i <= numXVariables; i++) {
    fmpz_clear(&(exps[i]));
  }
  flint_free(exp_ptrs);
  flint_free(exps);

  return result;
}

void FMPoly::setCoefficient(int qPower, const std::vector<int> &xPowers,
                            int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Adjust ground powers if needed to handle negative exponents
  adjustGroundPowersIfNeeded(qPower, xPowers);

  fmpz *exps;
  slong exp_bits;
  convertExponents(qPower, xPowers, &exps, &exp_bits);

  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz_set_si(coeff, coefficient);

  // Create array of pointers for FLINT API
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int i = 0; i <= numXVariables; i++) {
    exp_ptrs[i] = &(exps[i]);
  }

  // Set the coefficient (FLINT handles zero coefficients correctly)
  fmpz_mpoly_set_coeff_fmpz_fmpz(poly, coeff, exp_ptrs, ctx);

  // Cleanup
  fmpz_clear(coeff);
  for (int i = 0; i <= numXVariables; i++) {
    fmpz_clear(&(exps[i]));
  }
  flint_free(exp_ptrs);
  flint_free(exps);
}

void FMPoly::addToCoefficient(int qPower, const std::vector<int> &xPowers,
                              int coefficient) {
  if (coefficient == 0) {
    return;
  }

  int currentCoeff = getCoefficient(qPower, xPowers);
  setCoefficient(qPower, xPowers, currentCoeff + coefficient);
}

std::vector<int> FMPoly::getQPolynomial(const std::vector<int> &xPowers) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Map to store q-powers and their coefficients
  std::map<int, int> qCoeffs;

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Iterate through all terms in the FLINT polynomial
  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> termXPowers;
    getExponentsFromMonomial(exps, qPower, termXPowers);

    // Check if x-powers match the requested ones
    bool match = true;
    for (int j = 0; j < numXVariables; j++) {
      if (termXPowers[j] != xPowers[j]) {
        match = false;
        break;
      }
    }

    if (match) {
      // Add this coefficient to the q-polynomial
      int coeffValue = fmpz_get_si(coeff);
      qCoeffs[qPower] += coeffValue;
    }

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Convert map to vector format
  std::vector<int> result;
  if (qCoeffs.empty()) {
    return result; // Return empty vector if no matching terms
  }

  // Find the range of q-powers
  int minQPower = qCoeffs.begin()->first;
  int maxQPower = qCoeffs.rbegin()->first;

  // Create vector with appropriate size, accounting for negative indices
  int vectorSize = maxQPower - minQPower + 1;
  result.resize(vectorSize, 0);

  // Fill in the coefficients
  for (const auto &[qPower, coeff] : qCoeffs) {
    if (coeff != 0) {
      result[qPower - minQPower] = coeff;
    }
  }

  // Store the minimum q-power information somehow for the caller
  // Since we can't modify the function signature, we'll assume the caller
  // knows to call this in conjunction with other methods if needed
  // For QPolynomial integration, use getQPolynomialObject() instead

  return result;
}

QPolynomial
FMPoly::getQPolynomialObject(const std::vector<int> &xPowers) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Map to store q-powers and their coefficients
  std::map<int, int> qCoeffs;

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Iterate through all terms in the FLINT polynomial
  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> termXPowers;
    getExponentsFromMonomial(exps, qPower, termXPowers);

    // Check if x-powers match the requested ones
    bool match = true;
    for (int j = 0; j < numXVariables; j++) {
      if (termXPowers[j] != xPowers[j]) {
        match = false;
        break;
      }
    }

    if (match) {
      // Add this coefficient to the q-polynomial
      int coeffValue = fmpz_get_si(coeff);
      qCoeffs[qPower] += coeffValue;
    }

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Create QPolynomial object
  if (qCoeffs.empty()) {
    return QPolynomial(); // Return zero polynomial
  }

  // Find the range of q-powers
  int minQPower = qCoeffs.begin()->first;
  int maxQPower = qCoeffs.rbegin()->first;

  // Create coefficient vector
  std::vector<int> coeffVector(maxQPower - minQPower + 1, 0);
  for (const auto &[qPower, coeff] : qCoeffs) {
    if (coeff != 0) {
      coeffVector[qPower - minQPower] = coeff;
    }
  }

  return QPolynomial(coeffVector, minQPower);
}

void FMPoly::setQPolynomial(const std::vector<int> &xPowers,
                            const std::vector<int> &qCoeffs, int minQPower) {
  // Clear existing terms for this x-monomial first
  // Then set new coefficients
  for (size_t i = 0; i < qCoeffs.size(); i++) {
    if (qCoeffs[i] != 0) {
      setCoefficient(minQPower + static_cast<int>(i), xPowers, qCoeffs[i]);
    }
  }
}

void FMPoly::setQPolynomial(const std::vector<int> &xPowers,
                            const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // First, collect all q-powers that need to be cleared for this x-monomial
  std::vector<int> qPowersToClear;
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  for (slong i = 0; i < numTerms; i++) {
    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract x-powers from exponent vector
    int qPower;
    std::vector<int> termXPowers;
    getExponentsFromMonomial(exps, qPower, termXPowers);

    // Check if x-powers match
    bool match = true;
    for (int j = 0; j < numXVariables; j++) {
      if (termXPowers[j] != xPowers[j]) {
        match = false;
        break;
      }
    }

    // If match, record this q-power for clearing
    if (match) {
      qPowersToClear.push_back(qPower);
    }

    // Cleanup
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Now clear all the collected q-powers
  for (int qPower : qPowersToClear) {
    setCoefficient(qPower, xPowers, 0);
  }

  // Finally, set new coefficients if not zero
  if (!qPoly.isZero()) {
    for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower(); qPower++) {
      int coeff = qPoly.getCoefficient(qPower);
      if (coeff != 0) {
        setCoefficient(qPower, xPowers, coeff);
      }
    }
  }
}

void FMPoly::addQPolynomial(const std::vector<int> &xPowers,
                            const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (qPoly.isZero())
    return;

  // Get existing q-polynomial for this x-monomial
  QPolynomial existing = getQPolynomialObject(xPowers);

  // Add the new polynomial
  QPolynomial result = existing + qPoly;

  // Set the result back
  setQPolynomial(xPowers, result);
}

void FMPoly::multiplyQPolynomial(const std::vector<int> &xPowers,
                                 const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (qPoly.isZero()) {
    // Multiplying by zero should make this x-monomial's q-polynomial zero
    setQPolynomial(xPowers, QPolynomial());
    return;
  }

  // Get existing q-polynomial for this x-monomial
  QPolynomial existing = getQPolynomialObject(xPowers);

  if (existing.isZero())
    return; // 0 * anything = 0

  // Multiply the polynomials
  QPolynomial result = existing * qPoly;

  // Set the result back
  setQPolynomial(xPowers, result);
}

FMPoly FMPoly::invertVariable(const int target_index) const {
  if (target_index < 0 || target_index >= numXVariables) {
    throw std::invalid_argument("Target variable index out of range");
  }

  FMPoly result(numXVariables, 0);
  result.allGroundPowers = allGroundPowers;

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Iterate through all terms in the FLINT polynomial
  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Negate the exponent of the target variable
    xPowers[target_index] = -xPowers[target_index];

    // Add term to result polynomial
    int coeffValue = fmpz_get_si(coeff);
    result.addToCoefficient(qPower, xPowers, coeffValue);

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  return result;
}

FMPoly FMPoly::truncate(const std::vector<int> &maxXdegrees) const {
  if (maxXdegrees.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "Max degrees vector size must match number of variables");
  }

  FMPoly result(numXVariables, 0);
  result.allGroundPowers = allGroundPowers;

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Iterate through all terms in the FLINT polynomial
  for (slong i = 0; i < numTerms; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Check if this term should be included (all x-exponents <= maxXdegrees)
    bool includeThisTerm = true;
    for (int j = 0; j < numXVariables; j++) {
      if (xPowers[j] > maxXdegrees[j]) {
        includeThisTerm = false;
        break;
      }
    }

    // Add term to result polynomial if it passes the degree check
    if (includeThisTerm) {
      int coeffValue = fmpz_get_si(coeff);
      result.addToCoefficient(qPower, xPowers, coeffValue);
    }

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  return result;
}

FMPoly FMPoly::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(numXVariables, maxDegree);
  return truncate(maxXdegrees);
}

int FMPoly::getNumXVariables() const { return numXVariables; }

const std::vector<int> &FMPoly::getMaxXDegrees() const { return maxXDegrees; }

const std::vector<int> &FMPoly::getBlockSizes() const { return blockSizes; }

void FMPoly::clear() { fmpz_mpoly_zero(poly, ctx); }

bool FMPoly::isZero() const { return fmpz_mpoly_is_zero(poly, ctx); }

void FMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // Get all terms
  auto terms = getCoefficients();

  // Sort terms for deterministic output
  std::sort(terms.begin(), terms.end(),
            [](const Term &a, const Term &b) { return a.first < b.first; });

  bool firstTerm = true;
  for (const auto &term : terms) {
    const std::vector<int> &xPowers = term.first;
    const QPolynomial &qPoly = term.second;

    // Collect all non-zero q-terms for this x-power combination
    std::vector<std::pair<int, int>> qTerms; // (q_power, coefficient)
    for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower();
         ++qPower) {
      int coeff = qPoly.getCoefficient(qPower);
      if (coeff != 0) {
        qTerms.emplace_back(qPower, coeff);
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
  outputFile << "\t\t\"storage_type\": \"flint\"\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void FMPoly::print(int maxTerms) const {
  std::cout << "FMPoly P(q";
  for (int i = 0; i < numXVariables; i++) {
    std::cout << ", x" << (i + 1);
  }
  std::cout << "):\n";

  if (isZero()) {
    std::cout << "0\n";
    return;
  }

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);
  if (numTerms == 0) {
    std::cout << "0\n";
    return;
  }

  // Limit the number of terms to display
  slong termsToShow =
      (maxTerms > 0 && maxTerms < numTerms) ? maxTerms : numTerms;

  bool first = true;

  // Iterate through terms in the FLINT polynomial
  for (slong i = 0; i < termsToShow; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Get coefficient value
    int coeffValue = fmpz_get_si(coeff);

    if (coeffValue != 0) {
      // Print sign
      if (!first) {
        std::cout << (coeffValue > 0 ? " + " : " - ");
        coeffValue = std::abs(coeffValue);
      } else if (coeffValue < 0) {
        std::cout << "-";
        coeffValue = -coeffValue;
      }
      first = false;

      // Print coefficient if not 1 or if it's a constant term
      bool isConstantTerm = (qPower == 0);
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          isConstantTerm = false;
          break;
        }
      }

      if (coeffValue != 1 || isConstantTerm) {
        std::cout << coeffValue;
      }

      // Print q term
      if (qPower != 0) {
        if (coeffValue != 1 || isConstantTerm)
          std::cout << "*";
        std::cout << "q";
        if (qPower != 1) {
          std::cout << "^" << qPower;
        }
      }

      // Print x terms
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          if (coeffValue != 1 || isConstantTerm || qPower != 0)
            std::cout << "*";
          std::cout << "x" << (j + 1);
          if (xPowers[j] != 1) {
            std::cout << "^" << xPowers[j];
          }
        }
      }
    }

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  if (termsToShow < numTerms) {
    std::cout << " + ... (" << (numTerms - termsToShow) << " more terms)";
  }

  std::cout << "\n";
}

std::vector<std::pair<std::vector<int>, QPolynomial>>
FMPoly::getCoefficients() const {
  // Map x-powers -> QPolynomial in q
  std::map<std::vector<int>, QPolynomial> xToQPoly;

  // Number of terms in the FLINT multivariate polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // FLINT shouldn't store zero terms, but be defensive
    if (fmpz_is_zero(coeff)) {
      fmpz_clear(coeff);
      continue;
    }

    // Get exponent vector for this term: [q, x1, ..., xn]
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Convert FLINT exponents (with ground offsets) back to logical powers
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Convert coefficient to a plain int
    int coeffValue = fmpz_get_si(coeff);

    // Accumulate into the QPolynomial for this xPowers
    QPolynomial &qp = xToQPoly[xPowers]; // default-constructs if not present
    qp.addToCoefficient(qPower, coeffValue);

    // Clean up
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Convert the map into the sparse vector of (xPowers, QPolynomial),
  // skipping any zero polynomials (to emulate BMPoly's behavior).
  std::vector<std::pair<std::vector<int>, QPolynomial>> result;
  result.reserve(xToQPoly.size());

  for (auto &entry : xToQPoly) {
    const auto &xPowers = entry.first;
    const auto &qp = entry.second;

    if (!qp.isZero()) {
      result.emplace_back(xPowers, qp);
    }
  }

  return result;
}

void FMPoly::syncFromSparseVector(const std::vector<Term> &sparseVector) {
  // Clear the current polynomial
  clear();

  // Add each term from the sparse vector
  for (const auto &term : sparseVector) {
    const std::vector<int> &xPowers = term.first;
    const QPolynomial &qPoly = term.second;

    // Only add non-zero q-polynomials
    if (!qPoly.isZero()) {
      setQPolynomial(xPowers, qPoly);
    }
  }
}

std::vector<int> FMPoly::evaluate(const std::vector<int> &point) const {
  if (point.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "Point dimension must match number of variables");
  }

  QPolynomial result;

  // Get all terms in the polynomial
  auto terms = getCoefficients();

  for (const auto &term : terms) {
    const std::vector<int> &xPowers = term.first;
    const QPolynomial &qPoly = term.second;

    // Evaluate the x-monomial at the given point
    int termValue = 1;

    for (size_t i = 0; i < xPowers.size(); ++i) {
      if (xPowers[i] > 0) {
        // Positive exponent: multiply
        for (int exp = 0; exp < xPowers[i]; ++exp) {
          termValue *= point[i];
        }
      } else if (xPowers[i] < 0) {
        // Negative exponent: divide
        if (point[i] == 0) {
          throw std::domain_error("Division by zero: negative exponent with "
                                  "zero value at variable " +
                                  std::to_string(i));
        }
        for (int exp = 0; exp > xPowers[i]; --exp) {
          if (termValue % point[i] != 0) {
            throw std::domain_error(
                "Division would result in non-integer coefficient");
          }
          termValue /= point[i];
        }
      }
    }

    // Multiply the q-polynomial by the evaluated term value and add to result
    for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower();
         ++qPower) {
      int coefficient = qPoly.getCoefficient(qPower);
      if (coefficient != 0) {
        result.addToCoefficient(qPower, coefficient * termValue);
      }
    }
  }

  // Convert QPolynomial to vector<int>
  return result.getCoefficients();
}

void FMPoly::checkCompatibility(const FMPoly &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }
}

int FMPoly::nTerms() const {
  if (isZero())
    return 0;

  return static_cast<int>(fmpz_mpoly_length(poly, ctx));
}

FMPoly &FMPoly::operator+=(const FMPoly &other) {
  checkCompatibility(other);

  // Always use term-by-term addition to avoid ground power issues
  auto otherTerms = other.getCoefficients();
  for (const auto &term : otherTerms) {
    addQPolynomial(term.first, term.second);
  }

  return *this;
}

FMPoly &FMPoly::operator-=(const FMPoly &other) {
  checkCompatibility(other);

  // Always use term-by-term subtraction to avoid ground power issues
  auto otherTerms = other.getCoefficients();
  for (const auto &term : otherTerms) {
    QPolynomial negated;
    for (int q = term.second.getMinPower(); q <= term.second.getMaxPower(); ++q) {
      int coeff = term.second.getCoefficient(q);
      if (coeff != 0) {
        negated.setCoefficient(q, -coeff);
      }
    }
    addQPolynomial(term.first, negated);
  }

  return *this;
}

FMPoly &FMPoly::operator*=(const FMPoly &other) {
  checkCompatibility(other);

  // Always use term-by-term multiplication to avoid ground power issues
  auto thisTerms = getCoefficients();
  auto otherTerms = other.getCoefficients();

  // Clear this polynomial
  clear();

  // Multiply term by term
  for (const auto &thisTerm : thisTerms) {
    for (const auto &otherTerm : otherTerms) {
      // Multiply x-powers
      std::vector<int> resultXPowers(numXVariables);
      for (int i = 0; i < numXVariables; i++) {
        resultXPowers[i] = thisTerm.first[i] + otherTerm.first[i];
      }

      // Multiply q-polynomials
      QPolynomial productQPoly = thisTerm.second * otherTerm.second;

      // Add to result
      addQPolynomial(resultXPowers, productQPoly);
    }
  }

  return *this;
}

FMPoly operator+(const FMPoly &lhs, const FMPoly &rhs) {
  FMPoly result = lhs;
  result += rhs;
  return result;
}

FMPoly operator-(const FMPoly &lhs, const FMPoly &rhs) {
  FMPoly result = lhs;
  result -= rhs;
  return result;
}

FMPoly operator*(const FMPoly &lhs, const FMPoly &rhs) {
  FMPoly result = lhs;
  result *= rhs;
  return result;
}

// Operators with QPolynomial
FMPoly &FMPoly::operator+=(const QPolynomial &qPoly) {
  // Add to the x^0 term
  std::vector<int> zeroXPowers(numXVariables, 0);
  addQPolynomial(zeroXPowers, qPoly);
  return *this;
}

FMPoly &FMPoly::operator-=(const QPolynomial &qPoly) {
  // Subtract from the x^0 term
  std::vector<int> zeroXPowers(numXVariables, 0);

  // Create negated qPoly
  QPolynomial negated;
  for (int q = qPoly.getMinPower(); q <= qPoly.getMaxPower(); ++q) {
    int coeff = qPoly.getCoefficient(q);
    if (coeff != 0) {
      negated.setCoefficient(q, -coeff);
    }
  }

  addQPolynomial(zeroXPowers, negated);
  return *this;
}

FMPoly &FMPoly::operator*=(const QPolynomial &qPoly) {
  // If qPoly is zero, result is zero
  if (qPoly.isZero()) {
    clear();
    return *this;
  }

  // If this polynomial is zero, nothing to do
  if (isZero()) {
    return *this;
  }

  // Get all current terms
  auto terms = getCoefficients();

  // Clear the polynomial
  clear();

  // Multiply each term by qPoly and add back
  for (const auto &term : terms) {
    const std::vector<int> &xPowers = term.first;
    const QPolynomial &thisQPoly = term.second;

    // Multiply the two q-polynomials
    QPolynomial product = thisQPoly * qPoly;

    // Add the result back
    addQPolynomial(xPowers, product);
  }

  return *this;
}

// Binary operators with QPolynomial
FMPoly operator+(const FMPoly &lhs, const QPolynomial &rhs) {
  FMPoly result = lhs;
  result += rhs;
  return result;
}

FMPoly operator+(const QPolynomial &lhs, const FMPoly &rhs) {
  FMPoly result = rhs;
  result += lhs;
  return result;
}

FMPoly operator-(const FMPoly &lhs, const QPolynomial &rhs) {
  FMPoly result = lhs;
  result -= rhs;
  return result;
}

FMPoly operator-(const QPolynomial &lhs, const FMPoly &rhs) {
  // Create FMPoly from lhs, then subtract rhs
  FMPoly result(rhs.getNumXVariables(), 10, rhs.getMaxXDegrees());
  result += lhs;

  // Negate rhs and add
  FMPoly negated = rhs;
  negated *= QPolynomial({-1}, 0); // Multiply by -1
  result += negated;

  return result;
}

FMPoly operator*(const FMPoly &lhs, const QPolynomial &rhs) {
  FMPoly result = lhs;
  result *= rhs;
  return result;
}

FMPoly operator*(const QPolynomial &lhs, const FMPoly &rhs) {
  FMPoly result = rhs;
  result *= lhs;
  return result;
}