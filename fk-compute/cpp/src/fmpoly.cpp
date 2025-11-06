#include "fk/fmpoly.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

FMPoly::FMPoly(int numVariables, int degree, const std::vector<int> &maxDegrees)
    : numXVariables(numVariables), qOffset(0) {

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
    : numXVariables(newNumVariables), qOffset(source.qOffset) {

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
    : numXVariables(other.numXVariables), qOffset(other.qOffset),
      maxXDegrees(other.maxXDegrees), blockSizes(other.blockSizes) {

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
    qOffset = other.qOffset;
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
  fmpz_set_si(&((*exps)[0]), qPower + qOffset);

  // Set x variable exponents
  for (int i = 0; i < numXVariables; i++) {
    if (xPowers[i] < 0) {
      throw std::invalid_argument("Negative x exponents not yet supported");
    }
    fmpz_set_si(&((*exps)[i + 1]), xPowers[i]);
  }
}

bool FMPoly::getExponentsFromMonomial(const fmpz *exps, int &qPower,
                                      std::vector<int> &xPowers) const {
  xPowers.resize(numXVariables);

  // Extract q power (handle offset)
  qPower = fmpz_get_si(&exps[0]) - qOffset;

  // Extract x powers
  for (int i = 0; i < numXVariables; i++) {
    xPowers[i] = fmpz_get_si(&exps[i + 1]);
  }

  return true;
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
  // Check if we need to adjust qOffset for negative q powers
  if (qPower + qOffset < 0) {
    // We need to increase the offset to make all exponents non-negative
    int newOffset = -qPower;
    int offsetDiff = newOffset - qOffset;

    if (offsetDiff > 0) {
      // Create a new polynomial with adjusted exponents
      fmpz_mpoly_t newPoly;
      fmpz_mpoly_init(newPoly, ctx);

      // Get number of terms in current polynomial
      slong numTerms = fmpz_mpoly_length(poly, ctx);

      // Copy all existing terms with adjusted q exponents
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

        // Adjust q exponent by the offset difference
        fmpz_add_si(&exps[0], &exps[0], offsetDiff);

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

      // Replace the old polynomial with the new one
      fmpz_mpoly_swap(poly, newPoly, ctx);
      fmpz_mpoly_clear(newPoly, ctx);

      // Update the offset
      qOffset = newOffset;
    }
  }

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

  return result;
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

FMPoly FMPoly::invertVariable(const int target_index) const {
  if (target_index < 0 || target_index >= numXVariables) {
    throw std::invalid_argument("Target variable index out of range");
  }

  FMPoly result(numXVariables, 0);
  result.qOffset = qOffset;

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
  result.qOffset = qOffset;

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

int FMPoly::getNumXVariables() const { return numXVariables; }

const std::vector<int> &FMPoly::getMaxXDegrees() const { return maxXDegrees; }

const std::vector<int> &FMPoly::getBlockSizes() const { return blockSizes; }

void FMPoly::clear() { fmpz_mpoly_zero(poly, ctx); }

bool FMPoly::isZero() const { return fmpz_mpoly_is_zero(poly, ctx); }

void FMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // This would require iterating through all terms in the FLINT polynomial
  // For now, just write the metadata

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

  // For now, just print that it's a FLINT polynomial
  std::cout << "[FLINT polynomial with " << fmpz_mpoly_length(poly, ctx)
            << " terms]\n";
}

std::vector<int> FMPoly::evaluate(const std::vector<int> &point) const {
  if (point.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "Point dimension must match number of variables");
  }

  // This would involve evaluating the FLINT polynomial at the given x-values
  // and returning the resulting univariate polynomial in q
  std::vector<int> result;
  return result;
}

void FMPoly::checkCompatibility(const FMPoly &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }
}

FMPoly &FMPoly::operator+=(const FMPoly &other) {
  checkCompatibility(other);
  fmpz_mpoly_add(poly, poly, other.poly, ctx);
  return *this;
}

FMPoly &FMPoly::operator-=(const FMPoly &other) {
  checkCompatibility(other);
  fmpz_mpoly_sub(poly, poly, other.poly, ctx);
  return *this;
}

FMPoly &FMPoly::operator*=(const FMPoly &other) {
  checkCompatibility(other);
  fmpz_mpoly_mul(poly, poly, other.poly, ctx);
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
