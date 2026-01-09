#include "fk/fmpoly_class.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace {
fmpz_mpoly_ctx_t g_ctx;
bool g_ctx_initialized = false;
} // namespace

struct FMPolyCtxGuard {
  ~FMPolyCtxGuard() {
    if (g_ctx_initialized) {
      fmpz_mpoly_ctx_clear(g_ctx);
    }
  }
};

static FMPolyCtxGuard g_ctx_guard;

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
  // Initialize FLINT context with numXVariables + 1 variables (q, x1, x2, ..., xn)
  // Each FMPoly must have its own context properly initialized for its number of variables
  fmpz_mpoly_ctx_init(ctx, numXVariables + 1, ORD_LEX);
  fmpz_mpoly_init(poly, ctx);
}

/*
void FMPoly::setupContext() {
  // Initialize FLINT context with numXVariables + 1 variables (q, x1, x2, ...,
  // xn)
  fmpz_mpoly_ctx_init(ctx, numXVariables + 1, ORD_LEX);
  fmpz_mpoly_init(poly, ctx);
}
*/

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

  // Check if coefficient fits in int, warn if not
  if (!fmpz_fits_si(coeff)) {
    std::cerr << "WARNING: Coefficient too large to fit in int, value will be truncated!\n";
    std::cerr << "  Actual value: ";
    fmpz_fprint(stderr, coeff);
    std::cerr << "\n  Monomial: q^" << qPower;
    for (size_t i = 0; i < xPowers.size(); i++) {
      std::cerr << " * x" << (i+1) << "^" << xPowers[i];
    }
    std::cerr << "\n";
  }

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

  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Fast path: only check adjustment if exponent would be out of range
  bool needsAdjustmentCheck = (qPower < allGroundPowers[0]);
  if (!needsAdjustmentCheck) {
    for (int i = 0; i < numXVariables; i++) {
      if (xPowers[i] < allGroundPowers[i + 1]) {
        needsAdjustmentCheck = true;
        break;
      }
    }
  }

  // Only call adjustment if needed
  if (needsAdjustmentCheck) {
    adjustGroundPowersIfNeeded(qPower, xPowers);
  }

  // Convert exponents with ground power offset
  fmpz *exps;
  slong exp_bits;
  convertExponents(qPower, xPowers, &exps, &exp_bits);

  // Create coefficient to add
  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz_set_si(coeff, coefficient);

  // Create array of pointers for FLINT API
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int i = 0; i <= numXVariables; i++) {
    exp_ptrs[i] = &(exps[i]);
  }

  // Get current coefficient and add to it
  fmpz_t current;
  fmpz_init(current);
  fmpz_mpoly_get_coeff_fmpz_fmpz(current, poly, exp_ptrs, ctx);
  fmpz_add(current, current, coeff);

  // Set the new coefficient
  fmpz_mpoly_set_coeff_fmpz_fmpz(poly, current, exp_ptrs, ctx);

  // Cleanup
  fmpz_clear(coeff);
  fmpz_clear(current);
  for (int i = 0; i <= numXVariables; i++) {
    fmpz_clear(&(exps[i]));
  }
  flint_free(exp_ptrs);
  flint_free(exps);
}

QPolynomial
FMPoly::getQPolynomialObject(const std::vector<int> &xPowers) const {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Create result QPolynomial
  QPolynomial result;

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
      // Add this coefficient to the q-polynomial using fmpz directly
      result.addToCoefficientFmpz(qPower, coeff);
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
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
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
    // Use fmpz directly to avoid integer overflow
    fmpz_t coeff;
    fmpz_init(coeff);

    for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower();
         qPower++) {
      qPoly.getCoefficientFmpz(coeff, qPower);
      if (!fmpz_is_zero(coeff)) {
        // Adjust ground powers if needed to handle negative exponents
        adjustGroundPowersIfNeeded(qPower, xPowers);

        fmpz *exps;
        slong exp_bits;
        convertExponents(qPower, xPowers, &exps, &exp_bits);

        // Create array of pointers for FLINT API
        fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
        for (int i = 0; i <= numXVariables; i++) {
          exp_ptrs[i] = &(exps[i]);
        }

        // Set the coefficient (FLINT handles zero coefficients correctly)
        fmpz_mpoly_set_coeff_fmpz_fmpz(poly, coeff, exp_ptrs, ctx);

        // Cleanup
        for (int i = 0; i <= numXVariables; i++) {
          fmpz_clear(&(exps[i]));
        }
        flint_free(exp_ptrs);
        flint_free(exps);
      }
    }

    fmpz_clear(coeff);
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

  // Use fmpz directly to avoid integer overflow
  fmpz_t coeff;
  fmpz_init(coeff);

  for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower();
       qPower++) {
    qPoly.getCoefficientFmpz(coeff, qPower);
    if (!fmpz_is_zero(coeff)) {
      // Fast path: only check adjustment if exponent would be out of range
      bool needsAdjustmentCheck = (qPower < allGroundPowers[0]);
      if (!needsAdjustmentCheck) {
        for (int i = 0; i < numXVariables; i++) {
          if (xPowers[i] < allGroundPowers[i + 1]) {
            needsAdjustmentCheck = true;
            break;
          }
        }
      }

      // Only call adjustment if needed
      if (needsAdjustmentCheck) {
        adjustGroundPowersIfNeeded(qPower, xPowers);
      }

      // Convert exponents with ground power offset
      fmpz *exps;
      slong exp_bits;
      convertExponents(qPower, xPowers, &exps, &exp_bits);

      // Create array of pointers for FLINT API
      fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
      for (int i = 0; i <= numXVariables; i++) {
        exp_ptrs[i] = &(exps[i]);
      }

      // Get current coefficient and add to it
      fmpz_t current;
      fmpz_init(current);
      fmpz_mpoly_get_coeff_fmpz_fmpz(current, poly, exp_ptrs, ctx);
      fmpz_add(current, current, coeff);

      // Set the new coefficient
      fmpz_mpoly_set_coeff_fmpz_fmpz(poly, current, exp_ptrs, ctx);

      // Cleanup
      fmpz_clear(current);
      for (int i = 0; i <= numXVariables; i++) {
        fmpz_clear(&(exps[i]));
      }
      flint_free(exp_ptrs);
      flint_free(exps);
    }
  }

  fmpz_clear(coeff);
}

FMPoly FMPoly::truncate(const std::vector<int> &maxXdegrees) const {
  if (maxXdegrees.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "Max degrees vector size must match number of variables");
  }

  // Result: same x-structure & ground powers, but coefficients truncated
  FMPoly result(numXVariables, 0, maxXdegrees);
  result.allGroundPowers = allGroundPowers;

  slong numTerms = fmpz_mpoly_length(poly, ctx);
  if (numTerms == 0) {
    return result;
  }

  // Reuse a single coeff and exponent buffer for all terms
  fmpz_t coeff;
  fmpz_init(coeff);

  // exps holds (q, x1, ..., xn) in FLINT's stored form (nonnegative)
  fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_init(&exps[j]);
    exp_ptrs[j] = &exps[j];
  }

  // Pre-compute max stored values to avoid repeated additions in the loop
  // stored_exp[j+1] = x_j - allGroundPowers[j+1], so max_stored = maxXdegrees[j] - allGroundPowers[j+1]
  std::vector<long> max_stored_x(numXVariables);
  for (int j = 0; j < numXVariables; ++j) {
    max_stored_x[j] = maxXdegrees[j] - allGroundPowers[j + 1];
  }

  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of term i
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);
    if (fmpz_is_zero(coeff)) {
      continue;
    }

    // Get stored exponents of term i
    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Fast bounds check with pre-computed limits
    bool include = true;
    for (int j = 0; j < numXVariables; ++j) {
      // Compare stored exponent directly against max_stored_x
      if (fmpz_cmp_si(&exps[j + 1], max_stored_x[j]) > 0) {
        include = false;
        break;
      }
    }

    if (!include) {
      continue;
    }

    // Keep this term: set coefficient in result with same stored exponents
    // (same ctx, same allGroundPowers, so stored exps are valid as-is)
    fmpz_mpoly_set_coeff_fmpz_fmpz(result.poly, coeff, exp_ptrs, ctx);
  }

  // Cleanup
  fmpz_clear(coeff);
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_clear(&exps[j]);
  }
  flint_free(exp_ptrs);
  flint_free(exps);

  return result;
}

FMPoly FMPoly::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(numXVariables, maxDegree);
  return truncate(maxXdegrees);
}

int FMPoly::getNumXVariables() const { return numXVariables; }

void FMPoly::clear() { fmpz_mpoly_zero(poly, ctx); }

bool FMPoly::isZero() const { return fmpz_mpoly_is_zero(poly, ctx); }

void FMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Collect all terms grouped by x-powers
  std::map<std::vector<int>, std::vector<std::pair<int, std::string>>> xToQTerms;

  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Skip zero coefficients (shouldn't happen in FLINT, but be defensive)
    if (fmpz_is_zero(coeff)) {
      fmpz_clear(coeff);
      continue;
    }

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Convert coefficient to string to preserve full precision
    char *coeffStr = fmpz_get_str(NULL, 10, coeff);
    xToQTerms[xPowers].emplace_back(qPower, std::string(coeffStr));
    flint_free(coeffStr);

    // Clean up
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Output terms in sorted order
  bool firstTerm = true;
  for (const auto &entry : xToQTerms) {
    const std::vector<int> &xPowers = entry.first;
    const auto &qTerms = entry.second;

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
        // Output coefficient as string to preserve precision
        outputFile << "{\"q\": " << qTerms[i].first
                   << ", \"c\": \"" << qTerms[i].second << "\"}";
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

    if (!fmpz_is_zero(coeff)) {
      // Determine sign and get absolute value
      int sgn = fmpz_sgn(coeff);
      fmpz_t absCoeff;
      fmpz_init(absCoeff);
      fmpz_abs(absCoeff, coeff);

      // Print sign
      if (!first) {
        std::cout << (sgn > 0 ? " + " : " - ");
      } else if (sgn < 0) {
        std::cout << "-";
      }
      first = false;

      // Check if this is a constant term
      bool isConstantTerm = (qPower == 0);
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          isConstantTerm = false;
          break;
        }
      }

      // Print coefficient if not 1 or if it's a constant term
      bool coeffIsOne = fmpz_is_one(absCoeff);
      if (!coeffIsOne || isConstantTerm) {
        // Convert to string to handle big coefficients
        char *coeffStr = fmpz_get_str(NULL, 10, absCoeff);
        std::cout << coeffStr;
        flint_free(coeffStr);
      }

      // Print q term
      if (qPower != 0) {
        if (!coeffIsOne || isConstantTerm)
          std::cout << "*";
        std::cout << "q";
        if (qPower != 1) {
          std::cout << "^" << qPower;
        }
      }

      // Print x terms
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          if (!coeffIsOne || isConstantTerm || qPower != 0)
            std::cout << "*";
          std::cout << "x" << (j + 1);
          if (xPowers[j] != 1) {
            std::cout << "^" << xPowers[j];
          }
        }
      }

      fmpz_clear(absCoeff);
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

  // Allocate buffers ONCE outside the loop
  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_init(&exps[j]);
    exp_ptrs[j] = &exps[j];
  }

  // Loop through terms, reusing buffers
  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of this term
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // FLINT shouldn't store zero terms, but be defensive
    if (fmpz_is_zero(coeff)) {
      continue;
    }

    // Get exponent vector for this term: [q, x1, ..., xn]
    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Convert FLINT exponents (with ground offsets) back to logical powers
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Accumulate into the QPolynomial for this xPowers using fmpz directly
    QPolynomial &qp = xToQPoly[xPowers]; // default-constructs if not present
    qp.addToCoefficientFmpz(qPower, coeff);
  }

  // Clean up buffers ONCE after the loop
  fmpz_clear(coeff);
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_clear(&exps[j]);
  }
  flint_free(exp_ptrs);
  flint_free(exps);

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

void FMPoly::checkCompatibility(const FMPoly &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }
}

FMPoly &FMPoly::operator+=(const FMPoly &other) {
  checkCompatibility(other);

  // Pre-adjust ground powers once based on other's ground powers
  // This avoids multiple adjustments during term addition
  std::vector<int> minPowers = allGroundPowers;
  for (int i = 0; i <= numXVariables; i++) {
    if (other.allGroundPowers[i] < minPowers[i]) {
      minPowers[i] = other.allGroundPowers[i];
    }
  }

  // Extract coefficients ONCE - this is expensive so we reuse it
  auto otherTerms = other.getCoefficients();

  // Do a single adjustment if needed
  if (minPowers != allGroundPowers) {
    if (!otherTerms.empty()) {
      const auto &firstTerm = otherTerms.front();
      int minQPower = firstTerm.second.getMinPower();

      // Find minimum q power across all terms
      for (const auto &term : otherTerms) {
        if (term.second.getMinPower() < minQPower) {
          minQPower = term.second.getMinPower();
        }
      }

      // Find minimum x powers
      std::vector<int> minXPowers(numXVariables,
                                  std::numeric_limits<int>::max());
      for (const auto &term : otherTerms) {
        for (int i = 0; i < numXVariables; i++) {
          if (term.first[i] < minXPowers[i]) {
            minXPowers[i] = term.first[i];
          }
        }
      }

      // Trigger adjustment with minimum powers
      adjustGroundPowersIfNeeded(minQPower, minXPowers);
    }
  }

  // Now add all terms - reusing otherTerms from above
  for (const auto &term : otherTerms) {
    addQPolynomial(term.first, term.second);
  }

  return *this;
}

FMPoly &FMPoly::operator*=(const FMPoly &other) {
  checkCompatibility(other);

  // Quick zero checks: 0 * anything = 0
  if (isZero() || other.isZero()) {
    clear();
    // For a zero polynomial the ground powers don't really matter;
    // leaving them as-is is fine.
    return *this;
  }

  // Multiply underlying FLINT polynomials using the existing context.
  // We use a temporary to avoid aliasing issues (e.g. self-multiplication).
  fmpz_mpoly_t prod;
  fmpz_mpoly_init(prod, ctx);
  fmpz_mpoly_mul(prod, poly, other.poly, ctx);

  // Swap result into this->poly and clean up
  fmpz_mpoly_swap(poly, prod, ctx);
  fmpz_mpoly_clear(prod, ctx);

  // Update ground powers: product exponents are sums of absolute exponents,
  // and stored exponents are sums of stored exponents, so the new ground
  // vector is just the componentwise sum.
  for (int i = 0; i <= numXVariables; ++i) {
    allGroundPowers[i] += other.allGroundPowers[i];
  }

  // maxXDegrees / blockSizes are just "metadata" compatibility; we can
  // leave them unchanged, since FLINT doesn't use them.
  return *this;
}

FMPoly FMPoly::operator*(const QPolynomial &qPoly) const {
  // Result has the same x-structure metadata as *this
  FMPoly result(numXVariables, 0, maxXDegrees);

  // Zero short-circuits
  if (isZero() || qPoly.isZero()) {
    return result; // already zero poly with same shape
  }

  // Build an FMPoly that represents "qPoly in the first variable, no
  // x-dependence"
  FMPoly factor(numXVariables, 0, maxXDegrees);
  std::vector<int> zeroXPowers(numXVariables, 0);
  factor.addQPolynomial(zeroXPowers, qPoly);

  // Now use the fast FLINT-based FMPoly multiplication
  result = *this;
  result *= factor; // this is your fmpz_mpoly_mul-based operator*=

  return result;
}
