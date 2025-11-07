#include "fk/fmpoly.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sstream>

// QPolynomial Implementation

QPolynomial::QPolynomial() : minPower(0) {
    fmpz_poly_init(poly);
}

QPolynomial::QPolynomial(const std::vector<int> &coeffs, int minQPower)
    : minPower(minQPower) {
    fmpz_poly_init(poly);
    setFromCoefficients(coeffs, minQPower);
}

QPolynomial::QPolynomial(const QPolynomial &other) : minPower(other.minPower) {
    fmpz_poly_init(poly);
    fmpz_poly_set(poly, other.poly);
}

QPolynomial &QPolynomial::operator=(const QPolynomial &other) {
    if (this != &other) {
        fmpz_poly_set(poly, other.poly);
        minPower = other.minPower;
    }
    return *this;
}

QPolynomial::~QPolynomial() {
    fmpz_poly_clear(poly);
}

int QPolynomial::getCoefficient(int power) const {
    int index = power - minPower;
    if (index < 0 || index >= fmpz_poly_length(poly)) {
        return 0;
    }

    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_poly_get_coeff_fmpz(coeff, poly, index);
    int result = fmpz_get_si(coeff);
    fmpz_clear(coeff);
    return result;
}

void QPolynomial::setCoefficient(int power, int coeff) {
    if (coeff == 0) {
        // If setting to zero and polynomial is empty, just return
        if (fmpz_poly_is_zero(poly)) return;

        int index = power - minPower;
        if (index >= 0 && index < fmpz_poly_length(poly)) {
            fmpz_poly_set_coeff_si(poly, index, 0);
        }
        return;
    }

    // If polynomial is currently zero, initialize with this power
    if (fmpz_poly_is_zero(poly)) {
        minPower = power;
        fmpz_poly_set_coeff_si(poly, 0, coeff);
        return;
    }

    int index = power - minPower;

    if (index < 0) {
        // Need to shift polynomial to accommodate negative index
        int shift = -index;
        fmpz_poly_t temp;
        fmpz_poly_init(temp);

        // Create polynomial with leading zeros
        fmpz_poly_fit_length(temp, fmpz_poly_length(poly) + shift);

        // Copy existing coefficients shifted right
        for (slong i = 0; i < fmpz_poly_length(poly); i++) {
            fmpz_t coeff_val;
            fmpz_init(coeff_val);
            fmpz_poly_get_coeff_fmpz(coeff_val, poly, i);
            fmpz_poly_set_coeff_fmpz(temp, i + shift, coeff_val);
            fmpz_clear(coeff_val);
        }

        // Set the new coefficient
        fmpz_poly_set_coeff_si(temp, 0, coeff);

        fmpz_poly_swap(poly, temp);
        fmpz_poly_clear(temp);

        minPower = power;
    } else {
        fmpz_poly_set_coeff_si(poly, index, coeff);
    }
}

void QPolynomial::addToCoefficient(int power, int coeff) {
    if (coeff == 0) return;

    int currentCoeff = getCoefficient(power);
    setCoefficient(power, currentCoeff + coeff);
}

std::vector<int> QPolynomial::getCoefficients() const {
    std::vector<int> result;
    slong length = fmpz_poly_length(poly);

    for (slong i = 0; i < length; i++) {
        fmpz_t coeff;
        fmpz_init(coeff);
        fmpz_poly_get_coeff_fmpz(coeff, poly, i);
        result.push_back(fmpz_get_si(coeff));
        fmpz_clear(coeff);
    }

    return result;
}

void QPolynomial::setFromCoefficients(const std::vector<int> &coeffs, int minQPower) {
    fmpz_poly_zero(poly);
    minPower = minQPower;

    for (size_t i = 0; i < coeffs.size(); i++) {
        if (coeffs[i] != 0) {
            fmpz_poly_set_coeff_si(poly, i, coeffs[i]);
        }
    }
}

int QPolynomial::getMinPower() const {
    if (isZero()) return 0;

    slong length = fmpz_poly_length(poly);
    for (slong i = 0; i < length; i++) {
        fmpz_t coeff;
        fmpz_init(coeff);
        fmpz_poly_get_coeff_fmpz(coeff, poly, i);
        if (!fmpz_is_zero(coeff)) {
            fmpz_clear(coeff);
            return minPower + i;
        }
        fmpz_clear(coeff);
    }
    return minPower;
}

int QPolynomial::getMaxPower() const {
    if (isZero()) return minPower - 1;

    slong degree = fmpz_poly_degree(poly);
    return minPower + degree;
}

int QPolynomial::getDegree() const {
    if (isZero()) return -1;
    return fmpz_poly_degree(poly);
}

bool QPolynomial::isZero() const {
    return fmpz_poly_is_zero(poly);
}

void QPolynomial::clear() {
    fmpz_poly_zero(poly);
    minPower = 0;
}

int QPolynomial::evaluate(int q) const {
    if (isZero()) return 0;

    fmpz_t q_fmpz, result;
    fmpz_init(q_fmpz);
    fmpz_init(result);

    fmpz_set_si(q_fmpz, q);
    fmpz_poly_evaluate_fmpz(result, poly, q_fmpz);

    if (minPower != 0) {
        fmpz_t q_power;
        fmpz_init(q_power);
        fmpz_pow_ui(q_power, q_fmpz, abs(minPower));

        if (minPower > 0) {
            fmpz_mul(result, result, q_power);
        } else {
            fmpz_divexact(result, result, q_power);
        }

        fmpz_clear(q_power);
    }

    int ret = fmpz_get_si(result);
    fmpz_clear(q_fmpz);
    fmpz_clear(result);
    return ret;
}

void QPolynomial::print() const {
    if (isZero()) {
        std::cout << "0" << std::endl;
        return;
    }

    std::ostringstream oss;
    bool first = true;
    slong length = fmpz_poly_length(poly);

    for (slong i = length - 1; i >= 0; i--) {
        fmpz_t coeff;
        fmpz_init(coeff);
        fmpz_poly_get_coeff_fmpz(coeff, poly, i);

        if (!fmpz_is_zero(coeff)) {
            int c = fmpz_get_si(coeff);
            int power = minPower + i;

            if (!first && c > 0) oss << " + ";
            else if (c < 0) oss << " - ";

            if (abs(c) != 1 || power == 0) {
                oss << abs(c);
            }

            if (power != 0) {
                oss << "q";
                if (power != 1) {
                    oss << "^" << power;
                }
            }

            first = false;
        }
        fmpz_clear(coeff);
    }

    std::cout << oss.str() << std::endl;
}

std::vector<int> QPolynomial::getCoefficients() const {
    std::vector<int> result;

    // Length of the underlying FLINT polynomial: highest internal index + 1
    slong length = fmpz_poly_length(poly);

    // Reserve to avoid reallocations
    result.reserve(static_cast<std::size_t>(length));

    // For each internal index i, read the coefficient of q^(minPower + i)
    for (slong i = 0; i < length; ++i) {
        fmpz_t coeff;
        fmpz_init(coeff);

        // Get coefficient of q^i in the FLINT polynomial
        fmpz_poly_get_coeff_fmpz(coeff, poly, i);

        // Convert to a plain int and append
        result.push_back(fmpz_get_si(coeff));

        fmpz_clear(coeff);
    }

    // Semantics: result[i] is the coefficient of q^(minPower + i)
    // (call getMinPower() if you need to know which exponent result[0] corresponds to)
    return result;
}

QPolynomial &QPolynomial::operator+=(const QPolynomial &other) {
    if (other.isZero()) return *this;
    if (this->isZero()) {
        *this = other;
        return *this;
    }

    // Find the new minimum power
    int newMinPower = std::min(this->minPower, other.minPower);

    // Calculate how much each polynomial needs to be shifted
    int thisShift = this->minPower - newMinPower;
    int otherShift = other.minPower - newMinPower;

    // Create temporary polynomials with proper alignment
    fmpz_poly_t thisAligned, otherAligned;
    fmpz_poly_init(thisAligned);
    fmpz_poly_init(otherAligned);

    // Align this polynomial
    if (thisShift > 0) {
        fmpz_poly_shift_left(thisAligned, this->poly, thisShift);
    } else {
        fmpz_poly_set(thisAligned, this->poly);
    }

    // Align other polynomial
    if (otherShift > 0) {
        fmpz_poly_shift_left(otherAligned, other.poly, otherShift);
    } else {
        fmpz_poly_set(otherAligned, other.poly);
    }

    // Perform addition
    fmpz_poly_add(this->poly, thisAligned, otherAligned);
    this->minPower = newMinPower;

    fmpz_poly_clear(thisAligned);
    fmpz_poly_clear(otherAligned);

    return *this;
}

QPolynomial &QPolynomial::operator-=(const QPolynomial &other) {
    if (other.isZero()) return *this;

    // Handle case where this is zero
    if (this->isZero()) {
        *this = other;
        // Negate all coefficients
        fmpz_poly_neg(this->poly, this->poly);
        return *this;
    }

    // Find the new minimum power
    int newMinPower = std::min(this->minPower, other.minPower);

    // Calculate how much each polynomial needs to be shifted
    int thisShift = this->minPower - newMinPower;
    int otherShift = other.minPower - newMinPower;

    // Create temporary polynomials with proper alignment
    fmpz_poly_t thisAligned, otherAligned;
    fmpz_poly_init(thisAligned);
    fmpz_poly_init(otherAligned);

    // Align this polynomial
    if (thisShift > 0) {
        fmpz_poly_shift_left(thisAligned, this->poly, thisShift);
    } else {
        fmpz_poly_set(thisAligned, this->poly);
    }

    // Align other polynomial
    if (otherShift > 0) {
        fmpz_poly_shift_left(otherAligned, other.poly, otherShift);
    } else {
        fmpz_poly_set(otherAligned, other.poly);
    }

    // Perform subtraction
    fmpz_poly_sub(this->poly, thisAligned, otherAligned);
    this->minPower = newMinPower;

    fmpz_poly_clear(thisAligned);
    fmpz_poly_clear(otherAligned);

    return *this;
}

QPolynomial &QPolynomial::operator*=(const QPolynomial &other) {
    if (this->isZero() || other.isZero()) {
        this->clear();
        return *this;
    }

    fmpz_poly_mul(this->poly, this->poly, other.poly);
    this->minPower += other.minPower;

    return *this;
}

QPolynomial operator+(const QPolynomial &lhs, const QPolynomial &rhs) {
    QPolynomial result = lhs;
    result += rhs;
    return result;
}

QPolynomial operator-(const QPolynomial &lhs, const QPolynomial &rhs) {
    QPolynomial result = lhs;
    result -= rhs;
    return result;
}

QPolynomial operator*(const QPolynomial &lhs, const QPolynomial &rhs) {
    QPolynomial result = lhs;
    result *= rhs;
    return result;
}



// FMPoly Implementation

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
    : numXVariables(other.numXVariables), allGroundPowers(other.allGroundPowers),
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

void FMPoly::adjustGroundPowersIfNeeded(int qPower, const std::vector<int> &xPowers) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument("X powers vector size must match number of variables");
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
    fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
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

QPolynomial FMPoly::getQPolynomialObject(const std::vector<int> &xPowers) const {
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

void FMPoly::setQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Clear existing coefficients for this x-monomial
  // We need to iterate through existing terms and clear those matching xPowers
  // For simplicity, we'll clear by setting to zero for the range we're about to set

  if (qPoly.isZero()) {
    // If setting to zero polynomial, we should clear existing terms
    // This is complex to implement efficiently, so for now we'll skip
    return;
  }

  // Get the range of powers in the QPolynomial
  int minPower = qPoly.getMinPower();
  int maxPower = qPoly.getMaxPower();

  // Clear existing terms in this range
  for (int qPower = minPower; qPower <= maxPower; qPower++) {
    setCoefficient(qPower, xPowers, 0);
  }

  // Set new coefficients
  for (int qPower = minPower; qPower <= maxPower; qPower++) {
    int coeff = qPoly.getCoefficient(qPower);
    if (coeff != 0) {
      setCoefficient(qPower, xPowers, coeff);
    }
  }
}

void FMPoly::addQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (qPoly.isZero()) return;

  // Get existing q-polynomial for this x-monomial
  QPolynomial existing = getQPolynomialObject(xPowers);

  // Add the new polynomial
  QPolynomial result = existing + qPoly;

  // Set the result back
  setQPolynomial(xPowers, result);
}

void FMPoly::multiplyQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly) {
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

  if (existing.isZero()) return; // 0 * anything = 0

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

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);
  if (numTerms == 0) {
    std::cout << "0\n";
    return;
  }

  // Limit the number of terms to display
  slong termsToShow = (maxTerms > 0 && maxTerms < numTerms) ? maxTerms : numTerms;

  bool first = true;

  // Iterate through terms in the FLINT polynomial
  for (slong i = 0; i < termsToShow; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
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
        if (coeffValue != 1 || isConstantTerm) std::cout << "*";
        std::cout << "q";
        if (qPower != 1) {
          std::cout << "^" << qPower;
        }
      }

      // Print x terms
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          if (coeffValue != 1 || isConstantTerm || qPower != 0) std::cout << "*";
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
    QPolynomial &qp = xToQPoly[xPowers];  // default-constructs if not present
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
