#include "fk/qpolynomial.hpp"
#include <algorithm>
#include <flint/fmpz.h>
#include <iostream>
#include <sstream>

QPolynomial::QPolynomial() : minPower(0) { fmpz_poly_init(poly); }

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

QPolynomial::~QPolynomial() { fmpz_poly_clear(poly); }

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
    if (fmpz_poly_is_zero(poly))
      return;

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
  if (coeff == 0)
    return;

  int currentCoeff = getCoefficient(power);
  setCoefficient(power, currentCoeff + coeff);
}


void QPolynomial::setFromCoefficients(const std::vector<int> &coeffs,
                                      int minQPower) {
  fmpz_poly_zero(poly);
  minPower = minQPower;

  for (size_t i = 0; i < coeffs.size(); i++) {
    if (coeffs[i] != 0) {
      fmpz_poly_set_coeff_si(poly, i, coeffs[i]);
    }
  }
}

int QPolynomial::getMinPower() const {
  if (isZero())
    return 0;

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
  if (isZero())
    return minPower - 1;

  slong degree = fmpz_poly_degree(poly);
  return minPower + degree;
}

int QPolynomial::getDegree() const {
  if (isZero())
    return -1;
  return fmpz_poly_degree(poly);
}

bool QPolynomial::isZero() const { return fmpz_poly_is_zero(poly); }

void QPolynomial::clear() {
  fmpz_poly_zero(poly);
  minPower = 0;
}

int QPolynomial::evaluate(int q) const {
  if (isZero())
    return 0;

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

      if (!first && c > 0)
        oss << " + ";
      else if (c < 0)
        oss << " - ";

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
  // (call getMinPower() if you need to know which exponent result[0]
  // corresponds to)
  return result;
}

int QPolynomial::nTerms() const {
  if (isZero())
    return 0;

  int count = 0;
  slong length = fmpz_poly_length(poly);

  for (slong i = 0; i < length; ++i) {
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_poly_get_coeff_fmpz(coeff, poly, i);

    if (!fmpz_is_zero(coeff)) {
      ++count;
    }

    fmpz_clear(coeff);
  }

  return count;
}

QPolynomial &QPolynomial::operator+=(const QPolynomial &other) {
  if (other.isZero())
    return *this;
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
  if (other.isZero())
    return *this;

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