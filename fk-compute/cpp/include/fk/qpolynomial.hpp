#pragma once

#include <flint/fmpz_poly.h>
#include <string>
#include <vector>

/**
 * QPolynomial: A FLINT-based univariate polynomial class for q polynomials
 *
 * Represents polynomials of the form:
 * P(q) = Σᵢ coeffᵢ × q^(i + minPower)
 *
 * Features:
 * - FLINT backend for optimal performance
 * - Support for negative q exponents via offset tracking
 * - Easy-to-use wrapper around fmpz_poly_t
 */
class QPolynomial {
private:
  fmpz_poly_t poly; // FLINT univariate polynomial
  int minPower;     // Minimum power of q (for negative exponents)

public:
  /**
   * Default constructor - creates zero polynomial
   */
  QPolynomial();

  /**
   * Constructor from coefficient vector
   * @param coeffs Vector of coefficients
   * @param minQPower Minimum power of q (default 0)
   */
  QPolynomial(const std::vector<int> &coeffs, int minQPower = 0);

  /**
   * Copy constructor
   */
  QPolynomial(const QPolynomial &other);

  /**
   * Assignment operator
   */
  QPolynomial &operator=(const QPolynomial &other);

  /**
   * Destructor
   */
  ~QPolynomial();

  /**
   * Get coefficient for q^power
   * @param power Power of q
   * @return Coefficient value
   */
  int getCoefficient(int power) const;

  /**
   * Set coefficient for q^power
   * @param power Power of q
   * @param coeff Coefficient value
   */
  void setCoefficient(int power, int coeff);

  /**
   * Add to coefficient for q^power
   * @param power Power of q
   * @param coeff Value to add
   */
  void addToCoefficient(int power, int coeff);

  /**
   * Get coefficients as vector
   * @return Vector where index i represents coefficient of q^(i + minPower)
   */
  std::vector<int> getCoefficients() const;

  /**
   * Set from coefficient vector
   * @param coeffs Vector of coefficients
   * @param minQPower Minimum power of q
   */
  void setFromCoefficients(const std::vector<int> &coeffs, int minQPower = 0);

  /**
   * Get the numer of terms of poly
   */
  int nTerms() const;

  /**
   * Get minimum power of q
   */
  int getMinPower() const;

  /**
   * Get maximum power of q
   */
  int getMaxPower() const;

  /**
   * Get degree of polynomial (-1 for zero polynomial)
   */
  int getDegree() const;

  /**
   * Check if polynomial is zero
   */
  bool isZero() const;

  /**
   * Clear polynomial (set to zero)
   */
  void clear();

  /**
   * Evaluate polynomial at a given value
   * @param q Value to evaluate at
   * @return Result of evaluation
   */
  int evaluate(int q) const;

  /**
   * Print polynomial in human-readable format
   */
  void print() const;

  /**
   * Arithmetic operators
   */
  QPolynomial &operator+=(const QPolynomial &other);
  QPolynomial &operator-=(const QPolynomial &other);
  QPolynomial &operator*=(const QPolynomial &other);

  /**
   * Friend functions for binary operations
   */
  friend QPolynomial operator+(const QPolynomial &lhs, const QPolynomial &rhs);
  friend QPolynomial operator-(const QPolynomial &lhs, const QPolynomial &rhs);
  friend QPolynomial operator*(const QPolynomial &lhs, const QPolynomial &rhs);

  /**
   * Get access to the underlying FLINT polynomial (for advanced operations)
   */
  const fmpz_poly_t &getFlintPoly() const { return poly; }
};