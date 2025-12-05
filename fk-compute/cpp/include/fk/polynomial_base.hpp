#pragma once

#include "fk/bilvector.hpp"
#include <memory>
#include <string>
#include <vector>

/**
 * Abstract base class for multivariate polynomial implementations
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * This base class defines the common interface that all polynomial implementations
 * (MultivariablePolynomial, FMPoly, BMPoly) must provide, enabling them to be
 * used interchangeably via polymorphism.
 *
 * Design principles:
 * - Pure virtual interface for essential operations
 * - Support for negative exponents in both q and x variables
 * - Sparse or dense storage left to concrete implementations
 * - Template parameters allow different coefficient types and q-polynomial types
 *
 * @tparam CoeffType Type for individual coefficients (default: int)
 * @tparam QPolyType Type for q-polynomials (default: bilvector<CoeffType>)
 */
template <typename CoeffType = int, typename QPolyType = bilvector<CoeffType>>
class PolynomialBase {
public:
  virtual ~PolynomialBase() = default;

  // Type alias for terms (x-powers, q-polynomial pairs)
  using Term = std::pair<std::vector<int>, QPolyType>;

  // ===== Coefficient Access =====

  /**
   * Set coefficient for a specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient New coefficient value
   */
  virtual void setCoefficient(int qPower, const std::vector<int> &xPowers,
                              CoeffType coefficient) = 0;

  /**
   * Add to coefficient for a specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient Value to add
   */
  virtual void addToCoefficient(int qPower, const std::vector<int> &xPowers,
                                CoeffType coefficient) = 0;

  /**
   * Get all coefficients as a vector of (x-powers, q-polynomial) pairs
   * @return Vector of terms, where each term is a pair of x-power vector and q-polynomial
   */
  virtual std::vector<Term> getCoefficients() const = 0;

  // ===== Utilities =====

  /**
   * Clear all coefficients (make polynomial zero)
   */
  virtual void clear() = 0;

  /**
   * Print polynomial in human-readable format
   * @param maxTerms Maximum number of terms to print
   */
  virtual void print(int maxTerms = 10) const = 0;

  /**
   * Export polynomial to JSON format
   * @param fileName Output file name (without .json extension)
   */
  virtual void exportToJson(const std::string &fileName) const = 0;
};
