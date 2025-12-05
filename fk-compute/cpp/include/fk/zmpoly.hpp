#pragma once

#include "fk/polynomial_base.hpp"
#include "fk/qpolynomial.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * Custom hash function for vector<int> keys
 */
struct VectorHashZM {
  std::size_t operator()(const std::vector<int> &v) const;
};

/**
 * ZMPoly: A polynomial in variables q, x₁, x₂, ..., xₙ with arbitrary precision integer coefficients
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - Arbitrary precision integer coefficients (using FLINT's fmpz_t via QPolynomial)
 * - Arbitrary positive/negative q powers
 * - Configurable number of x variables
 * - Sparse storage using unordered_map (only non-zero coefficients)
 * - Support for negative x exponents
 * - Efficient indexing and arithmetic operations
 * - Same interface as MultivariablePolynomial for drop-in replacement
 */
class ZMPoly : public PolynomialBase<int, QPolynomial> {
private:
  int numXVariables;            // Number of x variables (components)
  std::vector<int> maxXDegrees; // Max degrees (advisory only, not enforced)
  std::vector<int> blockSizes; // Block sizes (advisory only, for compatibility)

  // Sparse storage: map from x-exponent vector to q-polynomial
  std::unordered_map<std::vector<int>, QPolynomial, VectorHashZM> coeffs_;

  // Prune zero coefficients from the map
  void pruneZeros();

public:
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable
   * (advisory only)
   */
  ZMPoly(int numVariables, int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Constructor to increase the number of variables from another polynomial
   * @param source The source polynomial with fewer variables
   * @param newNumVariables The new number of variables (must be >= source's numVariables)
   * @param targetVariableIndex The index (0-based) where the source's variable should be mapped
   * @param degree Maximum degree for each new x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable (advisory only)
   */
  ZMPoly(const ZMPoly &source,
         int newNumVariables,
         int targetVariableIndex,
         int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Set coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient New coefficient value
   */
  void setCoefficient(int qPower, const std::vector<int> &xPowers,
                      int coefficient) override;

  /**
   * Add to coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient Value to add
   */
  void addToCoefficient(int qPower, const std::vector<int> &xPowers,
                        int coefficient) override;

  /**
   * Truncate multivariable polynomial to given degrees
   */
  ZMPoly truncate(const std::vector<int> & maxXdegrees) const;

  /**
   * Truncate multivariable polynomial to the same degree for all x variables
   * @param maxDegree Maximum degree to keep for all x variables
   */
  ZMPoly truncate(int maxDegree) const;


  /**
   * Get coefficients as a vector of (x-powers, q-polynomial) pairs
   * Compatible with other polynomial classes' interface
   */
  using Term = std::pair<std::vector<int>, QPolynomial>;
  std::vector<Term> getCoefficients() const override;

  /**
   * Clear all coefficients
   */
  void clear() override;

  /**
   * Export to JSON format
   * @param fileName Output file name
   */
  void exportToJson(const std::string &fileName) const override;

  /**
   * Print polynomial in human-readable format (for debugging)
   * @param maxTerms Maximum number of terms to print
   */
  void print(int maxTerms = 10) const override;

  /**
   * Add another polynomial to this one
   */
  ZMPoly &operator+=(const ZMPoly &other);

  /**
   * Multiply this polynomial by another
   */
  ZMPoly &operator*=(const ZMPoly &other);

  /**
   * Multiply this polynomial by a q-polynomial (QPolynomial)
   * Each q-coefficient is multiplied by the QPolynomial
   */
  ZMPoly operator*(const QPolynomial &qPoly) const;
};
