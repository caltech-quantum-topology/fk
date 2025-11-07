#pragma once

#include "fk/bilvector.hpp"
#include <string>
#include <vector>

/**
 * BMPoly: Basic Multivariable Polynomial
 *
 * A polynomial in variables q, x₁, x₂, ..., xₙ with vector-based internal representation.
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - Dense vector-based storage using std::vector<bilvector<int>>
 * - Support for negative x exponents via ground degree tracking
 * - Same public API as MultivariablePolynomial for drop-in replacement
 * - Efficient indexing using multi-index to linear index conversion
 */
class BMPoly {
private:
  int numXVariables;                      // Number of x variables (components)
  std::vector<int> maxXDegrees;          // Max degrees for each x variable
  std::vector<int> groundXDegrees;       // Ground degrees (min exponents) for negative support
  std::vector<int> blockSizes;           // Block sizes for multi-index conversion

  // Dense storage: vector where index represents x-multi-index, value is q-polynomial
  std::vector<bilvector<int>> coeffs_;

  // Helper methods for index conversion and bounds management
  int multiIndexToLinear(const std::vector<int> &xPowers) const;
  std::vector<int> linearToMultiIndex(int linearIndex) const;
  void expandStorageIfNeeded(const std::vector<int> &xPowers);
  void updateGroundDegrees(const std::vector<int> &xPowers);
  void recalculateBlockSizes();
  void resizeStorage();

  // Internal helpers
  void pruneZeros();
  void checkIndexBounds(const std::vector<int> &xPowers) const;

public:
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable
   */
  BMPoly(int numVariables, int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Constructor to increase the number of variables from another polynomial
   * @param source The source polynomial with fewer variables
   * @param newNumVariables The new number of variables (must be >= source's numVariables)
   * @param targetVariableIndex The index (0-based) where the source's variable should be mapped
   * @param degree Maximum degree for each new x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable (advisory only)
   */
  BMPoly(const BMPoly &source,
         int newNumVariables,
         int targetVariableIndex,
         int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Get coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @return Coefficient value
   */
  int getCoefficient(int qPower, const std::vector<int> &xPowers) const;

  /**
   * Set coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient New coefficient value
   */
  void setCoefficient(int qPower, const std::vector<int> &xPowers,
                      int coefficient);

  /**
   * Add to coefficient for specific term
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient Value to add
   */
  void addToCoefficient(int qPower, const std::vector<int> &xPowers,
                        int coefficient);

  /**
   * Get access to the bilvector for a specific x-multi-index
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @return Reference to bilvector containing q-coefficients
   */
  bilvector<int> &getQPolynomial(const std::vector<int> &xPowers);

  const bilvector<int> &getQPolynomial(const std::vector<int> &xPowers) const;

  /**
   * Invert variable at target_index
   */
  BMPoly invertVariable(const int target_index);

  /**
   * Truncate multivariable polynomial to given degrees
   */
  BMPoly truncate(const std::vector<int> &maxXdegrees);

  /**
   * Get coefficients as sparse representation (same format as MultivariablePolynomial)
   * Returns vector of pairs: (x-powers, q-polynomial) for non-zero terms only
   */
  using Term = std::pair<std::vector<int>, bilvector<int>>;
  const std::vector<Term> getCoefficients() const;

  /**
   * Get coefficients as dense vector (original BMPoly format)
   * Returns the internal dense vector representation
   */
  std::vector<bilvector<int>> &getCoefficientsDense();

  /**
   * Sync sparse vector back to internal representation
   * This must be called after any operations that modify the sparse vector
   */
  void syncFromSparseVector(const std::vector<std::pair<std::vector<int>, bilvector<int>>> &sparseVector);

  /**
   * Sync dense vector back to internal representation
   * This must be called after any operations that modify the dense vector
   */
  void syncFromDenseVector(const std::vector<bilvector<int>> &denseVector);

  /**
   * Get number of x variables
   */
  int getNumXVariables() const;

  /**
   * Get maximum degrees for each x variable (advisory only)
   */
  const std::vector<int> &getMaxXDegrees() const;

  /**
   * Get block sizes (for compatibility, not used internally)
   */
  const std::vector<int> &getBlockSizes() const;

  /**
   * Clear all coefficients
   */
  void clear();

  /**
   * Check if polynomial is zero
   */
  bool isZero() const;

  /**
   * Export to JSON format
   * @param fileName Output file name
   */
  void exportToJson(const std::string &fileName) const;

  /**
   * Print polynomial in human-readable format (for debugging)
   * @param maxTerms Maximum number of terms to print
   */
  void print(int maxTerms = 10) const;

  /**
   * Evaluate polynomial at a given point, returning a polynomial in q
   * @param point Vector of values for x₁, x₂, ..., xₙ (must match numXVariables)
   * @return bilvector<int> representing the resulting polynomial in q
   * @throws std::invalid_argument if point dimension doesn't match numXVariables
   * @throws std::domain_error if division by zero occurs (negative exponent with zero value)
   */
  bilvector<int> evaluate(const std::vector<int> &point) const;

  /**
   * Check if two polynomials are compatible for arithmetic operations
   */
  void checkCompatibility(const BMPoly &other) const;

  /**
   * Add another polynomial to this one
   */
  BMPoly &operator+=(const BMPoly &other);

  /**
   * Subtract another polynomial from this one
   */
  BMPoly &operator-=(const BMPoly &other);

  /**
   * Multiply this polynomial by another
   */
  BMPoly &operator*=(const BMPoly &other);

  /**
   * Add a pure q-polynomial (bilvector) to this polynomial.
   * The q-polynomial is attached to the monomial with all x-exponents = 0.
   */
  BMPoly &operator+=(const bilvector<int> &qPoly);

  /**
   * Subtract a pure q-polynomial (bilvector) from this polynomial.
   * The q-polynomial is attached to the monomial with all x-exponents = 0.
   */
  BMPoly &operator-=(const bilvector<int> &qPoly);

  /**
   * Multiply this polynomial by a pure q-polynomial (bilvector).
   * Each q-polynomial coefficient (for each x-monomial) is multiplied by qPoly.
   */
  BMPoly &operator*=(const bilvector<int> &qPoly);

  /**
   * Friend functions for binary operations
   */
  friend BMPoly operator+(const BMPoly &lhs, const BMPoly &rhs);
  friend BMPoly operator-(const BMPoly &lhs, const BMPoly &rhs);
  friend BMPoly operator*(const BMPoly &lhs, const BMPoly &rhs);

  friend BMPoly operator+(const BMPoly &lhs, const bilvector<int> &rhs);
  friend BMPoly operator+(const bilvector<int> &lhs, const BMPoly &rhs);

  friend BMPoly operator-(const BMPoly &lhs, const bilvector<int> &rhs);
  friend BMPoly operator-(const bilvector<int> &lhs, const BMPoly &rhs);

  friend BMPoly operator*(const BMPoly &lhs, const bilvector<int> &rhs);
  friend BMPoly operator*(const bilvector<int> &lhs, const BMPoly &rhs);
};
