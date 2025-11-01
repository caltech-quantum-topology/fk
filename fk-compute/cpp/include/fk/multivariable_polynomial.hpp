#pragma once

#include "fk/bilvector.hpp"
#include <string>
#include <unordered_map>
#include <vector>

/**
 * Custom hash function for vector<int> keys
 */
struct VectorHash {
  std::size_t operator()(const std::vector<int>& v) const;
};

/**
 * MultivariablePolynomial: A polynomial in variables q, x₁, x₂, ..., xₙ
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - Arbitrary positive/negative q powers
 * - Configurable number of x variables
 * - Sparse storage using unordered_map (only non-zero coefficients)
 * - Support for negative x exponents
 * - Efficient indexing and arithmetic operations
 */
class MultivariablePolynomial {
private:
  int numXVariables;            // Number of x variables (components)
  std::vector<int> maxXDegrees; // Max degrees (advisory only, not enforced)
  std::vector<int> blockSizes;  // Block sizes (advisory only, for compatibility)

  // Sparse storage: map from x-exponent vector to q-polynomial
  std::unordered_map<std::vector<int>, bilvector<int>, VectorHash> coeffs_;

  // Prune zero coefficients from the map
  void pruneZeros();

public:
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable (advisory only)
   */
  MultivariablePolynomial(int numVariables, int degree = 10,
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
   * Get read-only access to coefficient map (replaces getCoefficients)
   */
  const std::unordered_map<std::vector<int>, bilvector<int>, VectorHash>&
  getCoefficientMap() const;

  /**
   * Backward compatibility: Convert sparse representation to dense vector
   * This is needed for existing code that expects vector<bilvector<int>>
   */
  std::vector<bilvector<int>> getCoefficients() const;

  /**
   * Backward compatibility: Non-const version
   * WARNING: This creates a mutable copy that must be manually synced back
   */
  std::vector<bilvector<int>>& getCoefficients();

  /**
   * Sync dense vector back to sparse representation
   * This must be called after any operations that modify the dense vector
   */
  void syncFromDenseVector(const std::vector<bilvector<int>>& denseVector);

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
   * Export to JSON format with new sparse format
   * @param fileName Output file name
   */
  void exportToJson(const std::string &fileName) const;

  /**
   * Print polynomial in human-readable format (for debugging)
   * @param maxTerms Maximum number of terms to print
   */
  void print(int maxTerms = 10) const;

  /**
   * Check if two polynomials are compatible for arithmetic operations
   */
  void checkCompatibility(const MultivariablePolynomial &other) const;

  /**
   * Add another polynomial to this one
   */
  MultivariablePolynomial &operator+=(const MultivariablePolynomial &other);

  /**
   * Subtract another polynomial from this one
   */
  MultivariablePolynomial &operator-=(const MultivariablePolynomial &other);

  /**
   * Multiply this polynomial by another
   */
  MultivariablePolynomial &operator*=(const MultivariablePolynomial &other);

  /**
   * Friend functions for binary operations
   */
  friend MultivariablePolynomial operator+(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs);

  friend MultivariablePolynomial operator-(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs);

  friend MultivariablePolynomial operator*(const MultivariablePolynomial &lhs,
                                           const MultivariablePolynomial &rhs);
};