#pragma once

#include "fk/bilvector.hpp"
#include "fk/polynomial_base.hpp"
#include "fk/qpolynomial.hpp"
#include <flint/fmpz_mpoly.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * Custom hash function for vector<int> keys
 */
struct HMPolyVectorHash {
  std::size_t operator()(const std::vector<int> &v) const;
};

/**
 * HMPoly: Hybrid Multi Poly - Best of Both Worlds
 *
 * Combines the strengths of both implementations:
 * - Uses unordered_map for fast individual coefficient access (O(1))
 * - Uses FLINT for bulk operations like multiplication (optimized)
 *
 * Internal representation: unordered_map (like MultivariablePolynomial)
 * - Fast addToCoefficient, setCoefficient, operator+=
 * - Private conversion to/from FLINT for bulk operations
 *
 * Drop-in replacement for both MultivariablePolynomial and FMPoly
 */
class HMPoly : public PolynomialBase<int, bilvector<int>> {
private:
  int numXVariables;            // Number of x variables
  std::vector<int> maxXDegrees; // Max degrees (advisory only)
  std::vector<int> blockSizes;  // Block sizes (for compatibility)

  // Primary storage: sparse map (fast individual access)
  std::unordered_map<std::vector<int>, bilvector<int>, HMPolyVectorHash>
      coeffs_;

  // FLINT context (created on-demand for bulk operations)
  mutable fmpz_mpoly_ctx_t ctx;
  mutable bool ctx_initialized;

  // Helper methods for FLINT conversion
  void initFlintContext() const;
  void clearFlintContext() const;
  void toFlint(fmpz_mpoly_t poly) const;
  void fromFlint(const fmpz_mpoly_t poly);

  // Helper methods
  void pruneZeros();

public:
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable
   */
  HMPoly(int numVariables, int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Constructor to increase the number of variables from another polynomial
   * @param source The source polynomial with fewer variables
   * @param newNumVariables The new number of variables (must be >= source's
   * numVariables)
   * @param targetVariableIndex The index (0-based) where the source's variable
   * should be mapped
   * @param degree Maximum degree for each new x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable
   * (advisory only)
   */
  HMPoly(const HMPoly &source, int newNumVariables, int targetVariableIndex,
         int degree = 10, const std::vector<int> &maxDegrees = {});

  /**
   * Copy constructor
   */
  HMPoly(const HMPoly &other);

  /**
   * Assignment operator
   */
  HMPoly &operator=(const HMPoly &other);

  /**
   * Destructor
   */
  ~HMPoly();

  /**
   * Set coefficient for specific term (uses map - O(1))
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient New coefficient value
   */
  void setCoefficient(int qPower, const std::vector<int> &xPowers,
                      int coefficient) override;

  /**
   * Add to coefficient for specific term (uses map - O(1))
   * @param qPower Power of q
   * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
   * @param coefficient Value to add
   */
  void addToCoefficient(int qPower, const std::vector<int> &xPowers,
                        int coefficient) override;

  /**
   * Get coefficients as a vector of (x-powers, q-polynomial) pairs
   */
  using Term = std::pair<std::vector<int>, bilvector<int>>;
  std::vector<Term> getCoefficients() const override;

  /**
   * Truncate polynomial to given degrees
   * @param maxXdegrees Maximum degree for each x variable
   */
  HMPoly truncate(const std::vector<int> &maxXdegrees) const;

  /**
   * Truncate polynomial to same degree for all x variables
   * @param maxDegree Maximum degree to keep for all x variables
   */
  HMPoly truncate(int maxDegree) const;

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
   * Print polynomial in human-readable format
   * @param maxTerms Maximum number of terms to print
   */
  void print(int maxTerms = 10) const override;

  /**
   * Add another polynomial to this one (uses map - fast)
   */
  HMPoly &operator+=(const HMPoly &other);

  /**
   * Multiply this polynomial by another (uses FLINT - optimized for bulk)
   */
  HMPoly &operator*=(const HMPoly &other);

  /**
   * Multiply this polynomial by a q-polynomial (bilvector)
   * Each q-coefficient is multiplied by the bilvector
   */
  HMPoly operator*(const bilvector<int> &qPoly) const;

  /**
   * Get number of x variables
   */
  int getNumXVariables() const { return numXVariables; }

  /**
   * Get access to the internal coefficient map (for advanced operations)
   */
  const std::unordered_map<std::vector<int>, bilvector<int>,
                           HMPolyVectorHash> &
  getCoefficientMap() const {
    return coeffs_;
  }
};
