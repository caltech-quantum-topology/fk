#pragma once

#include "fk/bilvector.hpp"
#include "fk/polynomial_base.hpp"
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * Custom hash function for vector<int> keys
 */
struct VectorHash {
  std::size_t operator()(const std::vector<int> &v) const;
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
class MultivariablePolynomial : public PolynomialBase<int> {
private:
  int numXVariables;            // Number of x variables (components)
  std::vector<int> maxXDegrees; // Max degrees (advisory only, not enforced)
  std::vector<int> blockSizes; // Block sizes (advisory only, for compatibility)

  // Prune zero coefficients from the map
  void pruneZeros();

public:
  // Choose a reasonable upper bound for number of x variables
  static constexpr int kMaxXVariables = 5; // adjust if you need more

  struct ExponentKey {
    std::array<int, kMaxXVariables> e{};

    bool operator==(const ExponentKey &other) const noexcept {
      return e == other.e;
    }
    
    // For std::sort and other ordered operations
    bool operator<(const ExponentKey &other) const noexcept {
      return e < other.e;  // std::array has lexicographic operator<
    }
  };

  ExponentKey makeKey(const std::vector<int> &xPowers) const;

  struct ExponentKeyHash {
    std::size_t operator()(const ExponentKey &ex) const noexcept {
      std::size_t h = 0;
      // We will only use first numXVariables entries in practice,
      // but hashing all kMaxXVariables is fine (or we can restrict later).
      for (int v : ex.e) {
        h = h * 1315423911u + std::hash<int>{}(v);
      }
      return h;
    }
  };

  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable
   * (advisory only)
   */
  MultivariablePolynomial(int numVariables, int degree = 10,
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
  MultivariablePolynomial(const MultivariablePolynomial &source,
                          int newNumVariables, int targetVariableIndex,
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
  MultivariablePolynomial truncate(const std::vector<int> &maxXdegrees) const;

  /**
   * Truncate multivariable polynomial to the same degree for all x variables
   * @param maxDegree Maximum degree to keep for all x variables
   */
  MultivariablePolynomial truncate(int maxDegree) const;

  /**
   * Get coefficients as a vector of (x-powers, q-polynomial) pairs
   * Compatible with FMPoly's interface
   */
  using Term = std::pair<std::vector<int>, bilvector<int>>;
  std::vector<Term> getCoefficients() const override;

  /**
   * Get number of x variables
   */
  int getNumXVariables() const { return numXVariables; }

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
  MultivariablePolynomial &operator+=(const MultivariablePolynomial &other);

  /**
   * Multiply this polynomial by another
   */
  MultivariablePolynomial &operator*=(const MultivariablePolynomial &other);

  /**
   * Multiply this polynomial by a q-polynomial (bilvector)
   * Each q-coefficient is multiplied by the bilvector
   */
  MultivariablePolynomial operator*(const bilvector<int> &qPoly) const;

private:
  // Sparse storage: map from x-exponent vector to q-polynomial
  std::unordered_map<ExponentKey, bilvector<int>, ExponentKeyHash> coeffs_;
};
