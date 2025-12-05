#pragma once

#include "fk/polynomial_base.hpp"
#include "fk/qpolynomial.hpp"
#include <flint/fmpz_mpoly.h>
#include <string>
#include <vector>

/**
 * FMPoly: A FLINT-based multivariate polynomial class
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - FLINT backend for optimal performance
 * - Dense polynomial representation
 * - Support for negative q exponents via offset tracking
 * - Drop-in replacement for MultivariablePolynomial
 */
class FMPoly : public PolynomialBase<int, QPolynomial> {
private:
  int numXVariables; // Number of x variables
  std::vector<int>
      allGroundPowers;  // Ground powers: [q_min, x1_min, x2_min, ..., xn_min]
  fmpz_mpoly_ctx_t ctx; // FLINT context for multivariate polynomials
  fmpz_mpoly_t poly;    // Main FLINT polynomial

  std::vector<int> maxXDegrees; // Max degrees (for compatibility)
  std::vector<int> blockSizes;  // Block sizes (for compatibility)

  // Helper methods
  void setupContext();
  void convertExponents(int qPower, const std::vector<int> &xPowers,
                        fmpz **exps, slong *exp_bits) const;
  bool getExponentsFromMonomial(const fmpz *exps, int &qPower,
                                std::vector<int> &xPowers) const;
  void adjustGroundPowersIfNeeded(int qPower, const std::vector<int> &xPowers);

  // Internal methods used by public interface
  int getCoefficient(int qPower, const std::vector<int> &xPowers) const;
  QPolynomial getQPolynomialObject(const std::vector<int> &xPowers) const;
  void setQPolynomial(const std::vector<int> &xPowers,
                      const std::vector<int> &qCoeffs, int minQPower = 0);
  void setQPolynomial(const std::vector<int> &xPowers,
                      const QPolynomial &qPoly);
  void addQPolynomial(const std::vector<int> &xPowers,
                      const QPolynomial &qPoly);
  void checkCompatibility(const FMPoly &other) const;
  bool isZero() const;

public:
  using Term = std::pair<std::vector<int>, QPolynomial>;
  std::vector<Term> getCoefficients() const override;
  /**
   * Constructor
   * @param numVariables Number of x variables
   * @param degree Maximum degree for each x variable (advisory)
   * @param maxDegrees Optional: different max degree for each variable
   */
  FMPoly(int numVariables, int degree = 10,
         const std::vector<int> &maxDegrees = {});

  /**
   * Constructor to increase the number of variables from another polynomial
   * @param source The source polynomial with fewer variables
   * @param newNumVariables The new number of variables (must be >= source's numVariables)
   * @param targetVariableIndex The index (0-based) where the source's variable should be mapped
   * @param degree Maximum degree for each new x variable (advisory only)
   * @param maxDegrees Optional: different max degree for each variable (advisory only)
   */
  FMPoly(const FMPoly &source, int newNumVariables, int targetVariableIndex,
         int degree = 10, const std::vector<int> &maxDegrees = {});

  /**
   * Copy constructor
   */
  FMPoly(const FMPoly &other);

  /**
   * Assignment operator
   */
  FMPoly &operator=(const FMPoly &other);

  /**
   * Destructor
   */
  ~FMPoly();

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
   * Truncate polynomial to given degrees
   * @param maxXdegrees Maximum degree for each x variable
   */
  FMPoly truncate(const std::vector<int> &maxXdegrees) const;

  /**
   * Truncate polynomial to same degree for all x variables
   * @param maxDegree Maximum degree to keep for all x variables
   */
  FMPoly truncate(int maxDegree) const;

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
   * Add another polynomial to this one
   */
  FMPoly &operator+=(const FMPoly &other);

  /**
   * Multiply this polynomial by another
   */
  FMPoly &operator*=(const FMPoly &other);

  /**
   * Multiply this polynomial by a q-polynomial (QPolynomial)
   * Each coefficient is multiplied by the QPolynomial
   */
  FMPoly operator*(const QPolynomial &qPoly) const;

  /**
   * Get number of x variables
   * (Kept for compatibility with existing code)
   */
  int getNumXVariables() const;

  /**
   * Get access to the underlying FLINT polynomial (for advanced operations)
   * Similar to MultivariablePolynomial::getCoefficientMap()
   */
  const fmpz_mpoly_t &getFlintPoly() const { return poly; }
  const fmpz_mpoly_ctx_t &getFlintContext() const { return ctx; }
};
