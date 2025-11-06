#pragma once

#include <string>
#include <vector>
#include <flint/fmpz_mpoly.h>

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
class FMPoly {
private:
    int numXVariables;                    // Number of x variables
    int qOffset;                          // Offset for handling negative q powers
    fmpz_mpoly_ctx_t ctx;                 // FLINT context for multivariate polynomials
    fmpz_mpoly_t poly;                    // Main FLINT polynomial

    std::vector<int> maxXDegrees;         // Max degrees (for compatibility)
    std::vector<int> blockSizes;          // Block sizes (for compatibility)

    // Helper methods
    void setupContext();
    void convertExponents(int qPower, const std::vector<int> &xPowers,
                         fmpz **exps, slong *exp_bits) const;
    bool getExponentsFromMonomial(const fmpz *exps, int &qPower,
                                 std::vector<int> &xPowers) const;

public:
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
     */
    FMPoly(const FMPoly &source, int newNumVariables,
           int targetVariableIndex, int degree = 10,
           const std::vector<int> &maxDegrees = {});

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
     * Get a q-polynomial (univariate in q) for a specific x-multi-index
     * Returns coefficients as a vector where index i represents q^(i + minQPower)
     */
    std::vector<int> getQPolynomial(const std::vector<int> &xPowers) const;

    /**
     * Set q-polynomial for a specific x-multi-index
     */
    void setQPolynomial(const std::vector<int> &xPowers,
                       const std::vector<int> &qCoeffs, int minQPower = 0);

    /**
     * Invert a specific variable
     */
    FMPoly invertVariable(const int target_index) const;

    /**
     * Truncate polynomial to given degrees
     */
    FMPoly truncate(const std::vector<int> &maxXdegrees) const;

    /**
     * Get number of x variables
     */
    int getNumXVariables() const;

    /**
     * Get maximum degrees for each x variable
     */
    const std::vector<int> &getMaxXDegrees() const;

    /**
     * Get block sizes (for compatibility)
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
     * Print polynomial in human-readable format
     * @param maxTerms Maximum number of terms to print
     */
    void print(int maxTerms = 10) const;

    /**
     * Evaluate polynomial at a given point, returning coefficients for q
     * @param point Vector of values for x₁, x₂, ..., xₙ
     * @return Vector representing coefficients of the resulting polynomial in q
     */
    std::vector<int> evaluate(const std::vector<int> &point) const;

    /**
     * Check if two polynomials are compatible for arithmetic operations
     */
    void checkCompatibility(const FMPoly &other) const;

    /**
     * Arithmetic operators
     */
    FMPoly &operator+=(const FMPoly &other);
    FMPoly &operator-=(const FMPoly &other);
    FMPoly &operator*=(const FMPoly &other);

    /**
     * Friend functions for binary operations
     */
    friend FMPoly operator+(const FMPoly &lhs, const FMPoly &rhs);
    friend FMPoly operator-(const FMPoly &lhs, const FMPoly &rhs);
    friend FMPoly operator*(const FMPoly &lhs, const FMPoly &rhs);

    /**
     * Get access to the underlying FLINT polynomial (for advanced operations)
     */
    const fmpz_mpoly_t &getFlintPoly() const { return poly; }
    const fmpz_mpoly_ctx_t &getFlintContext() const { return ctx; }
};