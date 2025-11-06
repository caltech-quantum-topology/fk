#pragma once

#include <string>
#include <vector>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_poly.h>

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
    fmpz_poly_t poly;        // FLINT univariate polynomial
    int minPower;            // Minimum power of q (for negative exponents)

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
    std::vector<int> allGroundPowers;     // Ground powers: [q_min, x1_min, x2_min, ..., xn_min]
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
    void adjustGroundPowersIfNeeded(int qPower, const std::vector<int> &xPowers);

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
     * Get a q-polynomial as QPolynomial object for a specific x-multi-index
     * Returns a QPolynomial object with proper handling of negative powers
     */
    QPolynomial getQPolynomialObject(const std::vector<int> &xPowers) const;

    /**
     * Set q-polynomial for a specific x-multi-index
     */
    void setQPolynomial(const std::vector<int> &xPowers,
                       const std::vector<int> &qCoeffs, int minQPower = 0);

    /**
     * Set q-polynomial for a specific x-multi-index using QPolynomial object
     */
    void setQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly);

    /**
     * Add a QPolynomial to existing q-polynomial for a specific x-multi-index
     */
    void addQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly);

    /**
     * Multiply existing q-polynomial by a QPolynomial for a specific x-multi-index
     */
    void multiplyQPolynomial(const std::vector<int> &xPowers, const QPolynomial &qPoly);

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
