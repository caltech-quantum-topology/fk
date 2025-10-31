#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "fk/bilvector.hpp"

/**
 * MultivariablePolynomial: A polynomial in variables q, x₁, x₂, ..., xₙ
 *
 * Represents polynomials of the form:
 * P(q, x₁, x₂, ..., xₙ) = Σᵢⱼ coeffᵢⱼ × q^j × x₁^a₁ᵢ × x₂^a₂ᵢ × ... × xₙ^aₙᵢ
 *
 * Features:
 * - Arbitrary positive/negative q powers
 * - Configurable number of x variables
 * - Sparse storage (only non-zero coefficients)
 * - Efficient indexing and arithmetic operations
 */
class MultivariablePolynomial {
private:
    int numXVariables;                          // Number of x variables (components)
    int maxDegree;                              // Maximum degree for x variables
    std::vector<int> maxXDegrees;              // Max degree for each x variable
    std::vector<int> blockSizes;               // For converting multi-index to linear index
    std::vector<bilvector<int>> coefficients;  // Storage for all coefficients

    // Convert multi-index (x₁^a₁, x₂^a₂, ...) to linear array index
    int multiIndexToLinear(const std::vector<int>& xPowers) const;

    // Convert linear index back to multi-index
    std::vector<int> linearToMultiIndex(int linearIndex) const;

public:
    /**
     * Constructor
     * @param numVariables Number of x variables
     * @param degree Maximum degree for each x variable
     * @param maxDegrees Optional: different max degree for each variable
     */
    MultivariablePolynomial(int numVariables, int degree,
                           const std::vector<int>& maxDegrees = {});

    /**
     * Get coefficient for specific term
     * @param qPower Power of q
     * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
     * @return Coefficient value
     */
    int getCoefficient(int qPower, const std::vector<int>& xPowers) const;

    /**
     * Set coefficient for specific term
     * @param qPower Power of q
     * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
     * @param coefficient New coefficient value
     */
    void setCoefficient(int qPower, const std::vector<int>& xPowers, int coefficient);

    /**
     * Add to coefficient for specific term
     * @param qPower Power of q
     * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
     * @param coefficient Value to add
     */
    void addToCoefficient(int qPower, const std::vector<int>& xPowers, int coefficient);

    /**
     * Get access to the bilvector for a specific x-multi-index
     * @param xPowers Vector of powers for x₁, x₂, ..., xₙ
     * @return Reference to bilvector containing q-coefficients
     */
    bilvector<int>& getQPolynomial(const std::vector<int>& xPowers);

    const bilvector<int>& getQPolynomial(const std::vector<int>& xPowers) const;

    /**
     * Direct access to underlying storage (for compatibility with existing code)
     */
    std::vector<bilvector<int>>& getCoefficients();

    const std::vector<bilvector<int>>& getCoefficients() const;

    /**
     * Get number of x variables
     */
    int getNumXVariables() const;

    /**
     * Get maximum degrees for each x variable
     */
    const std::vector<int>& getMaxXDegrees() const;

    /**
     * Get block sizes (for advanced usage)
     */
    const std::vector<int>& getBlockSizes() const;

    /**
     * Clear all coefficients
     */
    void clear();

    /**
     * Check if polynomial is zero
     */
    bool isZero() const;

    /**
     * Export to JSON format (compatible with existing output format)
     * @param fileName Output file name
     */
    void exportToJson(const std::string& fileName) const;

    /**
     * Print polynomial in human-readable format (for debugging)
     * @param maxTerms Maximum number of terms to print
     */
    void print(int maxTerms = 10) const;

    /**
     * Arithmetic operations
     */

    /**
     * Add another polynomial to this one
     * @param other Polynomial to add
     * @return Reference to this polynomial after addition
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    MultivariablePolynomial& operator+=(const MultivariablePolynomial& other);

    /**
     * Subtract another polynomial from this one
     * @param other Polynomial to subtract
     * @return Reference to this polynomial after subtraction
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    MultivariablePolynomial& operator-=(const MultivariablePolynomial& other);

    /**
     * Multiply this polynomial by another
     * @param other Polynomial to multiply by
     * @return Reference to this polynomial after multiplication
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    MultivariablePolynomial& operator*=(const MultivariablePolynomial& other);

    /**
     * Add two polynomials
     * @param lhs Left hand side polynomial
     * @param rhs Right hand side polynomial
     * @return New polynomial containing the sum
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    friend MultivariablePolynomial operator+(const MultivariablePolynomial& lhs, const MultivariablePolynomial& rhs);

    /**
     * Subtract two polynomials
     * @param lhs Left hand side polynomial
     * @param rhs Right hand side polynomial
     * @return New polynomial containing the difference
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    friend MultivariablePolynomial operator-(const MultivariablePolynomial& lhs, const MultivariablePolynomial& rhs);

    /**
     * Multiply two polynomials
     * @param lhs Left hand side polynomial
     * @param rhs Right hand side polynomial
     * @return New polynomial containing the product
     * @throws std::invalid_argument if polynomials have incompatible dimensions
     */
    friend MultivariablePolynomial operator*(const MultivariablePolynomial& lhs, const MultivariablePolynomial& rhs);

private:
    /**
     * Check if two polynomials are compatible for arithmetic operations
     * @param other Polynomial to check compatibility with
     * @throws std::invalid_argument if polynomials are incompatible
     */
    void checkCompatibility(const MultivariablePolynomial& other) const;
};
