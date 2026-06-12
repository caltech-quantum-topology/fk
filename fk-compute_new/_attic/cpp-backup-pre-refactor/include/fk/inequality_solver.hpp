#pragma once

#include "fk/polynomial_config.hpp"
#include <vector>
#include <set>

/**
 * Represents an integer solution point
 */
struct IntegerPoint {
    std::vector<int> coordinates;

    IntegerPoint(const std::vector<int>& coords) : coordinates(coords) {}

    bool operator<(const IntegerPoint& other) const {
        return coordinates < other.coordinates;
    }

    bool operator==(const IntegerPoint& other) const {
        return coordinates == other.coordinates;
    }
};

/**
 * Find all integer solutions to a system of linear inequalities
 *
 * @param inequalities Vector of degree-1 multivariable polynomials representing inequalities >= 0
 * @param bounds Search bounds for each variable [min, max] (inclusive)
 * @return Set of integer points that satisfy all inequalities
 *
 * Each polynomial should be of the form: a₁x₁ + a₂x₂ + ... + aₙxₙ + c >= 0
 * where the polynomial coefficients are stored with q^0 power (constant term).
 */
std::set<IntegerPoint> findIntegerSolutions(
    const std::vector<PolynomialType>& inequalities,
    const std::vector<std::pair<int, int>>& bounds
);

/**
 * Evaluate a linear polynomial at a given integer point
 *
 * @param poly The polynomial to evaluate
 * @param point The point at which to evaluate
 * @return The value of the polynomial at the point
 */
int evaluateLinearPolynomial(const PolynomialType& poly, const std::vector<int>& point);

/**
 * Check if a point satisfies all inequalities
 *
 * @param point The point to check
 * @param inequalities The system of inequalities
 * @return true if point satisfies all inequalities, false otherwise
 */
bool satisfiesAllInequalities(const std::vector<int>& point,
                            const std::vector<PolynomialType>& inequalities);