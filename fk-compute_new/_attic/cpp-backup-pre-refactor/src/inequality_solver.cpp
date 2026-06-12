#include "fk/inequality_solver.hpp"
#include <iostream>
#include <algorithm>
#include <functional>

int evaluateLinearPolynomial(const PolynomialType& poly, const std::vector<int>& point) {
    int result = 0;
    const auto coeffMap = poly.getCoefficients();

    for (const auto& [xPowers, qPoly] : coeffMap) {
        if (xPowers.size() != point.size()) {
            throw std::invalid_argument("Point dimension doesn't match polynomial variables");
        }

        // Calculate x₁^a₁ × x₂^a₂ × ... × xₙ^aₙ
        int termValue = 1;
        for (size_t i = 0; i < xPowers.size(); ++i) {
            if (xPowers[i] > 1) {
                throw std::invalid_argument("Polynomial is not linear (degree > 1)");
            }
            if (xPowers[i] < 0) {
                throw std::invalid_argument("Negative exponents not supported for integer evaluation");
            }
            if (xPowers[i] == 1) {
                termValue *= point[i];
            }
        }

        // Get the coefficient for q^0 (constant coefficient in q)
        int coefficient = qPoly[0];
        result += coefficient * termValue;
    }

    return result;
}

bool satisfiesAllInequalities(const std::vector<int>& point,
                            const std::vector<PolynomialType>& inequalities) {
    for (const auto& inequality : inequalities) {
        int value = evaluateLinearPolynomial(inequality, point);
        if (value < 0) {
            return false;
        }
    }
    return true;
}

std::set<IntegerPoint> findIntegerSolutions(
    const std::vector<PolynomialType>& inequalities,
    const std::vector<std::pair<int, int>>& bounds) {

    std::set<IntegerPoint> solutions;

    if (inequalities.empty() || bounds.empty()) {
        return solutions;
    }

    int numVars = bounds.size();

    // Validate that all polynomials have the same number of variables
    for (const auto& poly : inequalities) {
        if (poly.getNumXVariables() != numVars) {
            throw std::invalid_argument("All polynomials must have the same number of variables");
        }
    }

    // Generate all possible integer points within bounds
    std::vector<int> current(numVars);
    std::function<void(int)> generatePoints = [&](int varIndex) {
        if (varIndex == numVars) {
            // Check if current point satisfies all inequalities
            if (satisfiesAllInequalities(current, inequalities)) {
                solutions.insert(IntegerPoint(current));
            }
            return;
        }

        // Try all values for current variable within bounds
        for (int val = bounds[varIndex].first; val <= bounds[varIndex].second; ++val) {
            current[varIndex] = val;
            generatePoints(varIndex + 1);
        }
    };

    generatePoints(0);
    return solutions;
}