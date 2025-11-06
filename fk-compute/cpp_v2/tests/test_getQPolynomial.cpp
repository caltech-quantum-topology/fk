#include "fk/fmpoly.hpp"
#include <iostream>
#include <vector>

void testGetQPolynomial() {
    std::cout << "=== Testing FMPoly::getQPolynomial() ===\n\n";

    try {
        // Create a polynomial in q, x1, x2
        FMPoly poly(2, 5);

        // Set some coefficients for a specific x-monomial
        // For x1^1, x2^2: set coefficients for different q powers
        std::vector<int> xPowers = {1, 2};
        poly.setCoefficient(-2, xPowers, 3);  // 3*q^(-2)*x1^1*x2^2
        poly.setCoefficient(0, xPowers, 5);   // 5*q^0*x1^1*x2^2
        poly.setCoefficient(2, xPowers, 7);   // 7*q^2*x1^1*x2^2
        poly.setCoefficient(4, xPowers, 2);   // 2*q^4*x1^1*x2^2

        // Add some coefficients for different x-monomials to verify filtering
        poly.setCoefficient(1, {0, 0}, 10);   // 10*q^1*x1^0*x2^0
        poly.setCoefficient(3, {2, 1}, 4);    // 4*q^3*x1^2*x2^1

        std::cout << "Set up polynomial with:\n";
        std::cout << "- 3*q^(-2)*x1^1*x2^2\n";
        std::cout << "- 5*q^0*x1^1*x2^2\n";
        std::cout << "- 7*q^2*x1^1*x2^2\n";
        std::cout << "- 2*q^4*x1^1*x2^2\n";
        std::cout << "- 10*q^1*x1^0*x2^0 (different x-monomial)\n";
        std::cout << "- 4*q^3*x1^2*x2^1 (different x-monomial)\n\n";

        // Test getQPolynomial for x1^1*x2^2
        std::cout << "Getting Q-polynomial for x1^1*x2^2:\n";
        std::vector<int> qPoly = poly.getQPolynomial(xPowers);

        std::cout << "Returned vector size: " << qPoly.size() << "\n";
        std::cout << "Q-polynomial coefficients: [";
        for (size_t i = 0; i < qPoly.size(); i++) {
            std::cout << qPoly[i];
            if (i < qPoly.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";

        // The vector should contain coefficients for q^(-2) through q^4
        // Since it returns a vector, the interpretation depends on the minimum q-power
        std::cout << "\nExpected: coefficients for q powers from -2 to 4\n";
        std::cout << "Should be: [3, 0, 5, 0, 7, 0, 2] (for q^(-2), q^(-1), q^0, q^1, q^2, q^3, q^4)\n\n";

        // Test getQPolynomial for a different x-monomial
        std::vector<int> xPowers2 = {0, 0};
        std::cout << "Getting Q-polynomial for x1^0*x2^0:\n";
        std::vector<int> qPoly2 = poly.getQPolynomial(xPowers2);

        std::cout << "Returned vector size: " << qPoly2.size() << "\n";
        std::cout << "Q-polynomial coefficients: [";
        for (size_t i = 0; i < qPoly2.size(); i++) {
            std::cout << qPoly2[i];
            if (i < qPoly2.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
        std::cout << "Expected: [10] (for q^1)\n\n";

        // Test getQPolynomial for non-existent x-monomial
        std::vector<int> xPowers3 = {3, 3};
        std::cout << "Getting Q-polynomial for non-existent x1^3*x2^3:\n";
        std::vector<int> qPoly3 = poly.getQPolynomial(xPowers3);

        std::cout << "Returned vector size: " << qPoly3.size() << "\n";
        std::cout << "Expected: empty vector (size 0)\n\n";

        // Verify individual coefficients still work
        std::cout << "Verification using getCoefficient:\n";
        std::cout << "Coefficient of q^(-2)*x1^1*x2^2: " << poly.getCoefficient(-2, {1, 2}) << " (expected: 3)\n";
        std::cout << "Coefficient of q^0*x1^1*x2^2: " << poly.getCoefficient(0, {1, 2}) << " (expected: 5)\n";
        std::cout << "Coefficient of q^2*x1^1*x2^2: " << poly.getCoefficient(2, {1, 2}) << " (expected: 7)\n";
        std::cout << "Coefficient of q^4*x1^1*x2^2: " << poly.getCoefficient(4, {1, 2}) << " (expected: 2)\n";
        std::cout << "Coefficient of q^1*x1^0*x2^0: " << poly.getCoefficient(1, {0, 0}) << " (expected: 10)\n";

        std::cout << "\n✓ getQPolynomial test completed!\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    testGetQPolynomial();
    return 0;
}