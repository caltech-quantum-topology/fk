#include <iostream>
#include "fk/multivariable_polynomial.hpp"

int main() {
    // Create a polynomial in q, x1, x2 with max degree 3 for each x variable
    MultivariablePolynomial poly(2, 3);

    std::cout << "=== MultivariablePolynomial Demo ===\n\n";

    // Set some coefficients using the intuitive interface
    // P(q, x1, x2) = 5*q^2*x1^1*x2^2 + 3*q^(-1)*x1^0*x2^1 + 7*q^3*x1^2*x2^0

    poly.setCoefficient(2, {1, 2}, 5);   // 5*q^2*x1^1*x2^2
    poly.setCoefficient(-1, {0, 1}, 3);  // 3*q^(-1)*x1^0*x2^1
    poly.setCoefficient(3, {2, 0}, 7);   // 7*q^3*x1^2*x2^0

    // Add to existing coefficients
    poly.addToCoefficient(2, {1, 2}, 2); // Now coefficient is 7

    std::cout << "Created polynomial with coefficients:\n";
    std::cout << "- Coefficient of q^2*x1^1*x2^2: " << poly.getCoefficient(2, {1, 2}) << "\n";
    std::cout << "- Coefficient of q^(-1)*x1^0*x2^1: " << poly.getCoefficient(-1, {0, 1}) << "\n";
    std::cout << "- Coefficient of q^3*x1^2*x2^0: " << poly.getCoefficient(3, {2, 0}) << "\n\n";

    // Print the polynomial
    std::cout << "Polynomial structure:\n";
    poly.print(10);
    std::cout << "\n";

    // Access the q-polynomial for a specific x-term
    std::cout << "Q-polynomial for x1^1*x2^2:\n";
    auto& qPoly = poly.getQPolynomial({1, 2});
    for (int j = qPoly.getMaxNegativeIndex(); j <= qPoly.getMaxPositiveIndex(); j++) {
        if (qPoly[j] != 0) {
            std::cout << "  q^" << j << ": " << qPoly[j] << "\n";
        }
    }
    std::cout << "\n";

    // Export to JSON
    poly.exportToJson("example_polynomial");
    std::cout << "Exported to example_polynomial.json\n\n";

    // Show metadata
    std::cout << "Polynomial metadata:\n";
    std::cout << "- Number of x variables: " << poly.getNumXVariables() << "\n";
    std::cout << "- Max degrees: [";
    auto maxDegrees = poly.getMaxXDegrees();
    for (int i = 0; i < maxDegrees.size(); i++) {
        std::cout << maxDegrees[i];
        if (i < maxDegrees.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";

    return 0;
}
