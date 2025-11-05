#include "fk/fmpoly.hpp"
#include <iostream>

int main() {
    std::cout << "=== FMPoly (FLINT-based) Demo ===\n\n";

    try {
        // Create a polynomial in q, x1, x2 with max degree 3 for each x variable
        FMPoly poly(2, 3);

        std::cout << "Created FMPoly with 2 x-variables and max degree 3\n";

        // Set some coefficients using the same interface as MultivariablePolynomial
        // P(q, x1, x2) = 5*q^2*x1^1*x2^2 + 3*q^(-1)*x1^0*x2^1 + 7*q^3*x1^2*x2^0
        poly.setCoefficient(2, {1, 2}, 5);  // 5*q^2*x1^1*x2^2
        poly.setCoefficient(-1, {0, 1}, 3); // 3*q^(-1)*x1^0*x2^1
        poly.setCoefficient(3, {2, 0}, 7);  // 7*q^3*x1^2*x2^0

        std::cout << "Set coefficients:\n";
        std::cout << "- Coefficient of q^2*x1^1*x2^2: "
                  << poly.getCoefficient(2, {1, 2}) << "\n";
        std::cout << "- Coefficient of q^(-1)*x1^0*x2^1: "
                  << poly.getCoefficient(-1, {0, 1}) << "\n";
        std::cout << "- Coefficient of q^3*x1^2*x2^0: "
                  << poly.getCoefficient(3, {2, 0}) << "\n\n";

        // Add to existing coefficients
        poly.addToCoefficient(2, {1, 2}, 2); // Now coefficient should be 7
        std::cout << "After adding 2 to q^2*x1^1*x2^2 coefficient: "
                  << poly.getCoefficient(2, {1, 2}) << "\n\n";

        // Print the polynomial structure
        std::cout << "Polynomial structure:\n";
        poly.print(10);
        std::cout << "\n";

        // Test arithmetic operations
        FMPoly poly2(2, 3);
        poly2.setCoefficient(1, {1, 1}, 4); // 4*q^1*x1^1*x2^1
        poly2.setCoefficient(2, {0, 0}, 2); // 2*q^2*x1^0*x2^0

        std::cout << "Created second polynomial with:\n";
        std::cout << "- Coefficient of q^1*x1^1*x2^1: "
                  << poly2.getCoefficient(1, {1, 1}) << "\n";
        std::cout << "- Coefficient of q^2*x1^0*x2^0: "
                  << poly2.getCoefficient(2, {0, 0}) << "\n\n";

        // Test addition
        FMPoly sum = poly + poly2;
        std::cout << "Sum polynomial:\n";
        sum.print(10);
        std::cout << "\n";

        // Test multiplication
        FMPoly product = poly * poly2;
        std::cout << "Product polynomial:\n";
        product.print(10);
        std::cout << "\n";

        // Export to JSON
        poly.exportToJson("fmpoly_example");
        std::cout << "Exported to fmpoly_example.json\n\n";

        // Show metadata (compatibility with MultivariablePolynomial)
        std::cout << "Polynomial metadata:\n";
        std::cout << "- Number of x variables: " << poly.getNumXVariables() << "\n";
        std::cout << "- Max degrees: [";
        auto maxDegrees = poly.getMaxXDegrees();
        for (size_t i = 0; i < maxDegrees.size(); i++) {
            std::cout << maxDegrees[i];
            if (i < maxDegrees.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
        std::cout << "- Is zero: " << (poly.isZero() ? "yes" : "no") << "\n\n";

        // Test clear functionality
        FMPoly testPoly(1, 2);
        testPoly.setCoefficient(1, {1}, 5);
        std::cout << "Before clear - is zero: " << (testPoly.isZero() ? "yes" : "no") << "\n";
        testPoly.clear();
        std::cout << "After clear - is zero: " << (testPoly.isZero() ? "yes" : "no") << "\n\n";

        std::cout << "✓ FMPoly example completed successfully!\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}