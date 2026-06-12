#include "fk/fmpoly.hpp"
#include <iostream>
#include <vector>

void debugFMPoly() {
    std::cout << "=== Debugging FMPoly coefficient storage ===\n\n";

    try {
        // Test the exact sequence from the failing test
        std::cout << "=== Reproducing original test sequence ===\n";
        FMPoly poly(2, 5);

        std::cout << "Step 1: Setting 3*q^(-2)*x1^1*x2^2\n";
        poly.setCoefficient(-2, {1, 2}, 3);
        std::cout << "Check q^(-2): " << poly.getCoefficient(-2, {1, 2}) << "\n\n";

        std::cout << "Step 2: Setting 5*q^0*x1^1*x2^2\n";
        poly.setCoefficient(0, {1, 2}, 5);
        std::cout << "Check q^(-2): " << poly.getCoefficient(-2, {1, 2}) << "\n";
        std::cout << "Check q^0: " << poly.getCoefficient(0, {1, 2}) << "\n\n";

        std::cout << "Step 3: Setting 7*q^2*x1^1*x2^2\n";
        poly.setCoefficient(2, {1, 2}, 7);
        std::cout << "Check q^(-2): " << poly.getCoefficient(-2, {1, 2}) << "\n";
        std::cout << "Check q^0: " << poly.getCoefficient(0, {1, 2}) << "\n";
        std::cout << "Check q^2: " << poly.getCoefficient(2, {1, 2}) << "\n\n";

        std::cout << "Step 4: Setting 2*q^4*x1^1*x2^2\n";
        poly.setCoefficient(4, {1, 2}, 2);
        std::cout << "Check q^(-2): " << poly.getCoefficient(-2, {1, 2}) << "\n";
        std::cout << "Check q^0: " << poly.getCoefficient(0, {1, 2}) << "\n";
        std::cout << "Check q^2: " << poly.getCoefficient(2, {1, 2}) << "\n";
        std::cout << "Check q^4: " << poly.getCoefficient(4, {1, 2}) << "\n\n";

        std::cout << "getQPolynomial result:\n";
        std::vector<int> qPoly = poly.getQPolynomial({1, 2});
        std::cout << "Size: " << qPoly.size() << "\n";
        std::cout << "Coefficients: [";
        for (size_t i = 0; i < qPoly.size(); i++) {
            std::cout << qPoly[i];
            if (i < qPoly.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n\n";

        // Add different x-monomial
        std::cout << "Step 5: Setting 10*q^1*x1^0*x2^0 (different monomial)\n";
        poly.setCoefficient(1, {0, 0}, 10);
        std::cout << "Check q^(-2) for {1,2}: " << poly.getCoefficient(-2, {1, 2}) << "\n";
        std::cout << "Check q^1 for {0,0}: " << poly.getCoefficient(1, {0, 0}) << "\n\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    debugFMPoly();
    return 0;
}