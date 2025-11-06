#include <iostream>
#include "fk/multivariable_polynomial.hpp"

void print(std::string s){
    std::cout << s << std::endl;
}

int main() {
    // Create a polynomial in 3 variables with default settings
    MultivariablePolynomial poly(3);


    // Set some coefficients
    // poly = 2*q^1*x₁¹*x₂⁰*x₃⁰ + 3*q^(-1)*x₁⁰*x₂¹*x₃⁰
    poly.setCoefficient(1,{1, 0, 0},2);   // coefficient of q¹
    poly.setCoefficient(-1,{0, 1, 0},3);
    print("Poly1");
    poly.print();

    // Create another polynomial
    MultivariablePolynomial poly2(3);
    poly2.setCoefficient(0,{0, 0, 1},1);  // 1*q⁰*x₁⁰*x₂⁰*x₃¹
    print("Poly2");
    poly2.print();

    // Perform arithmetic
    auto sum = poly + poly2;
    auto product = poly * poly2;

    print("Sum");
    sum.print();

    print("Product");
    product.print();

    // Test the evaluate function
    print("\n=== Testing evaluate function ===");

    // Create a simple polynomial: 2*q^1*x₁ + 3*q^(-1)*x₂ + 5*q^0
    MultivariablePolynomial testPoly(2);
    testPoly.setCoefficient(1, {1, 0}, 2);   // 2*q^1*x₁
    testPoly.setCoefficient(-1, {0, 1}, 3);  // 3*q^(-1)*x₂
    testPoly.setCoefficient(0, {0, 0}, 5);   // 5*q^0

    print("Test polynomial:");
    testPoly.print();

    // Evaluate at point (2, 1) -> 2*q^1*2 + 3*q^(-1)*1 + 5*q^0 = 4*q^1 + 3*q^(-1) + 5*q^0
    std::vector<int> point = {2, 1};
    print("Evaluating at point (2, 1):");

    auto result = testPoly.evaluate(point);

    std::cout << "Result: ";
    for (int qPower = result.getMaxNegativeIndex(); qPower <= result.getMaxPositiveIndex(); ++qPower) {
        int coeff = result[qPower];
        if (coeff != 0) {
            if (qPower == 0) {
                std::cout << coeff;
            } else if (qPower == 1) {
                std::cout << coeff << "*q";
            } else if (qPower == -1) {
                std::cout << coeff << "*q^(-1)";
            } else {
                std::cout << coeff << "*q^(" << qPower << ")";
            }
            std::cout << " + ";
        }
    }
    std::cout << std::endl;

    // Test error handling - division by zero
    print("\nTesting division by zero error:");
    MultivariablePolynomial errorPoly(2);
    errorPoly.setCoefficient(0, {-1, 0}, 1);  // 1*q^0*x₁^(-1)

    try {
        std::vector<int> zeroPoint = {0, 1};
        auto errorResult = errorPoly.evaluate(zeroPoint);
    } catch (const std::domain_error& e) {
        std::cout << "Caught expected error: " << e.what() << std::endl;
    }

    return 0;
}
