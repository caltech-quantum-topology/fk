#include <iostream>
#include "fk/multivariable_polynomial.hpp"

int main() {
    // Create a polynomial in 3 variables with default settings
    MultivariablePolynomial poly(3);

    // Set some coefficients
    // poly = 2*q^1*x₁¹*x₂⁰*x₃⁰ + 3*q^(-1)*x₁⁰*x₂¹*x₃⁰
    poly.setoCoefficient(1,{1,0,0},2);
    poly.getCoefficient(-1,{0, 1, 0},3);
  /*
    poly.getCoefficient({1, 0, 0})[1] = 2;   // coefficient of q¹
    poly.getCoefficient({0, 1, 0})[-1] = 3;  // coefficient of q⁻¹

    // Create another polynomial
    MultivariablePolynomial poly2(3);
    poly2.getCoefficient({0, 0, 1})[0] = 1;  // 1*q⁰*x₁⁰*x₂⁰*x₃¹

    // Perform arithmetic
    auto sum = poly + poly2;
    auto product = poly * poly2;

    // Check results
    std::cout << "Sum has " << sum.getTermCount() << " terms" << std::endl;
    std::cout << "Product has " << product.getTermCount() << " terms" << std::endl;
  */
    return 0;
}
