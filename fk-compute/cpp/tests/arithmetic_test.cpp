#include "fk/multivariable_polynomial.hpp"
#include <iostream>

int main() {
  std::cout << "=== MultivariablePolynomial Arithmetic Test ===\n\n";

  // Create first polynomial: P1(q, x1, x2) = 3*q^2*x1^1*x2^0 + 2*q^0*x1^0*x2^1
  MultivariablePolynomial poly1(2, 3);
  poly1.setCoefficient(2, {1, 0}, 3); // 3*q^2*x1^1
  poly1.setCoefficient(0, {0, 1}, 2); // 2*q^0*x2^1

  std::cout << "Polynomial P1:\n";
  poly1.print();
  std::cout << "\n";

  // Create second polynomial: P2(q, x1, x2) = 1*q^1*x1^1*x2^0 + 4*q^0*x1^0*x2^1
  MultivariablePolynomial poly2(2, 3);
  poly2.setCoefficient(1, {1, 0}, 1); // 1*q^1*x1^1
  poly2.setCoefficient(0, {0, 1}, 4); // 4*q^0*x2^1

  std::cout << "Polynomial P2:\n";
  poly2.print();
  std::cout << "\n";

  // Test addition: P1 + P2
  std::cout << "=== Addition Test: P1 + P2 ===\n";
  MultivariablePolynomial sum = poly1 + poly2;
  std::cout << "Result:\n";
  sum.print();
  std::cout << "Expected: 3*q^2*x1^1 + 1*q^1*x1^1 + 6*q^0*x2^1\n\n";

  // Test subtraction: P1 - P2
  std::cout << "=== Subtraction Test: P1 - P2 ===\n";
  MultivariablePolynomial diff = poly1 - poly2;
  std::cout << "Result:\n";
  diff.print();
  std::cout << "Expected: 3*q^2*x1^1 + (-1)*q^1*x1^1 + (-2)*q^0*x2^1\n\n";

  // Test multiplication: P1 * P2
  std::cout << "=== Multiplication Test: P1 * P2 ===\n";
  MultivariablePolynomial product = poly1 * poly2;
  std::cout << "Result:\n";
  product.print();
  std::cout << "Expected terms:\n";
  std::cout << "  3*q^2*x1^1 * 1*q^1*x1^1 = 3*q^3*x1^2\n";
  std::cout << "  3*q^2*x1^1 * 4*q^0*x2^1 = 12*q^2*x1^1*x2^1\n";
  std::cout << "  2*q^0*x2^1 * 1*q^1*x1^1 = 2*q^1*x1^1*x2^1\n";
  std::cout << "  2*q^0*x2^1 * 4*q^0*x2^1 = 8*q^0*x2^2\n\n";

  // Test in-place operations
  std::cout << "=== In-place Operations Test ===\n";

  // Test +=
  MultivariablePolynomial poly3 = poly1; // Copy P1
  std::cout << "Before P1 += P2:\n";
  poly3.print();
  poly3 += poly2;
  std::cout << "After P1 += P2:\n";
  poly3.print();
  std::cout << "\n";

  // Test -=
  MultivariablePolynomial poly4 = poly1; // Copy P1
  std::cout << "Before P1 -= P2:\n";
  poly4.print();
  poly4 -= poly2;
  std::cout << "After P1 -= P2:\n";
  poly4.print();
  std::cout << "\n";

  // Test *=
  MultivariablePolynomial poly5 = poly1; // Copy P1
  std::cout << "Before P1 *= P2:\n";
  poly5.print();
  poly5 *= poly2;
  std::cout << "After P1 *= P2:\n";
  poly5.print();
  std::cout << "\n";

  // Test with negative coefficients
  std::cout << "=== Negative Coefficient Test ===\n";
  MultivariablePolynomial poly6(2, 2);
  poly6.setCoefficient(-1, {1, 0}, -3); // -3*q^(-1)*x1^1
  poly6.setCoefficient(1, {0, 1}, 5);   // 5*q^1*x2^1

  std::cout << "Polynomial with negative coefficients:\n";
  poly6.print();
  std::cout << "\n";

  MultivariablePolynomial poly7(2, 2);
  poly7.setCoefficient(2, {1, 0}, 2);   // 2*q^2*x1^1
  poly7.setCoefficient(-1, {0, 1}, -1); // -1*q^(-1)*x2^1

  std::cout << "Second polynomial:\n";
  poly7.print();
  std::cout << "\n";

  std::cout << "Sum of polynomials with negative coefficients:\n";
  MultivariablePolynomial negSum = poly6 + poly7;
  negSum.print();
  std::cout << "\n";

  // Test error handling
  std::cout << "=== Error Handling Test ===\n";
  try {
    MultivariablePolynomial incompatible(3, 2); // Different number of variables
    poly1 + incompatible;                       // This should throw
    std::cout << "ERROR: Should have thrown an exception!\n";
  } catch (const std::invalid_argument &e) {
    std::cout << "Correctly caught exception: " << e.what() << "\n";
  }

  std::cout << "\n=== Test Complete ===\n";
  return 0;
}