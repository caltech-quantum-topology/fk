#include "fk/multivariable_polynomial.hpp"
#include "fk/qalg_links.hpp"
#include <iostream>

int main() {
  // Create a polynomial in q, x1, x2 with max degree 3 for each x variable
  MultivariablePolynomial poly(2, 3);
  poly.setCoefficient(2, {1, 2}, 5);  // 5*q^2*x1^1*x2^2
  poly.setCoefficient(-1, {0, 1}, 3); // 3*q^(-1)*x1^0*x2^1
  poly.setCoefficient(3, {2, 0}, 7);  // 7*q^3*x1^2*x2^0
  poly.print(); 

  auto qb = QBinomialPositive(5,3);
  auto poly2 = poly * qb;
  poly2.print(50);
  auto poly3 = poly + qb;
  poly3.print(50);
  auto poly4 = poly - qb;
  poly4.print(50);

  qb.print();
  auto qb2 = qb.invertExponents();
  qb2.print();


  return 0;
}
