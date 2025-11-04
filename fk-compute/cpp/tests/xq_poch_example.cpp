#include "fk/multivariable_polynomial.hpp"
#include "fk/qalg_links.hpp"
#include <iostream>

void print_pterms(std::vector<bilvector<int>> polynomial_terms) {
  for ( auto i = 0; i < polynomial_terms.size(); ++i) {
      std::cout<<"x^"<<i<<" :";  
       polynomial_terms[i].print();
   }
}


int main() {

  std::vector<bilvector<int>> polynomial_terms(5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;

  print_pterms(polynomial_terms);

  MultivariablePolynomial poly(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  poly.print();
  print_pterms(polynomial_terms);

  auto xqpoch = qpochhammer_xq_q(5,1);
  xqpoch.print(100);
  MultivariablePolynomial test(xqpoch, 3,0);
  test.print(100);
  test = MultivariablePolynomial(xqpoch, 3,1);
  test.print(100);
  test = MultivariablePolynomial(xqpoch, 3,2);
  test.print(100);
  auto config = test.getCoefficientMap();
  auto newPoly = test.invertVariable(2);
  newPoly.print(100);

  MultivariablePolynomial test_inverse_qpochhammer = inverse_qpochhammer_xq_q(7,2,6);
  std::cout<<"======"<<std::endl;
  test_inverse_qpochhammer.print(200);

  return 0;
}
