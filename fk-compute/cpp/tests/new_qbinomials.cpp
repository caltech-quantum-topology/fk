#include "fk/multivariable_polynomial.hpp"
#include "fk/qalg_links.hpp"
#include <iostream>

int main() {

  std::cout<<"Case one: n > k > 0"<<std::endl;
  std::vector<bilvector<int>> polynomial_terms(5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }
  
  computePositiveQBinomial(polynomial_terms, 5, 3, false);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }

  polynomial_terms = std::vector<bilvector<int>> (5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;

  computePositiveQBinomial(polynomial_terms, 5, 3, true);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }

  MultivariablePolynomial poly(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  poly.print();
  
  poly.getQPolynomial({0}).print();
  poly.getQPolynomial({0}) *= QBinomial(5,3);
  poly.getQPolynomial({0}).print();
  poly.print();

  poly = MultivariablePolynomial(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  poly.print();
  
  poly.getQPolynomial({0}).print();
  poly.getQPolynomial({0}) *= QBinomial(5,3).invertExponents();
  poly.getQPolynomial({0}).print();
  poly.print();

  std::cout<<"Case two: 0 > n > k"<<std::endl;

  int n(-3),k(-5);

  polynomial_terms = std::vector<bilvector<int>> (5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;


  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }
  
  computeNegativeQBinomial(polynomial_terms, n, k, false);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }

  polynomial_terms = std::vector<bilvector<int>> (5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;

  computeNegativeQBinomial(polynomial_terms, n, k, true);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }


  poly =MultivariablePolynomial(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  
  poly.getQPolynomial({0}) *= QBinomial(n,k);
  poly.print(50);

  poly = MultivariablePolynomial(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  
  poly.getQPolynomial({0}) *= QBinomial(n,k).invertExponents();
  poly.print();



  std::cout<<"Case three: k > 0 > n"<<std::endl;

  polynomial_terms = std::vector<bilvector<int>> (5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;


  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }
  
  computeNegativeQBinomial(polynomial_terms, -5, 3, false);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }

  polynomial_terms = std::vector<bilvector<int>> (5, bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = 1;
  polynomial_terms[1][2] = 1;
  polynomial_terms[2][-1] = 2;

  computeNegativeQBinomial(polynomial_terms, -5, 3, true);

  for (std::size_t i = 0; i < polynomial_terms.size(); ++i) {
      std::cout << "x1^" << i << ": ";
      polynomial_terms[i].print();
      std::cout << "\n";
  }


  poly =MultivariablePolynomial(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  
  poly.getQPolynomial({0}) *= QBinomial(-5,3);
  poly.print(50);

  poly = MultivariablePolynomial(1,0);
  poly.setCoefficient(0,{0},1);
  poly.setCoefficient(2,{1},1);
  poly.setCoefficient(-1,{2},2);
  
  poly.getQPolynomial({0}) *= QBinomial(-5,3).invertExponents();
  poly.print(50);


  return 0;

}

