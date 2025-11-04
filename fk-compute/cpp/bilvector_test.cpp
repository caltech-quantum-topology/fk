#include <iostream>
#include <list>
#include <vector>
#include "fk/bilvector.hpp"
#include "fk/multivariable_polynomial.hpp"

// assumes bilvector<T>, makeLaurentPolynomial<T>, and
// operator+, operator*, operator+=, operator*= are defined in fk/bilvector.hpp

// Optional helper to pretty-print a Laurent polynomial
template <typename T>
void printPolynomial(const bilvector<T> &P,
                     int minExp,
                     int maxExp,
                     const std::string &name) {
  std::cout << name << "(q) = ";
  bool first = true;
  for (int e = minExp; e <= maxExp; ++e) {
    T c = P[e];
    if (c == 0) continue;
    if (!first) {
      std::cout << " + ";
    }
    first = false;
    std::cout << c << "*q^" << e;
  }
  if (first) {
    std::cout << "0";
  }
  std::cout << "\n";
}

int main() {
  // Example Laurent polynomial:
  //   P(q) = 3 q^{-2} - 5 q^{-1} + 2 + 7 q^3
  //
  // So we need exponents from -2 to 3.
  int minExp = -2;
  int maxExp = 3;

  bilvector<int> P = makeLaurentPolynomial<int>(minExp, maxExp, 0);

  P[-2] = 3;   // 3 q^{-2}
  P[-1] = -5;  // -5 q^{-1}
  P[0]  = 2;   // 2
  P[3]  = 7;   // 7 q^3

  std::cout << "Testing makeLaurentPolynomial\n";
  std::cout << "Constructed Laurent polynomial P(q) on range ["
            << minExp << ", " << maxExp << "]\n\n";

  // Print in a human-readable polynomial form
  printPolynomial(P, minExp, maxExp, "P");
  std::cout << "\n";

  // Print all coefficients including zeros
  std::cout << "Coefficients of P:\n";
  for (int e = minExp; e <= maxExp; ++e) {
    std::cout << "  coeff_P(q^" << e << ") = " << P[e] << "\n";
  }

  // ------------------------------------------------------------------------
  // Build a second Laurent polynomial Q and test addition/multiplication
  // ------------------------------------------------------------------------
  //
  // Let
  //   Q(q) = 1*q^{-1} + 1 + 1*q
  //
  // Exponent range: -1 to 1.
  bilvector<int> Q = makeLaurentPolynomial<int>(-1, 1, 0);
  Q[-1] = 1;
  Q[0]  = 1;
  Q[1]  = 1;

  std::cout << "\nConstructed Laurent polynomial Q(q) on range [-1, 1]\n";
  printPolynomial(Q, -1, 1, "Q");
  std::cout << "\n";

  // Test addition: R = P + Q
  bilvector<int> R = P + Q;
  int Rmin = R.getMaxNegativeIndex();
  int Rmax = R.getMaxPositiveIndex();

  std::cout << "Testing addition: R(q) = P(q) + Q(q)\n";
  printPolynomial(R, Rmin, Rmax, "R");
  std::cout << "Coefficients of R:\n";
  for (int e = Rmin; e <= Rmax; ++e) {
    std::cout << "  coeff_R(q^" << e << ") = " << R[e] << "\n";
  }
  std::cout << "\n";

  // Test multiplication: S = P * Q
  bilvector<int> S = P * Q;
  int Smin = S.getMaxNegativeIndex();
  int Smax = S.getMaxPositiveIndex();

  std::cout << "Testing multiplication: S(q) = P(q) * Q(q)\n";
  printPolynomial(S, Smin, Smax, "S");
  std::cout << "Coefficients of S:\n";
  for (int e = Smin; e <= Smax; ++e) {
    std::cout << "  coeff_S(q^" << e << ") = " << S[e] << "\n";
  }
  std::cout << "\n";

  // Also quickly test += and *=
  std::cout << "Testing operator+= and operator*=\n";
  bilvector<int> A = P;
  A += Q;
  int Amin = A.getMaxNegativeIndex();
  int Amax = A.getMaxPositiveIndex();
  printPolynomial(A, Amin, Amax, "A = P += Q");

  bilvector<int> B = P;
  B *= Q;
  int Bmin = B.getMaxNegativeIndex();
  int Bmax = B.getMaxPositiveIndex();
  printPolynomial(B, Bmin, Bmax, "B = P *= Q");
  std::cout << "\n";

  // ------------------------------------------------------------------------
  // Existing growth tests
  // ------------------------------------------------------------------------

  // Test automatic growth on the positive side
  std::cout << "\nTesting growth: setting coeff of q^5 = 42\n";
  P[5] = 42;
  std::cout << "  coeff_P(q^5) = " << P[5] << "\n";
  std::cout << "  maxPositiveIndex (recorded) = "
            << P.getMaxPositiveIndex() << "\n";

  // Test automatic growth on the negative side
  std::cout << "\nTesting growth: setting coeff of q^-4 = 11\n";
  P[-4] = 11;
  std::cout << "  coeff_P(q^-4) = " << P[-4] << "\n";
  std::cout << "  maxNegativeIndex (recorded) = "
            << P.getMaxNegativeIndex() << "\n";

  // Show a slightly extended range after growth
  std::cout << "\nExtended coefficients of P after growth (from -4 to 5):\n";
  for (int e = -4; e <= 5; ++e) {
    std::cout << "  coeff_P(q^" << e << ") = " << P[e] << "\n";
  }

  return 0;
}
