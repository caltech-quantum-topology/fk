#include "fk/fmpoly.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Testing QPolynomial class..." << std::endl;

    // Test 1: Default constructor
    QPolynomial p1;
    std::cout << "Test 1 - Default constructor: ";
    p1.print();

    // Test 2: Constructor with coefficients
    std::vector<int> coeffs = {1, -2, 3};
    QPolynomial p2(coeffs, 0);
    std::cout << "Test 2 - Constructor with [1, -2, 3], minPower=0: ";
    p2.print();

    // Test 3: Constructor with negative minPower
    QPolynomial p3(coeffs, -2);
    std::cout << "Test 3 - Constructor with [1, -2, 3], minPower=-2: ";
    p3.print();

    // Test 4: Set and get coefficients
    QPolynomial p4;
    p4.setCoefficient(2, 5);
    p4.setCoefficient(-1, 3);
    p4.setCoefficient(0, -2);
    std::cout << "Test 4 - Set coefficients q^2=5, q^(-1)=3, q^0=-2: ";
    p4.print();

    // Test 5: Addition
    QPolynomial p5 = p2 + p3;
    std::cout << "Test 5 - Addition p2 + p3: ";
    p5.print();

    // Test 5b: Addition with different minPowers
    QPolynomial p5b({1}, -3);    // q^(-3)
    QPolynomial p5c({2}, 2);     // 2q^2
    QPolynomial p5d = p5b + p5c; // Should be q^(-3) + 2q^2
    std::cout << "Test 5b - Addition q^(-3) + 2q^2: ";
    p5d.print();

    // Test 6: Multiplication
    QPolynomial p6({1, 1}, 0);  // 1 + q
    QPolynomial p7({1, -1}, 0); // 1 - q
    QPolynomial p8 = p6 * p7;   // Should be 1 - q^2
    std::cout << "Test 6 - Multiplication (1+q)*(1-q): ";
    p8.print();

    // Test 7: Evaluation
    int result = p2.evaluate(2);
    std::cout << "Test 7 - Evaluate p2 at q=2: " << result << std::endl;

    // Test 8: Get coefficient
    std::cout << "Test 8 - Get coefficient of q^1 in p2: " << p2.getCoefficient(1) << std::endl;

    std::cout << "All tests completed!" << std::endl;
    return 0;
}