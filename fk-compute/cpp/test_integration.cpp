#include "fk/fmpoly.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Testing FMPoly and QPolynomial integration..." << std::endl;

    // Create a 2-variable FMPoly
    FMPoly poly(2, 5);

    // Set some coefficients manually
    poly.setCoefficient(1, {1, 0}, 3);   // 3q^1 * x1^1 * x2^0
    poly.setCoefficient(-1, {1, 0}, 2);  // 2q^(-1) * x1^1 * x2^0
    poly.setCoefficient(2, {1, 0}, -1);  // -q^2 * x1^1 * x2^0

    std::cout << "Initial FMPoly:" << std::endl;
    poly.print();

    // Test 1: Get QPolynomial for x1^1, x2^0
    std::cout << "\nTest 1 - Get QPolynomial for x1^1*x2^0:" << std::endl;
    QPolynomial qpoly1 = poly.getQPolynomialObject({1, 0});
    std::cout << "QPolynomial: ";
    qpoly1.print();

    // Test 2: Create a new QPolynomial and set it
    std::cout << "\nTest 2 - Set QPolynomial for x1^0*x2^1:" << std::endl;
    QPolynomial qpoly2({1, 0, 2}, -1);  // q^(-1) + 2q^1
    std::cout << "Setting QPolynomial: ";
    qpoly2.print();
    poly.setQPolynomial({0, 1}, qpoly2);

    std::cout << "FMPoly after setting:" << std::endl;
    poly.print();

    // Test 3: Add QPolynomial to existing
    std::cout << "\nTest 3 - Add QPolynomial to existing x1^1*x2^0:" << std::endl;
    QPolynomial qpoly3({5}, 0);  // 5q^0
    std::cout << "Adding QPolynomial: ";
    qpoly3.print();
    poly.addQPolynomial({1, 0}, qpoly3);

    std::cout << "FMPoly after adding:" << std::endl;
    poly.print();

    // Test 4: Multiply QPolynomial with existing
    std::cout << "\nTest 4 - Multiply QPolynomial with existing x1^0*x2^1:" << std::endl;
    QPolynomial qpoly4({1, 1}, 0);  // 1 + q
    std::cout << "Multiplying by QPolynomial: ";
    qpoly4.print();
    poly.multiplyQPolynomial({0, 1}, qpoly4);

    std::cout << "FMPoly after multiplying:" << std::endl;
    poly.print();

    // Test 5: Verify the result by getting QPolynomial back
    std::cout << "\nTest 5 - Verify result for x1^0*x2^1:" << std::endl;
    QPolynomial result = poly.getQPolynomialObject({0, 1});
    std::cout << "Final QPolynomial for x1^0*x2^1: ";
    result.print();

    // Test 6: Verify original polynomial comparison
    std::cout << "\nTest 6 - Compare with original vector method:" << std::endl;
    std::vector<int> vecResult = poly.getQPolynomial({1, 0});
    QPolynomial objResult = poly.getQPolynomialObject({1, 0});

    std::cout << "Vector method coefficients: ";
    for (size_t i = 0; i < vecResult.size(); i++) {
        std::cout << vecResult[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "QPolynomial object: ";
    objResult.print();

    std::cout << "\nAll integration tests completed!" << std::endl;
    return 0;
}