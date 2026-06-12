#include <iostream>
#include "fk/multivariable_polynomial.hpp"
#include "fk/inequality_solver.hpp"

void printSolutions(const std::set<IntegerPoint>& solutions) {
    std::cout << "Found " << solutions.size() << " integer solutions:" << std::endl;
    for (const auto& point : solutions) {
        std::cout << "(";
        for (size_t i = 0; i < point.coordinates.size(); ++i) {
            std::cout << point.coordinates[i];
            if (i < point.coordinates.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "=== Linear Inequality Solver Example ===" << std::endl << std::endl;

    // Example 1: Simple 2D system
    // x + y >= 1
    // x - y >= -2
    // x >= 0
    // y >= 0
    std::cout << "Example 1: 2D system" << std::endl;
    std::cout << "Inequalities:" << std::endl;
    std::cout << "  x + y >= 1" << std::endl;
    std::cout << "  x - y >= -2" << std::endl;
    std::cout << "  x >= 0" << std::endl;
    std::cout << "  y >= 0" << std::endl;
    std::cout << "Search bounds: x ∈ [0, 5], y ∈ [0, 5]" << std::endl << std::endl;

    std::vector<MultivariablePolynomial> inequalities1;

    // x + y - 1 >= 0 (so x + y >= 1)
    MultivariablePolynomial ineq1(2);
    ineq1.setCoefficient(0, {1, 0}, 1);  // coefficient of x
    ineq1.setCoefficient(0, {0, 1}, 1);  // coefficient of y
    ineq1.setCoefficient(0, {0, 0}, -1); // constant term
    inequalities1.push_back(ineq1);

    // x - y + 2 >= 0 (so x - y >= -2)
    MultivariablePolynomial ineq2(2);
    ineq2.setCoefficient(0, {1, 0}, 1);  // coefficient of x
    ineq2.setCoefficient(0, {0, 1}, -1); // coefficient of y
    ineq2.setCoefficient(0, {0, 0}, 2);  // constant term
    inequalities1.push_back(ineq2);

    // x >= 0
    MultivariablePolynomial ineq3(2);
    ineq3.setCoefficient(0, {1, 0}, 1);  // coefficient of x
    inequalities1.push_back(ineq3);

    // y >= 0
    MultivariablePolynomial ineq4(2);
    ineq4.setCoefficient(0, {0, 1}, 1);  // coefficient of y
    inequalities1.push_back(ineq4);

    std::vector<std::pair<int, int>> bounds1 = {{0, 5}, {0, 5}};
    auto solutions1 = findIntegerSolutions(inequalities1, bounds1);
    printSolutions(solutions1);

    // Example 2: 3D system
    // x + y + z >= 2
    // x - y >= 0
    // z >= 1
    std::cout << "Example 2: 3D system" << std::endl;
    std::cout << "Inequalities:" << std::endl;
    std::cout << "  x + y + z >= 2" << std::endl;
    std::cout << "  x - y >= 0" << std::endl;
    std::cout << "  z >= 1" << std::endl;
    std::cout << "Search bounds: x ∈ [0, 3], y ∈ [0, 3], z ∈ [0, 3]" << std::endl << std::endl;

    std::vector<MultivariablePolynomial> inequalities2;

    // x + y + z - 2 >= 0
    MultivariablePolynomial ineq5(3);
    ineq5.setCoefficient(0, {1, 0, 0}, 1);  // coefficient of x
    ineq5.setCoefficient(0, {0, 1, 0}, 1);  // coefficient of y
    ineq5.setCoefficient(0, {0, 0, 1}, 1);  // coefficient of z
    ineq5.setCoefficient(0, {0, 0, 0}, -2); // constant term
    inequalities2.push_back(ineq5);

    // x - y >= 0
    MultivariablePolynomial ineq6(3);
    ineq6.setCoefficient(0, {1, 0, 0}, 1);  // coefficient of x
    ineq6.setCoefficient(0, {0, 1, 0}, -1); // coefficient of y
    inequalities2.push_back(ineq6);

    // z - 1 >= 0
    MultivariablePolynomial ineq7(3);
    ineq7.setCoefficient(0, {0, 0, 1}, 1);  // coefficient of z
    ineq7.setCoefficient(0, {0, 0, 0}, -1); // constant term
    inequalities2.push_back(ineq7);

    std::vector<std::pair<int, int>> bounds2 = {{0, 3}, {0, 3}, {0, 3}};
    auto solutions2 = findIntegerSolutions(inequalities2, bounds2);
    printSolutions(solutions2);

    // Example 3: Triangle constraint
    // x + y <= 4 (equivalent to -x - y + 4 >= 0)
    // x >= 1
    // y >= 1
    std::cout << "Example 3: Triangle region" << std::endl;
    std::cout << "Inequalities:" << std::endl;
    std::cout << "  x + y <= 4" << std::endl;
    std::cout << "  x >= 1" << std::endl;
    std::cout << "  y >= 1" << std::endl;
    std::cout << "Search bounds: x ∈ [0, 5], y ∈ [0, 5]" << std::endl << std::endl;

    std::vector<MultivariablePolynomial> inequalities3;

    // -x - y + 4 >= 0 (equivalent to x + y <= 4)
    MultivariablePolynomial ineq8(2);
    ineq8.setCoefficient(0, {1, 0}, -1); // coefficient of x
    ineq8.setCoefficient(0, {0, 1}, -1); // coefficient of y
    ineq8.setCoefficient(0, {0, 0}, 4);  // constant term
    inequalities3.push_back(ineq8);

    // x - 1 >= 0
    MultivariablePolynomial ineq9(2);
    ineq9.setCoefficient(0, {1, 0}, 1);  // coefficient of x
    ineq9.setCoefficient(0, {0, 0}, -1); // constant term
    inequalities3.push_back(ineq9);

    // y - 1 >= 0
    MultivariablePolynomial ineq10(2);
    ineq10.setCoefficient(0, {0, 1}, 1);  // coefficient of y
    ineq10.setCoefficient(0, {0, 0}, -1); // constant term
    inequalities3.push_back(ineq10);

    std::vector<std::pair<int, int>> bounds3 = {{0, 5}, {0, 5}};
    auto solutions3 = findIntegerSolutions(inequalities3, bounds3);
    printSolutions(solutions3);

    return 0;
}
