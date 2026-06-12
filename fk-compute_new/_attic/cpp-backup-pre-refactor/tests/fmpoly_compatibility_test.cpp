#include "fk/fmpoly.hpp"
#include "fk/multivariable_polynomial.hpp"
#include <iostream>
#include <cassert>

void testBasicCompatibility() {
    std::cout << "Testing basic interface compatibility...\n";

    // Create both polynomial types with same parameters
    FMPoly fmpoly(2, 3);
    MultivariablePolynomial mvpoly(2, 3);

    // Test setting coefficients
    fmpoly.setCoefficient(2, {1, 2}, 5);
    mvpoly.setCoefficient(2, {1, 2}, 5);

    fmpoly.setCoefficient(-1, {0, 1}, 3);
    mvpoly.setCoefficient(-1, {0, 1}, 3);

    // Test getting coefficients
    assert(fmpoly.getCoefficient(2, {1, 2}) == mvpoly.getCoefficient(2, {1, 2}));
    assert(fmpoly.getCoefficient(-1, {0, 1}) == mvpoly.getCoefficient(-1, {0, 1}));
    assert(fmpoly.getCoefficient(0, {0, 0}) == mvpoly.getCoefficient(0, {0, 0})); // Should be 0

    std::cout << "✓ Basic coefficient operations match\n";

    // Test adding to coefficients
    fmpoly.addToCoefficient(2, {1, 2}, 2);
    mvpoly.addToCoefficient(2, {1, 2}, 2);

    assert(fmpoly.getCoefficient(2, {1, 2}) == mvpoly.getCoefficient(2, {1, 2}));
    assert(fmpoly.getCoefficient(2, {1, 2}) == 7); // 5 + 2 = 7

    std::cout << "✓ addToCoefficient operations match\n";

    // Test metadata
    assert(fmpoly.getNumXVariables() == mvpoly.getNumXVariables());
    assert(fmpoly.getMaxXDegrees() == mvpoly.getMaxXDegrees());
    assert(fmpoly.getBlockSizes() == mvpoly.getBlockSizes());

    std::cout << "✓ Metadata compatibility confirmed\n";
}

void testArithmeticOperations() {
    std::cout << "\nTesting arithmetic operations...\n";

    // Create test polynomials
    FMPoly fm1(2, 3), fm2(2, 3);
    MultivariablePolynomial mv1(2, 3), mv2(2, 3);

    // Set same coefficients
    fm1.setCoefficient(2, {1, 1}, 3);
    mv1.setCoefficient(2, {1, 1}, 3);
    fm1.setCoefficient(1, {0, 2}, 2);
    mv1.setCoefficient(1, {0, 2}, 2);

    fm2.setCoefficient(1, {1, 0}, 4);
    mv2.setCoefficient(1, {1, 0}, 4);
    fm2.setCoefficient(3, {0, 0}, 1);
    mv2.setCoefficient(3, {0, 0}, 1);

    // Test addition
    FMPoly fm_sum = fm1 + fm2;
    MultivariablePolynomial mv_sum = mv1 + mv2;

    // Check that some key coefficients match
    assert(fm_sum.getCoefficient(2, {1, 1}) == mv_sum.getCoefficient(2, {1, 1}));
    assert(fm_sum.getCoefficient(1, {0, 2}) == mv_sum.getCoefficient(1, {0, 2}));
    assert(fm_sum.getCoefficient(1, {1, 0}) == mv_sum.getCoefficient(1, {1, 0}));
    assert(fm_sum.getCoefficient(3, {0, 0}) == mv_sum.getCoefficient(3, {0, 0}));

    std::cout << "✓ Addition operations produce compatible results\n";

    // Test subtraction
    FMPoly fm_diff = fm1 - fm2;
    MultivariablePolynomial mv_diff = mv1 - mv2;

    assert(fm_diff.getCoefficient(2, {1, 1}) == mv_diff.getCoefficient(2, {1, 1}));
    assert(fm_diff.getCoefficient(1, {1, 0}) == mv_diff.getCoefficient(1, {1, 0}));

    std::cout << "✓ Subtraction operations produce compatible results\n";
}

void testZeroOperations() {
    std::cout << "\nTesting zero operations...\n";

    FMPoly fmpoly(1, 2);
    MultivariablePolynomial mvpoly(1, 2);

    // Both should start as zero
    assert(fmpoly.isZero() == mvpoly.isZero());
    std::cout << "✓ Initial zero state matches\n";

    // Add some coefficients
    fmpoly.setCoefficient(1, {1}, 5);
    mvpoly.setCoefficient(1, {1}, 5);

    assert(fmpoly.isZero() == mvpoly.isZero());
    assert(!fmpoly.isZero()); // Should not be zero

    // Clear both
    fmpoly.clear();
    mvpoly.clear();

    assert(fmpoly.isZero() == mvpoly.isZero());
    assert(fmpoly.isZero()); // Should be zero again

    std::cout << "✓ Clear and zero detection match\n";
}

int main() {
    std::cout << "=== FMPoly Compatibility Test ===\n\n";

    try {
        testBasicCompatibility();
        testArithmeticOperations();
        testZeroOperations();

        std::cout << "\n🎉 All compatibility tests passed!\n";
        std::cout << "FMPoly can be used as a drop-in replacement for MultivariablePolynomial\n";

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "❌ Test failed with unknown exception" << std::endl;
        return 1;
    }
}