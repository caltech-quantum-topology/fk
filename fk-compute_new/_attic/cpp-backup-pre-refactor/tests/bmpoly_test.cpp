#include "fk/bmpoly.hpp"
#include <iostream>
#include <cassert>
#include <stdexcept>

// Test helper macros
#define TEST_ASSERT(condition, message) \
  do { \
    if (!(condition)) { \
      std::cerr << "ASSERTION FAILED: " << message << " at line " << __LINE__ << std::endl; \
      return false; \
    } \
  } while(0)

#define RUN_TEST(test_func) \
  do { \
    std::cout << "Running " << #test_func << "... "; \
    if (test_func()) { \
      std::cout << "PASSED" << std::endl; \
      passed++; \
    } else { \
      std::cout << "FAILED" << std::endl; \
      failed++; \
    } \
    total++; \
  } while(0)

// Test basic constructor and getters
bool test_basic_constructor() {
  BMPoly poly(2, 5);

  TEST_ASSERT(poly.getNumXVariables() == 2, "Number of variables should be 2");
  TEST_ASSERT(poly.getMaxXDegrees().size() == 2, "Max degrees vector size should be 2");
  TEST_ASSERT(poly.getMaxXDegrees()[0] == 5, "Max degree for variable 0 should be 5");
  TEST_ASSERT(poly.getMaxXDegrees()[1] == 5, "Max degree for variable 1 should be 5");
  TEST_ASSERT(poly.isZero(), "New polynomial should be zero");

  return true;
}

// Test constructor with custom max degrees
bool test_constructor_with_max_degrees() {
  std::vector<int> maxDegrees = {3, 7, 4};
  BMPoly poly(3, 0, maxDegrees);

  TEST_ASSERT(poly.getNumXVariables() == 3, "Number of variables should be 3");
  TEST_ASSERT(poly.getMaxXDegrees()[0] == 3, "Max degree for variable 0 should be 3");
  TEST_ASSERT(poly.getMaxXDegrees()[1] == 7, "Max degree for variable 1 should be 7");
  TEST_ASSERT(poly.getMaxXDegrees()[2] == 4, "Max degree for variable 2 should be 4");

  return true;
}

// Test coefficient operations
bool test_coefficient_operations() {
  BMPoly poly(2, 5);

  // Test setting coefficients
  poly.setCoefficient(2, {1, 0}, 3);  // 3*q^2*x1^1
  poly.setCoefficient(0, {0, 1}, 2);  // 2*q^0*x2^1
  poly.setCoefficient(-1, {1, 1}, 5); // 5*q^(-1)*x1^1*x2^1

  TEST_ASSERT(poly.getCoefficient(2, {1, 0}) == 3, "Coefficient should be 3");
  TEST_ASSERT(poly.getCoefficient(0, {0, 1}) == 2, "Coefficient should be 2");
  TEST_ASSERT(poly.getCoefficient(-1, {1, 1}) == 5, "Coefficient should be 5");
  TEST_ASSERT(poly.getCoefficient(1, {0, 0}) == 0, "Unset coefficient should be 0");

  // Test adding to coefficients
  poly.addToCoefficient(2, {1, 0}, 7);
  TEST_ASSERT(poly.getCoefficient(2, {1, 0}) == 10, "Coefficient should be 10 after adding 7");

  TEST_ASSERT(!poly.isZero(), "Polynomial with coefficients should not be zero");

  return true;
}

// Test negative exponent support
bool test_negative_exponents() {
  BMPoly poly(2, 3);

  // Set coefficients with negative x-exponents
  poly.setCoefficient(1, {-2, 1}, 4);   // 4*q^1*x1^(-2)*x2^1
  poly.setCoefficient(0, {1, -1}, 3);   // 3*q^0*x1^1*x2^(-1)
  poly.setCoefficient(-1, {-1, -2}, 2); // 2*q^(-1)*x1^(-1)*x2^(-2)

  TEST_ASSERT(poly.getCoefficient(1, {-2, 1}) == 4, "Negative x-exponent coefficient should be 4");
  TEST_ASSERT(poly.getCoefficient(0, {1, -1}) == 3, "Negative x-exponent coefficient should be 3");
  TEST_ASSERT(poly.getCoefficient(-1, {-1, -2}) == 2, "Negative x-exponent coefficient should be 2");

  return true;
}

// Test getQPolynomial access
bool test_q_polynomial_access() {
  BMPoly poly(2, 3);

  // Get reference and modify directly
  auto &qpoly = poly.getQPolynomial({1, 0});
  qpoly[2] = 5;
  qpoly[-1] = 3;

  TEST_ASSERT(poly.getCoefficient(2, {1, 0}) == 5, "Direct qpolynomial modification should work");
  TEST_ASSERT(poly.getCoefficient(-1, {1, 0}) == 3, "Direct qpolynomial modification should work");

  // Test const access
  const BMPoly &constPoly = poly;
  const auto &constQPoly = constPoly.getQPolynomial({1, 0});
  TEST_ASSERT(constQPoly[2] == 5, "Const access should return same values");

  return true;
}

// Test polynomial addition
bool test_addition() {
  BMPoly poly1(2, 5);
  BMPoly poly2(2, 5);

  // Set up first polynomial: 3*q^2*x1^1 + 2*q^0*x2^1
  poly1.setCoefficient(2, {1, 0}, 3);
  poly1.setCoefficient(0, {0, 1}, 2);

  // Set up second polynomial: 1*q^2*x1^1 + 4*q^0*x2^1 + 5*q^1*x1^1*x2^1
  poly2.setCoefficient(2, {1, 0}, 1);
  poly2.setCoefficient(0, {0, 1}, 4);
  poly2.setCoefficient(1, {1, 1}, 5);

  BMPoly sum = poly1 + poly2;

  TEST_ASSERT(sum.getCoefficient(2, {1, 0}) == 4, "Sum coefficient should be 4");
  TEST_ASSERT(sum.getCoefficient(0, {0, 1}) == 6, "Sum coefficient should be 6");
  TEST_ASSERT(sum.getCoefficient(1, {1, 1}) == 5, "Sum coefficient should be 5");

  // Test in-place addition
  poly1 += poly2;
  TEST_ASSERT(poly1.getCoefficient(2, {1, 0}) == 4, "In-place sum coefficient should be 4");
  TEST_ASSERT(poly1.getCoefficient(0, {0, 1}) == 6, "In-place sum coefficient should be 6");
  TEST_ASSERT(poly1.getCoefficient(1, {1, 1}) == 5, "In-place sum coefficient should be 5");

  return true;
}

// Test polynomial subtraction
bool test_subtraction() {
  BMPoly poly1(2, 5);
  BMPoly poly2(2, 5);

  // Set up polynomials
  poly1.setCoefficient(2, {1, 0}, 3);
  poly1.setCoefficient(0, {0, 1}, 2);

  poly2.setCoefficient(2, {1, 0}, 1);
  poly2.setCoefficient(0, {0, 1}, 4);
  poly2.setCoefficient(1, {1, 1}, 5);

  BMPoly diff = poly1 - poly2;

  TEST_ASSERT(diff.getCoefficient(2, {1, 0}) == 2, "Difference coefficient should be 2");
  TEST_ASSERT(diff.getCoefficient(0, {0, 1}) == -2, "Difference coefficient should be -2");
  TEST_ASSERT(diff.getCoefficient(1, {1, 1}) == -5, "Difference coefficient should be -5");

  return true;
}

// Test polynomial multiplication
bool test_multiplication() {
  BMPoly poly1(2, 5);
  BMPoly poly2(2, 5);

  // Set up first polynomial: 2*q^1*x1^1
  poly1.setCoefficient(1, {1, 0}, 2);

  // Set up second polynomial: 3*q^2*x2^1
  poly2.setCoefficient(2, {0, 1}, 3);

  BMPoly product = poly1 * poly2;

  // Result should be: 6*q^3*x1^1*x2^1
  TEST_ASSERT(product.getCoefficient(3, {1, 1}) == 6, "Product coefficient should be 6");
  TEST_ASSERT(product.getCoefficient(1, {1, 0}) == 0, "Other coefficients should be 0");
  TEST_ASSERT(product.getCoefficient(2, {0, 1}) == 0, "Other coefficients should be 0");

  return true;
}

// Test operations with bilvectors
bool test_bilvector_operations() {
  BMPoly poly(2, 5);
  bilvector<int> qpoly(0, 1, 20, 0);

  // Set up polynomial: 2*q^1*x1^1
  poly.setCoefficient(1, {1, 0}, 2);

  // Set up q-polynomial: 3*q^0 + 4*q^2
  qpoly[0] = 3;
  qpoly[2] = 4;

  // Test addition with bilvector (adds to x^0 term)
  BMPoly sum = poly + qpoly;
  TEST_ASSERT(sum.getCoefficient(0, {0, 0}) == 3, "Addition with bilvector should add to x^0 term");
  TEST_ASSERT(sum.getCoefficient(2, {0, 0}) == 4, "Addition with bilvector should add to x^0 term");
  TEST_ASSERT(sum.getCoefficient(1, {1, 0}) == 2, "Original terms should be preserved");

  // Test multiplication with bilvector
  BMPoly product = poly * qpoly;
  TEST_ASSERT(product.getCoefficient(1, {1, 0}) == 6, "q^1*x1^1 * q^0 = 6*q^1*x1^1");
  TEST_ASSERT(product.getCoefficient(3, {1, 0}) == 8, "q^1*x1^1 * q^2 = 8*q^3*x1^1");

  return true;
}

// Test variable inversion
bool test_variable_inversion() {
  BMPoly poly(2, 5);

  // Set up polynomial: 3*q^1*x1^2*x2^1 + 2*q^0*x1^0*x2^2
  poly.setCoefficient(1, {2, 1}, 3);
  poly.setCoefficient(0, {0, 2}, 2);

  // Invert variable 0 (x1)
  BMPoly inverted = poly.invertVariable(0);

  // Should get: 3*q^1*x1^(-2)*x2^1 + 2*q^0*x1^0*x2^2
  TEST_ASSERT(inverted.getCoefficient(1, {-2, 1}) == 3, "Inverted variable coefficient should be 3");
  TEST_ASSERT(inverted.getCoefficient(0, {0, 2}) == 2, "Non-inverted terms should be preserved");

  return true;
}

// Test polynomial truncation
bool test_truncation() {
  BMPoly poly(2, 10);

  // Set up polynomial with various degrees
  poly.setCoefficient(1, {0, 0}, 1);  // constant term
  poly.setCoefficient(1, {2, 1}, 2);  // degree 3 term
  poly.setCoefficient(1, {5, 2}, 3);  // degree 7 term
  poly.setCoefficient(1, {1, 4}, 4);  // degree 5 term

  // Truncate to max degrees [3, 3]
  BMPoly truncated = poly.truncate({3, 3});

  TEST_ASSERT(truncated.getCoefficient(1, {0, 0}) == 1, "Low degree terms should be preserved");
  TEST_ASSERT(truncated.getCoefficient(1, {2, 1}) == 2, "Medium degree terms should be preserved");
  TEST_ASSERT(truncated.getCoefficient(1, {5, 2}) == 0, "High degree terms should be removed");
  TEST_ASSERT(truncated.getCoefficient(1, {1, 4}) == 0, "High degree terms should be removed");

  return true;
}

// Test polynomial evaluation
bool test_evaluation() {
  BMPoly poly(2, 5);

  // Set up polynomial: 2*q^1*x1^1 + 3*q^0*x2^1 + 4*q^2*x1^1*x2^1
  poly.setCoefficient(1, {1, 0}, 2);
  poly.setCoefficient(0, {0, 1}, 3);
  poly.setCoefficient(2, {1, 1}, 4);

  // Evaluate at point [2, 3]
  std::vector<int> point = {2, 3};
  bilvector<int> result = poly.evaluate(point);

  // Expected: 2*q^1*2 + 3*q^0*3 + 4*q^2*2*3 = 4*q^1 + 9*q^0 + 24*q^2
  TEST_ASSERT(result[0] == 9, "q^0 coefficient should be 9");
  TEST_ASSERT(result[1] == 4, "q^1 coefficient should be 4");
  TEST_ASSERT(result[2] == 24, "q^2 coefficient should be 24");

  return true;
}

// Test compatibility checking
bool test_compatibility() {
  BMPoly poly1(2, 5);
  BMPoly poly2(3, 5);  // Different number of variables

  bool caught_exception = false;
  try {
    poly1 += poly2;  // Should throw
  } catch (const std::invalid_argument &) {
    caught_exception = true;
  }

  TEST_ASSERT(caught_exception, "Should throw exception for incompatible polynomials");

  return true;
}

// Test JSON export
bool test_json_export() {
  BMPoly poly(2, 5);

  poly.setCoefficient(1, {1, 0}, 2);
  poly.setCoefficient(0, {0, 1}, 3);
  poly.setCoefficient(-1, {1, 1}, 4);

  // Test that export doesn't crash
  try {
    poly.exportToJson("test_bmpoly_output.json");
  } catch (...) {
    TEST_ASSERT(false, "JSON export should not throw exception");
  }

  return true;
}

// Test constructor with source polynomial
bool test_copy_constructor() {
  BMPoly source(1, 5);
  source.setCoefficient(2, {3}, 7);
  source.setCoefficient(-1, {-2}, 5);

  // Create new polynomial mapping source to variable 1
  BMPoly target(source, 3, 1, 5);

  TEST_ASSERT(target.getNumXVariables() == 3, "Target should have 3 variables");
  TEST_ASSERT(target.getCoefficient(2, {0, 3, 0}) == 7, "Coefficient should be mapped correctly");
  TEST_ASSERT(target.getCoefficient(-1, {0, -2, 0}) == 5, "Negative exponent should be mapped correctly");

  return true;
}

// Test clear and isZero
bool test_clear_and_is_zero() {
  BMPoly poly(2, 5);

  TEST_ASSERT(poly.isZero(), "New polynomial should be zero");

  poly.setCoefficient(1, {1, 0}, 5);
  TEST_ASSERT(!poly.isZero(), "Polynomial with coefficients should not be zero");

  poly.clear();
  TEST_ASSERT(poly.isZero(), "Cleared polynomial should be zero");

  return true;
}

// Test getCoefficients and syncFromDenseVector
bool test_dense_vector_interface() {
  BMPoly poly(2, 3);

  poly.setCoefficient(1, {1, 0}, 5);
  poly.setCoefficient(2, {0, 1}, 3);

  // Get dense vector representation
  auto coeffs = poly.getCoefficients();
  TEST_ASSERT(coeffs.size() > 0, "Dense vector should not be empty");

  // Modify and sync back
  poly.syncFromDenseVector(coeffs);

  TEST_ASSERT(poly.getCoefficient(1, {1, 0}) == 5, "Coefficient should be preserved after sync");
  TEST_ASSERT(poly.getCoefficient(2, {0, 1}) == 3, "Coefficient should be preserved after sync");

  return true;
}

int main() {
  std::cout << "=== BMPoly Unit Tests ===\n\n";

  int total = 0, passed = 0, failed = 0;

  RUN_TEST(test_basic_constructor);
  RUN_TEST(test_constructor_with_max_degrees);
  RUN_TEST(test_coefficient_operations);
  RUN_TEST(test_negative_exponents);
  RUN_TEST(test_q_polynomial_access);
  RUN_TEST(test_addition);
  RUN_TEST(test_subtraction);
  RUN_TEST(test_multiplication);
  RUN_TEST(test_bilvector_operations);
  RUN_TEST(test_variable_inversion);
  RUN_TEST(test_truncation);
  RUN_TEST(test_evaluation);
  RUN_TEST(test_compatibility);
  RUN_TEST(test_json_export);
  RUN_TEST(test_copy_constructor);
  RUN_TEST(test_clear_and_is_zero);
  RUN_TEST(test_dense_vector_interface);

  std::cout << "\n=== Test Results ===\n";
  std::cout << "Total: " << total << std::endl;
  std::cout << "Passed: " << passed << std::endl;
  std::cout << "Failed: " << failed << std::endl;

  if (failed == 0) {
    std::cout << "\n✓ All tests passed!" << std::endl;
    return 0;
  } else {
    std::cout << "\n✗ Some tests failed!" << std::endl;
    return 1;
  }
}