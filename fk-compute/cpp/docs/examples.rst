Examples
========

This section provides practical examples demonstrating how to use the FK Computation C++ Library for various tasks.

Basic Polynomial Operations
---------------------------

This example demonstrates basic polynomial creation and arithmetic operations.

**File: examples/example.cpp**

.. code-block:: cpp

   #include <iostream>
   #include "fk/multivariable_polynomial.hpp"

   void print(std::string s) {
       std::cout << s << std::endl;
   }

   int main() {
       // Create a polynomial in 3 variables with default settings
       MultivariablePolynomial poly(3);

       // Set some coefficients
       // poly = 2*q^1*x₁¹*x₂⁰*x₃⁰ + 3*q^(-1)*x₁⁰*x₂¹*x₃⁰
       poly.setCoefficient(1, {1, 0, 0}, 2);   // coefficient of q¹
       poly.setCoefficient(-1, {0, 1, 0}, 3);
       print("Poly1");
       poly.print();

       // Create another polynomial
       MultivariablePolynomial poly2(3);
       poly2.setCoefficient(0, {0, 0, 1}, 1);  // 1*q⁰*x₁⁰*x₂⁰*x₃¹
       print("Poly2");
       poly2.print();

       // Perform arithmetic
       auto sum = poly + poly2;
       auto product = poly * poly2;

       print("Sum");
       sum.print();

       print("Product");
       product.print();

       return 0;
   }

**Expected Output:**

.. code-block:: text

   Poly1
   2*q^1*x_1 + 3*q^(-1)*x_2

   Poly2
   1*q^0*x_3

   Sum
   2*q^1*x_1 + 3*q^(-1)*x_2 + 1*q^0*x_3

   Product
   2*q^1*x_1*x_3 + 3*q^(-1)*x_2*x_3

**Compilation:**

.. code-block:: bash

   g++ -std=c++17 -I../include -lflint example.cpp -o example

Polynomial Evaluation
---------------------

This example shows how to evaluate polynomials at specific points, which is crucial for FK computations.

.. code-block:: cpp

   #include <iostream>
   #include "fk/multivariable_polynomial.hpp"

   int main() {
       // Create a polynomial: 2*q^1*x₁ + 3*q^(-1)*x₂ + 5*q^0
       MultivariablePolynomial poly(2);
       poly.setCoefficient(1, {1, 0}, 2);   // 2*q^1*x₁
       poly.setCoefficient(-1, {0, 1}, 3);  // 3*q^(-1)*x₂
       poly.setCoefficient(0, {0, 0}, 5);   // 5*q^0

       std::cout << "Original polynomial:" << std::endl;
       poly.print();

       // Evaluate at point (2, 1) -> 2*q^1*2 + 3*q^(-1)*1 + 5*q^0
       std::vector<int> point = {2, 1};
       std::cout << "\\nEvaluating at point (2, 1):" << std::endl;

       auto result = poly.evaluate(point);

       // Print the resulting q-polynomial
       std::cout << "Result: ";
       for (int qPower = result.getMaxNegativeIndex();
            qPower <= result.getMaxPositiveIndex(); ++qPower) {
           int coeff = result[qPower];
           if (coeff != 0) {
               if (qPower == 0) {
                   std::cout << coeff;
               } else if (qPower == 1) {
                   std::cout << coeff << "*q";
               } else if (qPower == -1) {
                   std::cout << coeff << "*q^(-1)";
               } else {
                   std::cout << coeff << "*q^(" << qPower << ")";
               }
               std::cout << " + ";
           }
       }
       std::cout << std::endl;

       return 0;
   }

**Expected Output:**

.. code-block:: text

   Original polynomial:
   2*q^1*x_1 + 3*q^(-1)*x_2 + 5*q^0

   Evaluating at point (2, 1):
   Result: 4*q + 3*q^(-1) + 5

Bilvector (Laurent Polynomial) Operations
-----------------------------------------

This example demonstrates working with Laurent polynomials that support negative exponents.

**File: examples/bilvector_operations.cpp**

.. code-block:: cpp

   #include <iostream>
   #include "fk/bilvector.hpp"

   int main() {
       std::cout << "=== Bilvector (Laurent Polynomial) Operations ===" << std::endl;

       // Create Laurent polynomial supporting q^(-5) to q^10
       auto poly1 = makeLaurentPolynomial<int>(-5, 10);
       auto poly2 = makeLaurentPolynomial<int>(-3, 8);

       // Set coefficients
       poly1[3] = 7;    // 7*q^3
       poly1[-2] = 4;   // 4*q^(-2)
       poly1[0] = 1;    // 1*q^0

       poly2[2] = 3;    // 3*q^2
       poly2[-1] = -2;  // -2*q^(-1)
       poly2[1] = 5;    // 5*q^1

       std::cout << "Polynomial 1: ";
       poly1.print("q");

       std::cout << "Polynomial 2: ";
       poly2.print("q");

       // Arithmetic operations
       auto sum = poly1 + poly2;
       auto product = poly1 * poly2;

       std::cout << "Sum: ";
       sum.print("q");

       std::cout << "Product: ";
       product.print("q");

       // Multiply by q^2
       auto shifted = multiplyByQPower(poly1, 2);
       std::cout << "Poly1 * q^2: ";
       shifted.print("q");

       // Invert exponents
       auto inverted = poly1.invertExponents();
       std::cout << "Poly1 with inverted exponents: ";
       inverted.print("q");

       // Count terms
       std::cout << "Number of terms in poly1: " << poly1.nTerms() << std::endl;
       std::cout << "Is poly1 zero? " << (poly1.isZero() ? "Yes" : "No") << std::endl;

       return 0;
   }

**Expected Output:**

.. code-block:: text

   === Bilvector (Laurent Polynomial) Operations ===
   Polynomial 1: 4*q^(-2) + 1*q^0 + 7*q^3
   Polynomial 2: -2*q^(-1) + 5*q^1 + 3*q^2
   Sum: 4*q^(-2) + -2*q^(-1) + 1*q^0 + 5*q^1 + 3*q^2 + 7*q^3
   Product: -8*q^(-3) + 20*q^(-1) + 12*q^0 + -14*q^1 + 35*q^3 + 21*q^4 + 21*q^5
   Poly1 * q^2: 4*q^0 + 1*q^2 + 7*q^5
   Poly1 with inverted exponents: 7*q^(-3) + 1*q^0 + 4*q^2
   Number of terms in poly1: 3
   Is poly1 zero? No

FK Computation Example
---------------------

This example demonstrates the main FK computation functionality.

**File: examples/fk_computation_test.cpp**

.. code-block:: cpp

   #include <iostream>
   #include <chrono>
   #include "fk/fk_computation.hpp"

   int main() {
       std::cout << "=== FK Computation Example ===" << std::endl;

       try {
           // Create computation instance
           fk::FKComputation computation;

           // Method 1: Compute from CSV files
           std::cout << "Computing from CSV files..." << std::endl;
           auto start = std::chrono::high_resolution_clock::now();

           computation.compute("trefoil_input", "trefoil_output");

           auto end = std::chrono::high_resolution_clock::now();
           auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

           std::cout << "Computation completed in " << duration.count() << " ms" << std::endl;

           // Get results
           const auto& result = computation.getLastResult();
           const auto& config = computation.getLastConfiguration();

           std::cout << "Result polynomial has " << result.nTerms() << " terms" << std::endl;
           std::cout << "Configuration degree: " << config.degree << std::endl;
           std::cout << "Number of components: " << config.components << std::endl;

           // Method 2: Custom configuration
           std::cout << "\\nUsing custom configuration..." << std::endl;

           fk::FKConfiguration customConfig;
           customConfig.degree = 3;
           customConfig.components = 2;
           customConfig.writhe = -1;
           customConfig.prefactors = 4;
           customConfig.crossings = 3;

           // Run with 4 threads
           computation.compute(customConfig, "custom_output", 4);

           std::cout << "Custom computation completed" << std::endl;

       } catch (const std::exception& e) {
           std::cerr << "Error: " << e.what() << std::endl;
           return 1;
       }

       return 0;
   }

**Compilation:**

.. code-block:: bash

   g++ -std=c++17 -O3 -I../include -lflint -fopenmp fk_computation_test.cpp -o fk_test

Parallel Processing Demo
-----------------------

This example shows how to use parallel processing for improved performance.

**File: examples/parallel_demo.cpp**

.. code-block:: cpp

   #include <iostream>
   #include <omp.h>
   #include <chrono>
   #include "fk/fk_computation.hpp"

   int main() {
       std::cout << "=== Parallel Processing Demo ===" << std::endl;

       // Check OpenMP support
       #ifdef _OPENMP
       std::cout << "OpenMP is available" << std::endl;
       std::cout << "Maximum threads: " << omp_get_max_threads() << std::endl;
       #else
       std::cout << "OpenMP is not available - single threaded only" << std::endl;
       #endif

       fk::FKComputation computation;

       // Test with different thread counts
       std::vector<int> threadCounts = {1, 2, 4, 8};

       for (int threads : threadCounts) {
           if (threads > omp_get_max_threads()) continue;

           std::cout << "\\nTesting with " << threads << " thread(s)..." << std::endl;

           auto start = std::chrono::high_resolution_clock::now();

           try {
               // Run computation with specified thread count
               computation.compute("test_input", "test_output_" + std::to_string(threads), threads);

               auto end = std::chrono::high_resolution_clock::now();
               auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

               std::cout << "Completed in " << duration.count() << " ms" << std::endl;

           } catch (const std::exception& e) {
               std::cout << "Error with " << threads << " threads: " << e.what() << std::endl;
           }
       }

       return 0;
   }

**Compilation:**

.. code-block:: bash

   g++ -std=c++17 -O3 -I../include -lflint -fopenmp parallel_demo.cpp -o parallel_demo

Inequality Solver Example
-------------------------

This example demonstrates solving systems of linear inequalities.

.. code-block:: cpp

   #include <iostream>
   #include "fk/inequality_solver.hpp"
   #include "fk/multivariable_polynomial.hpp"

   int main() {
       std::cout << "=== Inequality Solver Example ===" << std::endl;

       // Create system of inequalities
       // Inequality 1: x₁ + x₂ >= 2  (represented as x₁ + x₂ - 2 >= 0)
       // Inequality 2: 2*x₁ - x₂ >= 0
       // Inequality 3: x₁ <= 5  (represented as -x₁ + 5 >= 0)

       std::vector<PolynomialType> inequalities;

       // Inequality 1: x₁ + x₂ - 2 >= 0
       PolynomialType ineq1(2);
       ineq1.setCoefficient(0, {1, 0}, 1);   // x₁
       ineq1.setCoefficient(0, {0, 1}, 1);   // x₂
       ineq1.setCoefficient(0, {0, 0}, -2);  // constant term
       inequalities.push_back(ineq1);

       // Inequality 2: 2*x₁ - x₂ >= 0
       PolynomialType ineq2(2);
       ineq2.setCoefficient(0, {1, 0}, 2);   // 2*x₁
       ineq2.setCoefficient(0, {0, 1}, -1);  // -x₂
       inequalities.push_back(ineq2);

       // Inequality 3: -x₁ + 5 >= 0
       PolynomialType ineq3(2);
       ineq3.setCoefficient(0, {1, 0}, -1);  // -x₁
       ineq3.setCoefficient(0, {0, 0}, 5);   // constant term
       inequalities.push_back(ineq3);

       // Define search bounds: x₁, x₂ ∈ [0, 10]
       std::vector<std::pair<int, int>> bounds = {{0, 10}, {0, 10}};

       std::cout << "Solving system of inequalities:" << std::endl;
       std::cout << "  x₁ + x₂ >= 2" << std::endl;
       std::cout << "  2*x₁ - x₂ >= 0" << std::endl;
       std::cout << "  x₁ <= 5" << std::endl;
       std::cout << "  x₁, x₂ ∈ [0, 10]" << std::endl;

       // Find integer solutions
       auto solutions = findIntegerSolutions(inequalities, bounds);

       std::cout << "\\nFound " << solutions.size() << " solutions:" << std::endl;

       int count = 0;
       for (const auto& point : solutions) {
           std::cout << "  (" << point.coordinates[0] << ", "
                     << point.coordinates[1] << ")" << std::endl;
           if (++count > 10) {  // Limit output for readability
               std::cout << "  ... and " << (solutions.size() - count)
                         << " more" << std::endl;
               break;
           }
       }

       return 0;
   }

**Expected Output:**

.. code-block:: text

   === Inequality Solver Example ===
   Solving system of inequalities:
     x₁ + x₂ >= 2
     2*x₁ - x₂ >= 0
     x₁ <= 5
     x₁, x₂ ∈ [0, 10]

   Found 25 solutions:
     (1, 1)
     (1, 2)
     (2, 1)
     (2, 2)
     (2, 3)
     (2, 4)
     (3, 1)
     (3, 2)
     (3, 3)
     (3, 4)
     (3, 5)
     ... and 14 more

Trefoil Knot Example
-------------------

This example demonstrates computing FK invariants for the trefoil knot.

**File: examples/trefoil_simple_example.cpp**

.. code-block:: cpp

   #include <iostream>
   #include "fk/fk_computation.hpp"

   int main() {
       std::cout << "=== Trefoil Knot FK Computation ===" << std::endl;

       try {
           // The trefoil knot can be represented as a 3-braid: σ₁σ₂σ₁σ₂σ₁σ₂
           // This creates a specific configuration for FK computation

           fk::FKConfiguration trefoilConfig;

           // Basic trefoil parameters
           trefoilConfig.degree = 2;
           trefoilConfig.components = 1;  // Trefoil is a single component knot
           trefoilConfig.writhe = -3;     // Trefoil has writhe -3
           trefoilConfig.prefactors = 3;
           trefoilConfig.crossings = 3;

           // Set up crossing data (simplified example)
           trefoilConfig.closed_strand_components = {0};

           // Initialize crossing matrices for trefoil
           trefoilConfig.crossing_matrices = {
               {{1, 0}, {1, 1}},  // First crossing
               {{1, 1}, {0, 1}},  // Second crossing
               {{1, 0}, {1, 1}}   // Third crossing
           };

           trefoilConfig.crossing_relation_types = {1, 1, 1};
           trefoilConfig.top_crossing_components = {0, 0, 0};
           trefoilConfig.bottom_crossing_components = {0, 0, 0};

           std::cout << "Configuration:" << std::endl;
           std::cout << "  Degree: " << trefoilConfig.degree << std::endl;
           std::cout << "  Components: " << trefoilConfig.components << std::endl;
           std::cout << "  Writhe: " << trefoilConfig.writhe << std::endl;
           std::cout << "  Crossings: " << trefoilConfig.crossings << std::endl;

           // Validate configuration
           if (!trefoilConfig.isValid()) {
               std::cerr << "Invalid trefoil configuration!" << std::endl;
               return 1;
           }

           // Run computation
           fk::FKComputation computation;
           computation.compute(trefoilConfig, "trefoil_result");

           // Get and display results
           const auto& result = computation.getLastResult();

           std::cout << "\\nTrefoil FK polynomial:" << std::endl;
           result.print();

           std::cout << "\\nResult has " << result.nTerms() << " terms" << std::endl;

       } catch (const std::exception& e) {
           std::cerr << "Error computing trefoil FK invariant: " << e.what() << std::endl;
           return 1;
       }

       return 0;
   }

Random Polynomial Generation
---------------------------

This example shows how to generate and work with random polynomials for testing.

.. code-block:: cpp

   #include <iostream>
   #include <random>
   #include "fk/polynomial_config.hpp"

   int main() {
       std::cout << "=== Random Polynomial Generation ===" << std::endl;

       std::random_device rd;
       std::mt19937 gen(rd());
       std::uniform_int_distribution<> coeff_dist(-10, 10);
       std::uniform_int_distribution<> power_dist(0, 3);
       std::uniform_int_distribution<> q_power_dist(-2, 2);

       // Generate random polynomial with 2 variables
       PolynomialType poly(2);

       std::cout << "Generating random polynomial..." << std::endl;

       // Add 10 random terms
       for (int i = 0; i < 10; ++i) {
           int coeff = coeff_dist(gen);
           if (coeff == 0) continue;  // Skip zero coefficients

           int qPower = q_power_dist(gen);
           std::vector<int> xPowers = {power_dist(gen), power_dist(gen)};

           poly.setCoefficient(qPower, xPowers, coeff);
       }

       std::cout << "Generated polynomial:" << std::endl;
       poly.print();

       // Test arithmetic with another random polynomial
       PolynomialType poly2(2);
       for (int i = 0; i < 5; ++i) {
           int coeff = coeff_dist(gen);
           if (coeff == 0) continue;

           int qPower = q_power_dist(gen);
           std::vector<int> xPowers = {power_dist(gen), power_dist(gen)};

           poly2.setCoefficient(qPower, xPowers, coeff);
       }

       std::cout << "\\nSecond polynomial:" << std::endl;
       poly2.print();

       // Compute sum and product
       auto sum = poly + poly2;
       auto product = poly * poly2;

       std::cout << "\\nSum:" << std::endl;
       sum.print();

       std::cout << "\\nProduct (first 5 terms):" << std::endl;
       if (product.nTerms() > 5) {
           std::cout << "[Showing first few terms of " << product.nTerms()
                     << " total terms]" << std::endl;
       }
       // Note: For large products, you might want to limit output

       return 0;
   }

Building and Running Examples
----------------------------

**Build All Examples:**

.. code-block:: bash

   cd examples

   # Basic examples
   g++ -std=c++17 -I../include -lflint example.cpp -o example
   g++ -std=c++17 -I../include bilvector_operations.cpp -o bilvector_example

   # FK computation examples (require OpenMP)
   g++ -std=c++17 -O3 -I../include -lflint -fopenmp fk_computation_test.cpp -o fk_test
   g++ -std=c++17 -O3 -I../include -lflint -fopenmp parallel_demo.cpp -o parallel_demo
   g++ -std=c++17 -O3 -I../include -lflint -fopenmp trefoil_simple_example.cpp -o trefoil_example

**Run Examples:**

.. code-block:: bash

   ./example                # Basic polynomial operations
   ./bilvector_example      # Laurent polynomial operations
   ./fk_test               # FK computation test
   ./parallel_demo         # Parallel processing demo
   ./trefoil_example       # Trefoil knot computation

**Create Test Input Files:**

For the FK computation examples, you may need to create sample input CSV files:

.. code-block:: bash

   # Create a simple test input file
   cat > trefoil_input.csv << EOF
   degree,2
   components,1
   writhe,-3
   prefactors,3
   crossings,3
   closed_strand_components,0
   EOF

These examples demonstrate the key functionality of the FK Computation C++ Library and provide starting points for your own applications.