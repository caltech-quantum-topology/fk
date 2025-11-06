Quick Start Guide
=================

This guide will get you up and running with the FK Computation library in just a few minutes.

.. contents:: Table of Contents
   :local:
   :depth: 2

First Steps
-----------

After :doc:`installation`, verify everything works:

.. code-block:: bash

   # Build the library
   make all

   # Run a quick test
   make run-arithmetic

   # Try parallel processing
   ./parallel_demo --help

Basic Polynomial Operations
--------------------------

Let's start with basic polynomial operations:

.. code-block:: cpp

  #include <iostream>
  #include "fk/multivariable_polynomial.hpp"

  void print(std::string s){
      std::cout << s << std::endl;
  }

  int main() {
      // Create a polynomial in 3 variables with default settings
      MultivariablePolynomial poly(3);


      // Set some coefficients
      // poly = 2*q^1*x₁¹*x₂⁰*x₃⁰ + 3*q^(-1)*x₁⁰*x₂¹*x₃⁰
      poly.setCoefficient(1,{1, 0, 0},2);   // coefficient of q¹
      poly.setCoefficient(-1,{0, 1, 0},3);
      print("Poly1");
      poly.print();

      // Create another polynomial
      MultivariablePolynomial poly2(3);
      poly2.setCoefficient(0,{0, 0, 1},1);  // 1*q⁰*x₁⁰*x₂⁰*x₃¹
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

Compile and run:

.. code-block:: bash

   g++ -std=c++17 -Iinclude example.cpp src/multivariable_polynomial.cpp -o example
   ./example

Simple Solution Finding
-----------------------

Find integer solutions to a constraint system:

.. code-block:: cpp

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

Expected output:

.. code-block:: text

     === Linear Inequality Solver Example ===

  Example 1: 2D system
  Inequalities:
    x + y >= 1
    x - y >= -2
    x >= 0
    y >= 0
  Search bounds: x ∈ [0, 5], y ∈ [0, 5]

  Found 29 integer solutions:
  (0, 1)
  (0, 2)
  (1, 0)
  (1, 1)
  (1, 2)
  (1, 3)
  (2, 0)
  (2, 1)
  (2, 2)
  (2, 3)
  (2, 4)
  (3, 0)
  (3, 1)
  (3, 2)
  (3, 3)
  (3, 4)
  (3, 5)
  (4, 0)
  (4, 1)
  (4, 2)
  (4, 3)
  (4, 4)
  (4, 5)
  (5, 0)
  (5, 1)
  (5, 2)
  (5, 3)
  (5, 4)
  (5, 5)

  Example 2: 3D system
  Inequalities:
    x + y + z >= 2
    x - y >= 0
    z >= 1
  Search bounds: x ∈ [0, 3], y ∈ [0, 3], z ∈ [0, 3]

  Found 29 integer solutions:
  (0, 0, 2)
  (0, 0, 3)
  (1, 0, 1)
  (1, 0, 2)
  (1, 0, 3)
  (1, 1, 1)
  (1, 1, 2)
  (1, 1, 3)
  (2, 0, 1)
  (2, 0, 2)
  (2, 0, 3)
  (2, 1, 1)
  (2, 1, 2)
  (2, 1, 3)
  (2, 2, 1)
  (2, 2, 2)
  (2, 2, 3)
  (3, 0, 1)
  (3, 0, 2)
  (3, 0, 3)
  (3, 1, 1)
  (3, 1, 2)
  (3, 1, 3)
  (3, 2, 1)
  (3, 2, 2)
  (3, 2, 3)
  (3, 3, 1)
  (3, 3, 2)
  (3, 3, 3)

  Example 3: Triangle region
  Inequalities:
    x + y <= 4
    x >= 1
    y >= 1
  Search bounds: x ∈ [0, 5], y ∈ [0, 5]

  Found 6 integer solutions:
  (1, 1)
  (1, 2)
  (1, 3)
  (2, 1)
  (2, 2)
  (3, 1)


Parallel Processing
------------------

Speed up computations with parallel processing:

.. code-block:: cpp

   #include "fk/parallel_pool.hpp"
   #include <chrono>

   int main() {
       // Larger constraint system
       std::vector<std::vector<double>> main_inequalities = {
           {1, 1, 1, -10},  // x + y + z ≤ 10
           {1, 0, 0, 0},    // x ≥ 0
           {0, 1, 0, 0},    // y ≥ 0
           {0, 0, 1, 0}     // z ≥ 0
       };

       std::vector<std::vector<double>> supporting_inequalities = {
           {2, 1, 1, -15}   // 2x + y + z ≤ 15
       };

       int solution_count = 0;
       auto callback = [&solution_count](const std::vector<int>& solution) {
           solution_count++;
           if (solution_count % 100 == 0) {
               std::cout << "Found " << solution_count << " solutions..." << std::endl;
           }
       };

       // Time the computation
       auto start = std::chrono::steady_clock::now();

       // Run with automatic thread detection
       parallelPooling(main_inequalities, supporting_inequalities, callback);

       auto end = std::chrono::steady_clock::now();
       auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

       std::cout << "Found " << solution_count << " solutions in "
                 << duration.count() << "ms" << std::endl;

       return 0;
   }

Checkpointed Computation
-----------------------

For long-running computations, enable checkpointing:

.. code-block:: cpp

   #include "fk/parallel_pool.hpp"
   #include <csignal>

   // Global state for signal handling
   bool interrupted = false;

   void signalHandler(int signal) {
       std::cout << "\\nReceived interrupt signal, saving checkpoint..." << std::endl;
       interrupted = true;
   }

   int main() {
       // Setup signal handling
       std::signal(SIGINT, signalHandler);

       // Large problem definition
       std::vector<std::vector<double>> main_inequalities = {
           // ... large constraint system
       };

       std::vector<std::vector<double>> supporting_inequalities = {
           // ... additional constraints
       };

       auto callback = [](const std::vector<int>& solution) {
           // Process each solution
       };

       try {
           // Run with checkpointing every 1000 solutions
           parallelPoolingWithCheckpoints(
               main_inequalities,
               supporting_inequalities,
               callback,
               "computation.ckpt",  // checkpoint file
               true,                // resume if exists
               1000,                // checkpoint interval
               0                    // auto-detect threads
           );
       } catch (const std::exception& e) {
           std::cout << "Computation interrupted: " << e.what() << std::endl;
           std::cout << "Progress saved to checkpoint file" << std::endl;
       }

       return 0;
   }

To resume a computation:

.. code-block:: bash

   # The program will automatically resume from checkpoint.ckpt
   ./your_program

Advanced Data Structures
------------------------

The library provides specialized data structures:

**Bilvector (Bidirectional Vector)**

.. code-block:: cpp

   #include "fk/bilvector.hpp"

   // Create bilvector: 5 negative, 5 positive, component size 1, default 0
   bilvector<int> bv(5, 5, 1, 0);

   // Access with positive and negative indices
   bv[10] = 42;    // Automatically expands
   bv[-5] = 7;     // Negative indexing
   bv[0] = 1;      // Zero goes to positive side

   std::cout << bv[-5] << ", " << bv[0] << ", " << bv[10] << std::endl;
   // Output: 7, 1, 42

**Binary Tree with Substring Matching**

.. code-block:: cpp

   #include "fk/btree.hpp"

   btree<int> tree;

   // Insert vectors
   tree.insertVector({1, 2, 3, 4});
   tree.insertVector({5, 6, 7});

   // Check containment (supports substring matching)
   bool exact = tree.containsVector({1, 2, 3, 4});  // true
   bool substring = tree.containsVector({3, 4});     // true (substring)
   bool not_found = tree.containsVector({8, 9});     // false

Performance Tips
---------------

**1. Problem Sizing**

- Small problems (< 100 solutions): Use sequential version
- Medium problems (100-10,000 solutions): Use 2-8 threads
- Large problems (> 10,000 solutions): Use all available cores

**2. Memory Management**

.. code-block:: cpp

   // For large polynomials, clear when done
   polynomial.clear();
   polynomial.pruneZeros();  // Remove zero coefficients

**3. Thread Configuration**

.. code-block:: cpp

   // Detect optimal thread count
   int threads = std::thread::hardware_concurrency();
   if (expected_solutions < 1000) threads = 1;

   parallelPooling(main_ineq, supp_ineq, callback, threads);

**4. Checkpointing Strategy**

.. code-block:: cpp

   // Checkpoint frequency based on computation time
   size_t checkpoint_interval;
   if (expected_runtime_hours > 1) {
       checkpoint_interval = 100;   // Frequent checkpoints
   } else {
       checkpoint_interval = 1000;  // Less frequent
   }

Common Patterns
--------------

**Pattern 1: Batch Processing**

.. code-block:: cpp

   std::vector<std::vector<int>> solutions;
   auto collector = [&solutions](const std::vector<int>& sol) {
       solutions.push_back(sol);
   };

   parallelPooling(main_ineq, supp_ineq, collector);

   // Process solutions in batches
   for (size_t i = 0; i < solutions.size(); i += 1000) {
       processBatch(solutions.begin() + i,
                   solutions.begin() + std::min(i + 1000, solutions.size()));
   }

**Pattern 2: Progress Monitoring**

.. code-block:: cpp

   #include <atomic>

   std::atomic<size_t> progress{0};
   auto callback = [&progress](const std::vector<int>& solution) {
       if (++progress % 1000 == 0) {
           std::cout << "Progress: " << progress << " solutions" << std::endl;
       }
   };

**Pattern 3: Filtering Solutions**

.. code-block:: cpp

   auto filter_callback = [](const std::vector<int>& solution) {
       // Only process solutions meeting additional criteria
       int sum = std::accumulate(solution.begin(), solution.end(), 0);
       if (sum % 2 == 0) {  // Even sum
           processSolution(solution);
       }
   };

Next Steps
----------

Now that you're familiar with the basics:

1. Explore :doc:`tutorials/index` for detailed examples
2. Read :doc:`api/polynomials` for complete API reference
3. Check :doc:`performance/index` for optimization techniques
4. See :doc:`examples/index` for real-world applications

Common Issues
-------------

**Thread count seems wrong?**

.. code-block:: bash

   # Check actual CPU cores
   nproc
   cat /proc/cpuinfo | grep "cpu cores"

**Memory usage too high?**

.. code-block:: cpp

   // Use checkpointing for large problems
   parallelPoolingWithCheckpoints(/* ... */, "checkpoint.dat", true, 500);

**Compilation errors?**

.. code-block:: bash

   # Ensure C++17 support
   g++ --version
   # Should be 7.0+ or equivalent clang

**Need help?**

- Check the :doc:`installation` guide for troubleshooting
- Review :doc:`api/index` for detailed function documentation
- See :doc:`examples/index` for working code samples
