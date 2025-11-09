User Guide
==========

This comprehensive guide covers the key concepts and usage patterns for the FK Computation C++ Library.

Core Concepts
-------------

Polynomial Systems
^^^^^^^^^^^^^^^^^^

The library provides three different polynomial implementations, all with compatible interfaces:

**MultivariablePolynomial** (Default)
   Sparse polynomial representation using hash maps. Best for polynomials with few terms.

.. code-block:: cpp

   #include "fk/multivariable_polynomial.hpp"

   MultivariablePolynomial poly(3);  // 3 variables
   poly.setCoefficient(1, {2, 0, 1}, 5);  // 5*q^1*x₁²*x₃

**FMPoly** (FLINT-based)
   Dense polynomial representation using the FLINT library. Best for performance-critical applications.

.. code-block:: cpp

   #include "fk/fmpoly.hpp"

   FMPoly poly(3);  // 3 variables
   poly.setCoefficient(1, {2, 0, 1}, 5);  // Same interface

**BMPoly** (Basic implementation)
   Simple vector-based implementation. Good for educational purposes and debugging.

.. code-block:: cpp

   #include "fk/bmpoly.hpp"

   BMPoly poly(3);  // 3 variables
   poly.setCoefficient(1, {2, 0, 1}, 5);  // Same interface

Configuration Selection
^^^^^^^^^^^^^^^^^^^^^^^

Choose the polynomial backend at compile time using the configuration header:

.. code-block:: cpp

   #include "fk/polynomial_config.hpp"

   // PolynomialType is automatically set based on POLYNOMIAL_TYPE macro
   PolynomialType poly(3);

Or set the macro directly:

.. code-block:: bash

   g++ -DPOLYNOMIAL_TYPE=1 ...  # Use FMPoly

Laurent Polynomials (Bilvector)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For polynomials with negative exponents, use the bilvector data structure:

.. code-block:: cpp

   #include "fk/bilvector.hpp"

   // Create Laurent polynomial: supports negative powers
   auto poly = makeLaurentPolynomial<int>(-5, 10);  // q^(-5) to q^10

   poly[3] = 7;   // Coefficient of q^3
   poly[-2] = 4;  // Coefficient of q^(-2)

   poly.print("q");  // Print in readable form

Working with Polynomials
-------------------------

Basic Operations
^^^^^^^^^^^^^^^^

All polynomial types support the same basic interface:

.. code-block:: cpp

   #include "fk/polynomial_config.hpp"

   // Create polynomials
   PolynomialType poly1(2);  // 2 variables
   PolynomialType poly2(2);

   // Set coefficients
   poly1.setCoefficient(1, {1, 0}, 3);  // 3*q^1*x₁
   poly1.setCoefficient(0, {0, 1}, 2);  // 2*q^0*x₂

   poly2.setCoefficient(1, {0, 1}, 4);  // 4*q^1*x₂
   poly2.setCoefficient(-1, {1, 0}, 1); // 1*q^(-1)*x₁

   // Arithmetic operations
   auto sum = poly1 + poly2;
   auto product = poly1 * poly2;
   auto difference = poly1 - poly2;

   // Get coefficients
   int coeff = poly1.getCoefficient(1, {1, 0});  // Returns 3

   // Print polynomials
   poly1.print();
   sum.print();

Advanced Operations
^^^^^^^^^^^^^^^^^^^

**Evaluation at Points**

Evaluate polynomials at specific integer points:

.. code-block:: cpp

   PolynomialType poly(2);
   poly.setCoefficient(1, {1, 1}, 2);  // 2*q^1*x₁*x₂
   poly.setCoefficient(0, {0, 0}, 3);  // 3*q^0

   std::vector<int> point = {2, 3};  // x₁=2, x₂=3
   auto result = poly.evaluate(point);  // Returns Laurent polynomial in q

**Term Counting and Analysis**

.. code-block:: cpp

   int numTerms = poly.nTerms();
   bool isEmpty = poly.isZero();

   // Get polynomial properties
   int degree = poly.getDegree();
   int numVars = poly.getNumVariables();

**Serialization**

.. code-block:: cpp

   // Convert to/from vectors for storage
   auto coeffData = poly.getCoefficients();
   auto termData = poly.getTerms();

   // Reconstruct polynomial
   PolynomialType newPoly(numVars);
   newPoly.setFromData(coeffData, termData);

FK Computation Engine
---------------------

The main computation interface provides high-level orchestration of FK calculations:

Basic Usage
^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/fk_computation.hpp"

   // Create computation instance
   fk::FKComputation computation;

   // Run computation from CSV input file
   computation.compute("trefoil_input", "trefoil_output");

   // Get results
   const auto& result = computation.getLastResult();
   const auto& config = computation.getLastConfiguration();

Custom Configuration
^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/fk_computation.hpp"

   // Create custom configuration
   fk::FKConfiguration config;
   config.degree = 3;
   config.components = 2;
   config.writhe = -2;

   // Set up constraint data
   config.inequalities = loadInequalityData();
   config.variable_assignments = generateVariableAssignments();

   // Run computation with custom config
   fk::FKComputation computation;
   computation.compute(config, "output_file", 4);  // 4 threads

Parallel Processing
^^^^^^^^^^^^^^^^^^^

The library supports parallel processing for improved performance:

.. code-block:: cpp

   #include "fk/fk_computation.hpp"
   #include <omp.h>

   // Set OpenMP thread count
   omp_set_num_threads(8);

   // Run parallel computation
   fk::FKComputation computation;
   computation.compute("input", "output", 8);  // Use 8 worker threads

Input/Output Handling
---------------------

File Formats
^^^^^^^^^^^^

**Input CSV Format**

The library expects CSV input files with the following structure:

.. code-block:: text

   degree,3
   components,2
   writhe,-1
   prefactors,4
   crossings,3
   closed_strand_components,1,2
   crossing_matrices,[[1,0],[0,1]],[[1,1],[0,1]]
   # ... additional configuration data

**Output JSON Format**

Results are written in JSON format:

.. code-block:: json

   {
     "degree": 3,
     "components": 2,
     "result_polynomial": {
       "terms": [
         {"q_power": 1, "x_powers": [1, 0], "coefficient": 3},
         {"q_power": 0, "x_powers": [0, 1], "coefficient": -2}
       ]
     },
     "computation_metadata": {
       "num_threads": 4,
       "computation_time": 1.23
     }
   }

Parser Usage
^^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/fk_computation.hpp"

   // Parse input file
   fk::FKInputParser parser;
   fk::FKConfiguration config = parser.parseFromFile("input_data");

   // Validate configuration
   if (!config.isValid()) {
       throw std::runtime_error("Invalid input configuration");
   }

   // Use configuration
   fk::FKComputationEngine engine(config);

Result Writer
^^^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/fk_computation.hpp"

   fk::FKResultWriter writer;

   // Write to JSON
   writer.writeToJson(result, "output.json");

   // Write to human-readable text
   writer.writeToText(result, "output.txt");

Linear Algebra and Utilities
-----------------------------

The library includes mathematical utilities for common operations:

Matrix Operations
^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/linalg.hpp"

   // Matrix-vector operations
   std::vector<std::vector<double>> matrix = {{1, 2}, {3, 4}};
   std::vector<double> vector = {5, 6};

   auto result = matrixVectorMultiply(matrix, vector);

Inequality Solving
^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   #include "fk/inequality_solver.hpp"

   // Define system of inequalities
   std::vector<PolynomialType> inequalities;
   inequalities.push_back(createLinearInequality({1, -1, 2}));  // x₁ - x₂ + 2 >= 0

   // Define search bounds
   std::vector<std::pair<int, int>> bounds = {{0, 10}, {0, 10}};  // x₁, x₂ ∈ [0, 10]

   // Find integer solutions
   auto solutions = findIntegerSolutions(inequalities, bounds);

   for (const auto& point : solutions) {
       // Process each solution point
       std::cout << "Solution: ";
       for (int coord : point.coordinates) {
           std::cout << coord << " ";
       }
       std::cout << std::endl;
   }

Advanced Topics
---------------

Template Specialization
^^^^^^^^^^^^^^^^^^^^^^^^

The library uses templates extensively. You can specialize for custom types:

.. code-block:: cpp

   #include "fk/bilvector.hpp"

   // Specialize for custom coefficient type
   template<>
   struct bilvector<MyCustomType> {
       // Custom implementation
   };

Memory Management
^^^^^^^^^^^^^^^^^

**Large Polynomials**

For very large polynomials, consider:

.. code-block:: cpp

   // Use sparse representation for polynomials with few terms
   #define POLYNOMIAL_TYPE 0  // MultivariablePolynomial

   // Use dense representation for polynomials with many terms
   #define POLYNOMIAL_TYPE 1  // FMPoly

**Memory Monitoring**

.. code-block:: cpp

   #include <chrono>

   auto start = std::chrono::high_resolution_clock::now();

   // Perform computation
   computation.compute("input", "output");

   auto end = std::chrono::high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

Error Handling
^^^^^^^^^^^^^^

The library uses exceptions for error reporting:

.. code-block:: cpp

   try {
       fk::FKComputation computation;
       computation.compute("invalid_input", "output");
   } catch (const std::runtime_error& e) {
       std::cerr << "Computation error: " << e.what() << std::endl;
   } catch (const std::domain_error& e) {
       std::cerr << "Domain error: " << e.what() << std::endl;
   }

Best Practices
--------------

Performance Optimization
^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Choose the Right Backend**: Use FMPoly for dense polynomials, MultivariablePolynomial for sparse ones
2. **Enable Compiler Optimizations**: Always use ``-O3`` for production builds
3. **Use Parallel Processing**: Enable OpenMP for multi-core systems
4. **Profile Your Code**: Use profiling tools to identify bottlenecks

Code Organization
^^^^^^^^^^^^^^^^^

1. **Include Only What You Need**: Don't include the entire library if you only need specific components
2. **Use the Configuration Header**: Use ``polynomial_config.hpp`` for backend-agnostic code
3. **Handle Errors Gracefully**: Always check return values and catch exceptions
4. **Document Your Usage**: Complex polynomial expressions benefit from clear comments

Testing and Validation
^^^^^^^^^^^^^^^^^^^^^^^

1. **Validate Input Data**: Always check configuration validity before computation
2. **Test with Known Results**: Use simple cases with known outcomes for validation
3. **Compare Backends**: Run the same computation with different polynomial backends to verify consistency
4. **Monitor Performance**: Track computation times and memory usage for regression testing

Next Steps
----------

* Explore the :doc:`examples` for practical applications
* Review the :doc:`api/library_root` for detailed API documentation
* Look at the source code for advanced usage patterns
* Consider contributing improvements or reporting issues