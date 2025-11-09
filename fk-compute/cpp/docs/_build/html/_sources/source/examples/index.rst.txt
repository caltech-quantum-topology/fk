Examples and Use Cases
======================

.. contents:: Table of Contents
   :local:
   :depth: 2

This section provides practical examples and real-world use cases for the FK Computation library.

Getting Started Examples
------------------------

Basic Polynomial Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: polynomial_basic.cpp
   :language: cpp
   :caption: Basic polynomial operations with sparse storage

This example demonstrates:

- Creating polynomials with multiple variables
- Setting coefficients for different exponent combinations
- Performing arithmetic operations (+, -, *)
- Accessing results and term counts

Constraint Solving
~~~~~~~~~~~~~~~~~

.. literalinclude:: constraint_solving.cpp
   :language: cpp
   :caption: Finding integer solutions to linear constraints

**Problem**: Find all non-negative integer solutions to:

- x + y ≤ 5
- x ≥ 0, y ≥ 0
- Additional constraint: 2x + y ≤ 7

This example shows how to:

- Define constraint matrices
- Set up solution callbacks
- Enumerate all feasible integer points

Advanced Examples
-----------------

Parallel Processing with Progress Monitoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: parallel_monitoring.cpp
   :language: cpp
   :caption: Parallel computation with real-time progress tracking

Features demonstrated:

- Automatic thread detection
- Real-time progress reporting
- Performance timing and statistics
- Thread-safe solution collection

Fault-Tolerant Long-Running Computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: checkpointed_computation.cpp
   :language: cpp
   :caption: Checkpointed computation with graceful interruption handling

This example covers:

- Signal handling for graceful interruption
- Automatic checkpoint saving
- Resume from saved state
- Progress persistence across runs

Real-World Applications
----------------------

Mathematical Research: Polynomial Families
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: polynomial_families.cpp
   :language: cpp
   :caption: Generating and analyzing polynomial families

**Use Case**: Generate families of polynomials with specific properties and analyze their structure.

Applications:

- Studying polynomial invariants
- Pattern recognition in coefficient sequences
- Mathematical conjecture testing

High-Performance Computing: Large Solution Spaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: large_solution_spaces.cpp
   :language: cpp
   :caption: Efficient enumeration of large solution spaces

**Problem Scale**: Enumerate solutions to systems with:

- 6+ variables
- Complex constraint interactions
- Expected 100,000+ solutions

Techniques demonstrated:

- Memory-efficient solution processing
- Batch processing for large datasets
- Performance optimization strategies

Scientific Computing: Parameter Sweeps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: parameter_sweeps.cpp
   :language: cpp
   :caption: Systematic parameter space exploration

**Research Application**: Systematically explore parameter spaces for mathematical models.

Features:

- Parameterized problem generation
- Statistical analysis of results
- Data export for further analysis

Specialized Use Cases
--------------------

Custom Data Structures
~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: custom_structures.cpp
   :language: cpp
   :caption: Advanced usage of bilvector and btree

Demonstrates:

- Bidirectional indexing with bilvectors
- Substring matching with binary trees
- Custom data organization patterns

Integration with External Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: external_integration.cpp
   :language: cpp
   :caption: Integrating with visualization and analysis tools

Shows how to:

- Export results to common formats (JSON, CSV)
- Interface with plotting libraries
- Prepare data for statistical analysis

Performance Optimization Examples
---------------------------------

Memory-Conscious Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: memory_optimization.cpp
   :language: cpp
   :caption: Optimizing memory usage for large problems

Techniques:

- Streaming solution processing
- Memory pooling for frequent allocations
- Garbage collection strategies

CPU-Intensive Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: cpu_optimization.cpp
   :language: cpp
   :caption: CPU optimization techniques

Covers:

- Thread affinity setting
- Cache-friendly data access patterns
- SIMD optimization opportunities

Debugging and Testing Examples
------------------------------

Correctness Verification
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: correctness_verification.cpp
   :language: cpp
   :caption: Verifying algorithmic correctness

Methods demonstrated:

- Cross-validation between sequential and parallel versions
- Result reproducibility testing
- Mathematical property verification

Performance Regression Testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: regression_testing.cpp
   :language: cpp
   :caption: Automated performance regression detection

Tools and techniques:

- Baseline performance measurement
- Automated regression detection
- Performance trend analysis

Interactive Examples
--------------------

Command-Line Interface
~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: cli_interface.cpp
   :language: cpp
   :caption: Building interactive command-line tools

Features:

- Argument parsing and validation
- Interactive problem setup
- Real-time user feedback

Web Interface Integration
~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: web_interface.cpp
   :language: cpp
   :caption: Creating web-accessible computation services

Technologies:

- REST API endpoints
- JSON request/response handling
- Asynchronous computation management

Educational Examples
--------------------

Algorithm Visualization
~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: algorithm_visualization.cpp
   :language: cpp
   :caption: Visualizing algorithm progress and decision trees

Educational value:

- Step-by-step algorithm execution
- Decision tree exploration
- Visual feedback for learning

Complexity Analysis
~~~~~~~~~~~~~~~~~

.. literalinclude:: complexity_analysis.cpp
   :language: cpp
   :caption: Empirical complexity analysis and measurement

Demonstrates:

- Runtime scaling measurement
- Memory usage profiling
- Algorithmic complexity verification

Example Code Files
------------------

The following code files are referenced in the examples above:

.. toctree::
   :maxdepth: 1
   :caption: Source Code Files

   source_files/polynomial_basic
   source_files/constraint_solving
   source_files/parallel_monitoring
   source_files/checkpointed_computation
   source_files/polynomial_families
   source_files/large_solution_spaces
   source_files/parameter_sweeps
   source_files/custom_structures
   source_files/external_integration
   source_files/memory_optimization
   source_files/cpu_optimization
   source_files/correctness_verification
   source_files/regression_testing
   source_files/cli_interface
   source_files/web_interface
   source_files/algorithm_visualization
   source_files/complexity_analysis

Building and Running Examples
-----------------------------

Prerequisites
~~~~~~~~~~~~

.. code-block:: bash

   # Ensure library is built
   make all

   # Install any additional dependencies
   sudo apt install libcurl4-openssl-dev  # For web examples
   pip install matplotlib numpy           # For Python integration

Compilation
~~~~~~~~~~

Each example can be compiled individually:

.. code-block:: bash

   # Basic examples
   g++ -std=c++17 -Iinclude examples/polynomial_basic.cpp \
       src/multivariable_polynomial.cpp -o polynomial_basic

   # Parallel examples (requires pthread)
   g++ -std=c++17 -Iinclude -pthread examples/parallel_monitoring.cpp \
       src/parallel_pool.cpp src/solution_pool_1a_double_links.cpp \
       -o parallel_monitoring

   # All examples with Makefile
   make examples

Execution
~~~~~~~~

.. code-block:: bash

   # Run basic examples
   ./polynomial_basic

   # Run parallel examples with different configurations
   ./parallel_monitoring --threads 4 --problem medium
   ./checkpointed_computation --checkpoint-interval 1000

   # Performance examples
   ./memory_optimization --max-memory 4GB
   ./cpu_optimization --enable-simd

Example Output
~~~~~~~~~~~~~

Expected output from various examples:

**Polynomial Basic Example:**

.. code-block:: text

   Creating polynomials with 3 variables...
   Polynomial p1: 2 terms
   Polynomial p2: 1 terms
   Sum: 3 terms
   Product: 2 terms

   Coefficient at {1,0,0}, q^1: 2
   Coefficient at {0,1,0}, q^-1: 3

**Parallel Monitoring Example:**

.. code-block:: text

   === Parallel Solution Enumeration ===
   Hardware threads detected: 24
   Using 4 threads for medium problem

   Progress: 1000 solutions (1.2 seconds, 833/sec)
   Progress: 2000 solutions (2.1 seconds, 952/sec)
   Progress: 3000 solutions (3.0 seconds, 1000/sec)

   Computation completed:
   - Total solutions: 3,247
   - Execution time: 3.2 seconds
   - Average rate: 1,015 solutions/second
   - Thread efficiency: 94%

**Checkpointed Computation Example:**

.. code-block:: text

   === Fault-Tolerant Computation ===
   Checkpoint file: computation.ckpt (found existing)
   Resuming from checkpoint: 15,420 solutions processed

   Progress: 16000 solutions...
   Progress: 17000 solutions...
   ^C
   Interrupt signal received!
   Saving checkpoint: 17,389 solutions processed
   Computation can be resumed with: ./checkpointed_computation

Troubleshooting Examples
-----------------------

Common Issues
~~~~~~~~~~~~

**Compilation Errors:**

.. code-block:: bash

   # Missing headers
   error: 'std::thread' has not been declared
   # Solution: Add -pthread flag

   # Undefined references
   undefined reference to `parallelPooling`
   # Solution: Link required .cpp files

**Runtime Issues:**

.. code-block:: bash

   # Segmentation fault
   # Solution: Check array bounds, use debug build
   make debug
   gdb ./example_program

   # Performance issues
   # Solution: Profile and optimize
   perf record ./example_program
   perf report

Getting Help
~~~~~~~~~~~

For additional support:

1. **Documentation**: Refer to :doc:`../api/index` for detailed API documentation
2. **Performance**: See :doc:`../performance/index` for optimization guidance
3. **Community**: Check project issues and discussions
4. **Debugging**: Use :doc:`../contributing` for development tools

Contributing Examples
~~~~~~~~~~~~~~~~~~~~

To contribute new examples:

1. Follow the existing code structure and documentation style
2. Include comprehensive comments explaining key concepts
3. Provide expected output and performance characteristics
4. Test on multiple platforms and configurations
5. Update this index with appropriate links and descriptions