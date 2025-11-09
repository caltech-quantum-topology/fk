FK Computation C++ Library Documentation
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   getting_started
   user_guide
   examples
   build_instructions

Welcome to the FK Computation C++ Library documentation. This library provides high-performance computation of FK invariants for braids using advanced mathematical algorithms and optimized C++ implementations.

Overview
--------

The FK Computation C++ Library is the core computational engine for calculating FK invariants of braids. It features:

* **High-Performance Polynomial Operations**: Multiple polynomial implementations optimized for different use cases
* **Advanced Mathematical Algorithms**: Sophisticated inversion and ILP reduction techniques
* **Flexible Architecture**: Pluggable polynomial backends and configurable computation parameters
* **Parallel Processing**: Multi-threaded computation support for improved performance
* **Comprehensive API**: Well-designed C++ interfaces for integration into larger systems

Key Features
------------

Core Mathematical Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Polynomial Systems**: Multiple polynomial implementations including MultivariablePolynomial, FMPoly (FLINT-based), and BMPoly
* **Laurent Polynomials**: Bilvector implementation supporting negative exponents
* **Linear Algebra**: Optimized linear algebra operations for matrix computations
* **Inequality Solving**: Integer solution finding for systems of linear inequalities

Computation Engine
^^^^^^^^^^^^^^^^^^

* **FK Computation**: Main orchestrator for FK invariant calculations
* **Configuration Management**: Flexible parameter configuration and validation
* **Result Processing**: Comprehensive output formatting and serialization
* **Error Handling**: Robust error management and validation

Performance Optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^

* **FLINT Integration**: High-performance arithmetic using the FLINT library
* **Memory Efficiency**: Optimized data structures for large-scale computations
* **Parallel Processing**: OpenMP-based parallelization support
* **Template-Based Design**: Compile-time optimizations through C++ templates

Architecture
------------

The library is organized into several key modules:

**Core Computation** (``fk/``)
   Main FK computation algorithms and orchestration

**Polynomial Systems** (``fk/polynomial_*.hpp``)
   Various polynomial implementations with unified interfaces

**Mathematical Utilities** (``fk/linalg.hpp``, ``fk/qalg_links.hpp``)
   Linear algebra and mathematical helper functions

**Data Structures** (``fk/bilvector.hpp``, ``fk/btree.hpp``)
   Specialized data structures for efficient computation

**Solver Components** (``fk/inequality_solver.hpp``)
   Integer programming and constraint solving

Quick Start
-----------

Here's a simple example of using the library:

.. code-block:: cpp

   #include "fk/fk_computation.hpp"
   #include "fk/multivariable_polynomial.hpp"

   int main() {
       // Create FK computation instance
       fk::FKComputation computation;

       // Run computation from input file
       computation.compute("input_data", "output_results");

       // Get the result
       const auto& result = computation.getLastResult();

       return 0;
   }

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`