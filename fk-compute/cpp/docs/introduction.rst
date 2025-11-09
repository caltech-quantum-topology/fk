Introduction
============

The FK Computation C++ Library is a high-performance mathematical computing library designed for calculating FK invariants of braids. It serves as the computational core for the larger FK computation system, providing optimized algorithms and data structures for complex mathematical operations.

What are FK Invariants?
------------------------

FK invariants are mathematical objects associated with braids that capture important topological and algebraic properties. The computation of these invariants involves sophisticated mathematical techniques including:

* **Polynomial Arithmetic**: Operations on multivariable polynomials with support for negative exponents
* **Linear Programming**: Integer linear programming (ILP) for optimization problems
* **Topological Analysis**: Analysis of braid structure and crossing relationships
* **Algebraic Computation**: Advanced algebraic operations on polynomial rings

Mathematical Background
-----------------------

The library implements algorithms for computing FK invariants through several key mathematical concepts:

Polynomial Representations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The library supports multiple polynomial representations:

1. **Multivariable Polynomials**: Sparse representation using hash maps for coefficients
2. **FLINT-based Polynomials**: Dense representation using the FLINT library for performance
3. **Laurent Polynomials**: Support for negative exponents using the bilvector data structure

The general form of polynomials handled is:

.. math::

   P(q, x_1, x_2, \ldots, x_n) = \sum_{i,j_1,j_2,\ldots,j_n} c_{i,j_1,j_2,\ldots,j_n} \cdot q^i \cdot x_1^{j_1} \cdot x_2^{j_2} \cdots x_n^{j_n}

where :math:`i` can be negative (Laurent polynomials) and coefficients are integers.

Braid Theory Foundations
^^^^^^^^^^^^^^^^^^^^^^^^

The computation works with braids represented as sequences of crossings. Key concepts include:

* **Braid Words**: Sequences of generators representing braid crossings
* **Writhe**: The sum of crossing signs in a braid diagram
* **Components**: Connected components in the braid closure
* **Crossing Relations**: Mathematical relationships between different crossings

Algorithmic Approach
^^^^^^^^^^^^^^^^^^^^

The FK computation algorithm involves several stages:

1. **Preprocessing**: Analysis of braid structure and generation of constraint systems
2. **Variable Assignment**: Systematic enumeration of variable assignments satisfying constraints
3. **Polynomial Accumulation**: Building the result polynomial through term-by-term computation
4. **Optimization**: Application of mathematical optimizations to reduce computation time

Performance Characteristics
---------------------------

The library is designed for high-performance computation with several optimization strategies:

**Memory Efficiency**
   * Sparse representations for polynomials with few terms
   * Lazy allocation of data structures
   * Memory pooling for frequent allocations

**Computational Efficiency**
   * Template-based implementations for compile-time optimization
   * FLINT library integration for high-performance arithmetic
   * Parallel processing support using OpenMP

**Scalability**
   * Algorithms designed to handle large-degree polynomials
   * Configurable memory usage limits
   * Incremental computation capabilities

Library Philosophy
------------------

The FK Computation C++ Library is built on several design principles:

**Modularity**
   Clear separation of concerns with well-defined interfaces between components

**Flexibility**
   Pluggable backends allow switching between different polynomial implementations

**Performance**
   Optimized for computational efficiency without sacrificing code clarity

**Robustness**
   Comprehensive error handling and input validation

**Extensibility**
   Template-based design facilitates extension to new mathematical structures

Integration with Python
------------------------

While this is a C++ library, it's designed to integrate seamlessly with the Python FK computation package:

* **Compiled Extensions**: C++ code is compiled into Python extensions for maximum performance
* **Data Interchange**: Efficient serialization/deserialization between C++ and Python
* **Configuration**: Python-configurable parameters passed to C++ computation core
* **Result Processing**: C++ results can be processed by Python symbolic mathematics libraries

Next Steps
----------

To get started with the FK Computation C++ Library:

1. Read the :doc:`getting_started` guide for installation and basic usage
2. Explore the :doc:`user_guide` for detailed usage patterns
3. Browse the :doc:`api/library_root` for complete API reference
4. Check out the :doc:`examples` for practical code examples