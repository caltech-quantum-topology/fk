FK Computation C++ Library Documentation
=========================================

.. image:: https://img.shields.io/badge/C%2B%2B-17-blue.svg
   :alt: C++17
   :target: https://isocpp.org/

.. image:: https://img.shields.io/badge/License-MIT-green.svg
   :alt: License
   :target: https://opensource.org/licenses/MIT

.. image:: https://img.shields.io/badge/Thread--Safe-Yes-brightgreen.svg
   :alt: Thread Safe

Welcome to the FK Computation C++ Library documentation! This library provides high-performance
implementations for FK computations and related algebraic topology operations with advanced
parallel processing, checkpointing, and sparse polynomial capabilities.

.. note::
   This library is designed for mathematical research and high-performance scientific computing.
   All public APIs are thread-safe and optimized for modern multi-core systems.

Quick Start
-----------

.. code-block:: bash

   # Clone and build
   git clone <repository>
   cd fk-compute/cpp
   make all

   # Run examples
   make run-parallel-demo --benchmark
   make run-checkpoint-demo

Key Features
------------

🧮 **Mathematical Components**
   - Sparse multivariable polynomials with negative exponent support
   - Iterative solution space enumeration algorithms
   - Linear algebra operations for topological computations
   - Q-algebra and quantum operations

⚡ **High Performance**
   - Multi-threaded parallel processing with work stealing
   - Automatic CPU core detection and load balancing
   - Memory-efficient sparse data structures
   - Zero-copy operations and move semantics

🛡️ **Reliability**
   - Fault-tolerant computation with binary checkpointing
   - Thread-safe concurrent data structures
   - Comprehensive error handling and recovery
   - Backward compatibility with legacy formats

🔧 **Advanced Features**
   - Bidirectional vectors supporting negative indexing
   - B-trees with substring matching capabilities
   - Atomic operations and lock-free designs
   - Progress monitoring and performance benchmarking

Architecture Overview
--------------------

The FK Computation library is built around several core components:

.. mermaid::

   graph TD
       A[FK Main Class] --> B[Multivariable Polynomials]
       A --> C[Solution Pool Algorithms]
       A --> D[Parallel Processing]
       B --> E[Sparse Storage]
       B --> F[Bilvector]
       C --> G[Constraint Solving]
       C --> H[Checkpointing]
       D --> I[Work Queue]
       D --> J[Thread Pool]

Performance Highlights
---------------------

.. class:: performance-note

   **Benchmark Results**: The parallel implementation demonstrates excellent algorithmic
   correctness with identical solution counts between sequential and parallel versions
   across all test cases.

.. table:: Performance Characteristics
   :class: performance-table

   +------------------+------------+------------+-------------------+
   | Problem Size     | Sequential | Parallel   | Solutions Found   |
   +==================+============+============+===================+
   | Small (3 vars)   | 0ms        | 0ms        | 206 ✓            |
   +------------------+------------+------------+-------------------+
   | Medium (4 vars)  | 2ms        | 3ms        | 1,329 ✓          |
   +------------------+------------+------------+-------------------+
   | Large (5 vars)   | 74ms       | 81ms       | 10,999 ✓         |
   +------------------+------------+------------+-------------------+

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   tutorials/index
   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/polynomials
   api/solution_pools
   api/parallel
   api/checkpointing
   api/data_structures
   api/utilities

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   performance/index
   threading
   memory_management
   optimization

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   testing
   benchmarking
   roadmap

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: ../README.md
   :parser: myst_parser.sphinx_