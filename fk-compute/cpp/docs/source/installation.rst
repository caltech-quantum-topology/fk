Installation Guide
==================

This guide covers the installation and setup of the FK Computation C++ Library.

System Requirements
------------------

**Operating Systems**
   - Linux (Ubuntu 18.04+, CentOS 7+, Arch Linux)
   - macOS (10.14+)
   - Windows (with WSL2 or MinGW-w64)

**Compiler Requirements**
   - GCC 7.0+ or Clang 5.0+ with C++17 support
   - MSVC 2017+ for Windows native builds

**Dependencies**
   - Make build system
   - pthread library (for parallel processing)
   - Optional: clang-format for code formatting

Quick Installation
------------------

**Ubuntu/Debian**

.. code-block:: bash

   # Install build tools
   sudo apt update
   sudo apt install build-essential cmake git

   # Clone repository
   git clone <repository-url>
   cd fk-compute/cpp

   # Build library
   make all

**CentOS/RHEL**

.. code-block:: bash

   # Install build tools
   sudo yum groupinstall "Development Tools"
   sudo yum install cmake git

   # Clone and build
   git clone <repository-url>
   cd fk-compute/cpp
   make all

**macOS**

.. code-block:: bash

   # Install Xcode command line tools
   xcode-select --install

   # Or use Homebrew
   brew install gcc cmake git

   # Clone and build
   git clone <repository-url>
   cd fk-compute/cpp
   make all

Build Targets
-------------

The library provides several build targets for different use cases:

**Core Targets**

.. code-block:: bash

   make all                    # Build main executable and examples
   make fk_segments_links      # Main FK computation executable
   make polynomial_example     # Polynomial operations demo
   make parallel_demo          # Parallel processing benchmark

**Testing Targets**

.. code-block:: bash

   make run-test              # Run comprehensive test suite
   make run-arithmetic        # Test polynomial arithmetic
   make run-parallel-demo     # Benchmark parallel performance

**Development Targets**

.. code-block:: bash

   make debug                 # Build with debug symbols
   make check                 # Verify compilation
   make format                # Format code (requires clang-format)
   make clean                 # Remove build artifacts

Configuration Options
--------------------

**Compiler Flags**

The default build uses optimization flags suitable for production:

.. code-block:: makefile

   CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude

For debug builds:

.. code-block:: bash

   make debug  # Adds -g -DDEBUG flags

**Thread Configuration**

The library automatically detects available CPU cores. To override:

.. code-block:: cpp

   // Specify thread count explicitly
   parallelPooling(inequalities, supporting, callback, 8);

   // Use all available cores (default)
   parallelPooling(inequalities, supporting, callback, 0);

**Memory Configuration**

For large problems, you may need to adjust system limits:

.. code-block:: bash

   # Increase virtual memory limit
   ulimit -v unlimited

   # Increase stack size for deep recursion
   ulimit -s 16384

Installation Verification
------------------------

Verify your installation by running the test suite:

.. code-block:: bash

   # Quick verification
   make run-test

   # Performance verification
   make run-parallel-demo --benchmark

Expected output:

.. code-block:: text

   === FK Computation Test Suite ===
   ✓ Polynomial arithmetic tests passed
   ✓ Solution pool algorithms passed
   ✓ Parallel processing tests passed

   === Parallel Performance Test ===
   Hardware threads detected: 24
   Problem size: medium
   Sequential: 1,329 solutions in 2ms
   Parallel (24 threads): 1,329 solutions in 3ms
   ✓ Solution counts match!

Troubleshooting
---------------

**Common Build Issues**

*Compiler not found*:

.. code-block:: bash

   # Specify compiler explicitly
   make CXX=g++-9 all

*Missing C++17 support*:

.. code-block:: bash

   # Check compiler version
   g++ --version
   clang++ --version

   # Update if needed (Ubuntu)
   sudo apt install gcc-9 g++-9

*pthread linking errors*:

.. code-block:: bash

   # Ensure pthread is available
   sudo apt install libc6-dev

**Runtime Issues**

*Permission denied*:

.. code-block:: bash

   # Make executables runnable
   chmod +x fk_segments_links parallel_demo

*Segmentation fault*:

.. code-block:: bash

   # Run with debug build
   make debug
   gdb ./fk_segments_links

*Memory allocation failures*:

.. code-block:: bash

   # Check available memory
   free -h

   # Monitor memory usage
   valgrind --tool=massif ./parallel_demo

**Performance Issues**

*Poor parallel performance*:

- Ensure problem size is large enough (>1000 solutions)
- Check CPU core count: ``nproc``
- Monitor CPU usage: ``htop``
- Verify thread count in output

*High memory usage*:

- Use sparse polynomial representation (default)
- Monitor with: ``ps aux | grep fk_segments_links``

Development Installation
-----------------------

For development work, install additional tools:

.. code-block:: bash

   # Code formatting
   sudo apt install clang-format

   # Memory debugging
   sudo apt install valgrind

   # Performance profiling
   sudo apt install linux-tools-generic

   # Documentation building
   pip install sphinx sphinx-rtd-theme

Environment Setup
----------------

Add these to your ``.bashrc`` or ``.zshrc`` for convenience:

.. code-block:: bash

   # FK Computation aliases
   alias fk-build='cd /path/to/fk-compute/cpp && make all'
   alias fk-test='cd /path/to/fk-compute/cpp && make run-test'
   alias fk-bench='cd /path/to/fk-compute/cpp && make run-parallel-demo'

   # Development environment
   export FK_COMPUTE_HOME=/path/to/fk-compute/cpp
   export PATH=$FK_COMPUTE_HOME:$PATH

Next Steps
----------

After installation, proceed to:

- :doc:`quickstart` - Basic usage examples
- :doc:`tutorials/index` - Step-by-step tutorials
- :doc:`api/polynomials` - API documentation