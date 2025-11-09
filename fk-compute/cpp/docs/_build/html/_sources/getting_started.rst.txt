Getting Started
===============

This guide will help you set up and start using the FK Computation C++ Library.

Prerequisites
-------------

Before building the FK Computation C++ Library, ensure you have the following dependencies installed:

Required Dependencies
^^^^^^^^^^^^^^^^^^^^^

**C++ Compiler**
   * GCC 8.0+ or Clang 10.0+ with C++17 support
   * MSVC 2019+ on Windows

**FLINT Library**
   * FLINT 2.9.0 or later for high-performance arithmetic
   * Installation instructions: https://flintlib.org/

**OpenMP** (Optional)
   * For parallel processing support
   * Usually included with modern GCC/Clang installations

**CMake** (Recommended)
   * CMake 3.15+ for building
   * Alternative: direct compilation with Makefiles

Installation on Ubuntu/Debian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Install build tools
   sudo apt update
   sudo apt install build-essential cmake

   # Install FLINT
   sudo apt install libflint-dev

   # Install OpenMP (usually included with GCC)
   sudo apt install libomp-dev

Installation on macOS
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Using Homebrew
   brew install cmake
   brew install flint
   brew install libomp

Installation on Windows
^^^^^^^^^^^^^^^^^^^^^^^

1. Install Visual Studio 2019 or later with C++ support
2. Install vcpkg package manager
3. Install dependencies:

.. code-block:: cmd

   vcpkg install flint:x64-windows
   vcpkg install openmp:x64-windows

Building the Library
--------------------

Using CMake (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Clone or navigate to the project directory:

.. code-block:: bash

   cd /path/to/fk-compute/cpp

2. Create a build directory:

.. code-block:: bash

   mkdir build
   cd build

3. Configure the build:

.. code-block:: bash

   cmake ..

4. Build the library:

.. code-block:: bash

   make -j$(nproc)

Using Direct Compilation
^^^^^^^^^^^^^^^^^^^^^^^^

For simple builds without CMake:

.. code-block:: bash

   # Basic compilation
   g++ -std=c++17 -O3 -I./include -lflint -fopenmp \
       src/*.cpp examples/example.cpp -o fk_example

Build Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The library supports several configuration options:

**Polynomial Backend Selection**

Choose the polynomial implementation at compile time:

.. code-block:: bash

   # Use MultivariablePolynomial (default)
   cmake -DPOLYNOMIAL_TYPE=0 ..

   # Use FMPoly (FLINT-based)
   cmake -DPOLYNOMIAL_TYPE=1 ..

   # Use BMPoly (Basic implementation)
   cmake -DPOLYNOMIAL_TYPE=2 ..

**Debug vs Release**

.. code-block:: bash

   # Debug build
   cmake -DCMAKE_BUILD_TYPE=Debug ..

   # Release build (default)
   cmake -DCMAKE_BUILD_TYPE=Release ..

**OpenMP Support**

.. code-block:: bash

   # Enable OpenMP
   cmake -DUSE_OPENMP=ON ..

   # Disable OpenMP
   cmake -DUSE_OPENMP=OFF ..

First Example
-------------

Let's create a simple program to test the installation:

**hello_fk.cpp**

.. code-block:: cpp

   #include <iostream>
   #include "fk/multivariable_polynomial.hpp"

   int main() {
       std::cout << "FK Computation Library - Hello World!" << std::endl;

       // Create a simple polynomial
       MultivariablePolynomial poly(2);  // 2 variables

       // Set coefficient: 3*q^1*x₁*x₂
       poly.setCoefficient(1, {1, 1}, 3);

       // Set coefficient: 2*q^0 (constant term)
       poly.setCoefficient(0, {0, 0}, 2);

       std::cout << "Created polynomial:" << std::endl;
       poly.print();

       return 0;
   }

Compile and run:

.. code-block:: bash

   g++ -std=c++17 -I./include -lflint hello_fk.cpp -o hello_fk
   ./hello_fk

Expected output:

.. code-block:: text

   FK Computation Library - Hello World!
   Created polynomial:
   2*q^0 + 3*q^1*x_1*x_2

Running Examples
----------------

The library comes with several example programs demonstrating different features:

**Basic Polynomial Operations**

.. code-block:: bash

   cd examples
   g++ -std=c++17 -I../include -lflint example.cpp -o example
   ./example

**Bilvector (Laurent Polynomial) Operations**

.. code-block:: bash

   g++ -std=c++17 -I../include bilvector_operations.cpp -o bilvector_example
   ./bilvector_example

**FK Computation Example**

.. code-block:: bash

   g++ -std=c++17 -I../include -lflint -fopenmp \
       fk_computation_test.cpp -o fk_test
   ./fk_test

**Parallel Processing Demo**

.. code-block:: bash

   g++ -std=c++17 -I../include -lflint -fopenmp \
       parallel_demo.cpp -o parallel_demo
   ./parallel_demo

Troubleshooting
---------------

Common Issues and Solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**FLINT Not Found**

If you get errors about missing FLINT headers:

.. code-block:: bash

   # Make sure FLINT is installed
   pkg-config --cflags --libs flint

   # If FLINT is in a custom location
   export PKG_CONFIG_PATH=/path/to/flint/lib/pkgconfig:$PKG_CONFIG_PATH

**OpenMP Linking Issues**

If you get OpenMP linking errors:

.. code-block:: bash

   # For GCC
   g++ -fopenmp ...

   # For Clang on macOS
   g++ -Xpreprocessor -fopenmp -lomp ...

**Template Compilation Errors**

If you get template-related compilation errors:

1. Ensure you're using C++17 or later
2. Check that all header files are properly included
3. Verify template specializations match your usage

**Performance Issues**

If the library seems slow:

1. Make sure you're using a Release build (``-O3`` optimization)
2. Enable OpenMP for parallel processing
3. Consider using FMPoly backend for better performance with dense polynomials

Getting Help
^^^^^^^^^^^^

If you encounter issues:

1. Check the :doc:`api/library_root` for detailed API documentation
2. Look at the :doc:`examples` for usage patterns
3. Review the source code in the ``include/fk/`` directory
4. Check compiler warnings and error messages carefully

Next Steps
----------

Now that you have the library installed and working:

1. Read the :doc:`user_guide` for detailed usage information
2. Explore the :doc:`examples` to see practical applications
3. Browse the :doc:`api/library_root` for complete API reference
4. Try modifying the examples to experiment with different features