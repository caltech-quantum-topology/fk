Build Instructions
==================

This document provides comprehensive instructions for building the FK Computation C++ Library and its documentation.

Prerequisites
-------------

System Requirements
^^^^^^^^^^^^^^^^^^^

**Supported Platforms:**
- Linux (Ubuntu 18.04+, CentOS 7+, Arch Linux, etc.)
- macOS (10.14+)
- Windows (Windows 10 with WSL2 or native MSVC)

**Required Tools:**
- C++ compiler with C++17 support (GCC 8+, Clang 10+, MSVC 2019+)
- CMake 3.15 or later (recommended for building)
- Make or Ninja build system

**Required Libraries:**
- FLINT 2.9.0+ (Fast Library for Number Theory)
- OpenMP (optional, for parallel processing)

**Documentation Requirements:**
- Python 3.7+ (for Sphinx documentation)
- Doxygen 1.8.0+ (for C++ API documentation)
- Sphinx and extensions (see documentation setup below)

Installing Dependencies
-----------------------

Ubuntu/Debian
^^^^^^^^^^^^^

.. code-block:: bash

   # Update package list
   sudo apt update

   # Install build tools
   sudo apt install build-essential cmake ninja-build

   # Install FLINT library
   sudo apt install libflint-dev

   # Install OpenMP (usually included with GCC)
   sudo apt install libomp-dev

   # Install documentation tools
   sudo apt install doxygen python3-pip
   pip3 install sphinx breathe exhale sphinx-rtd-theme

CentOS/RHEL/Fedora
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # For CentOS/RHEL
   sudo yum groupinstall "Development Tools"
   sudo yum install cmake ninja-build

   # For Fedora
   sudo dnf groupinstall "Development Tools"
   sudo dnf install cmake ninja-build

   # Install FLINT (may need EPEL repository)
   sudo yum install epel-release  # CentOS/RHEL only
   sudo yum install flint-devel   # or dnf for Fedora

   # Install documentation tools
   sudo yum install doxygen python3-pip  # or dnf for Fedora
   pip3 install sphinx breathe exhale sphinx-rtd-theme

macOS
^^^^^

.. code-block:: bash

   # Install Homebrew if not already installed
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

   # Install build tools
   brew install cmake ninja

   # Install FLINT
   brew install flint

   # Install OpenMP for Clang
   brew install libomp

   # Install documentation tools
   brew install doxygen python3
   pip3 install sphinx breathe exhale sphinx-rtd-theme

Windows
^^^^^^^

**Option 1: WSL2 (Recommended)**

Install WSL2 and Ubuntu, then follow Ubuntu instructions above.

**Option 2: Native Windows with vcpkg**

.. code-block:: cmd

   # Install Visual Studio 2019 or later with C++ support
   # Install vcpkg package manager
   git clone https://github.com/Microsoft/vcpkg.git
   cd vcpkg
   .\\bootstrap-vcpkg.bat

   # Install dependencies
   .\\vcpkg install flint:x64-windows
   .\\vcpkg install openmp:x64-windows

   # Install Python and documentation tools
   # Download and install Python from python.org
   pip install sphinx breathe exhale sphinx-rtd-theme

Building the Library
--------------------

Using CMake (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic Build:**

.. code-block:: bash

   # Clone/navigate to project directory
   cd /path/to/fk-compute/cpp

   # Create build directory
   mkdir build
   cd build

   # Configure build
   cmake ..

   # Build library and examples
   make -j$(nproc)

   # Or with Ninja for faster builds
   cmake -GNinja ..
   ninja

**Configuration Options:**

.. code-block:: bash

   # Choose polynomial backend
   cmake -DPOLYNOMIAL_TYPE=0 ..  # MultivariablePolynomial (default)
   cmake -DPOLYNOMIAL_TYPE=1 ..  # FMPoly (FLINT-based)
   cmake -DPOLYNOMIAL_TYPE=2 ..  # BMPoly (Basic)

   # Build type
   cmake -DCMAKE_BUILD_TYPE=Release ..    # Optimized (default)
   cmake -DCMAKE_BUILD_TYPE=Debug ..      # Debug symbols
   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..  # Optimized + debug

   # OpenMP support
   cmake -DUSE_OPENMP=ON ..   # Enable parallel processing
   cmake -DUSE_OPENMP=OFF ..  # Disable OpenMP

   # Custom install prefix
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..

   # Verbose makefiles for debugging
   cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..

**Advanced Configuration:**

.. code-block:: bash

   # Cross-compilation for different architectures
   cmake -DCMAKE_TOOLCHAIN_FILE=toolchain.cmake ..

   # Static linking
   cmake -DBUILD_SHARED_LIBS=OFF ..

   # Custom FLINT location
   cmake -DFLINT_ROOT=/custom/path/to/flint ..

   # Enable additional warnings
   cmake -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic" ..

**Windows with vcpkg:**

.. code-block:: cmd

   # Configure with vcpkg toolchain
   cmake -DCMAKE_TOOLCHAIN_FILE=C:\\path\\to\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake ..

   # Build with MSBuild
   cmake --build . --config Release

Using Direct Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^

For simple builds without CMake:

.. code-block:: bash

   # Basic compilation flags
   CXX_FLAGS="-std=c++17 -O3 -DNDEBUG"
   INCLUDE_FLAGS="-I./include"
   LINK_FLAGS="-lflint"

   # With OpenMP
   OPENMP_FLAGS="-fopenmp"

   # Compile examples
   g++ $CXX_FLAGS $INCLUDE_FLAGS $OPENMP_FLAGS examples/example.cpp $LINK_FLAGS -o example

   # Compile with specific polynomial backend
   g++ $CXX_FLAGS -DPOLYNOMIAL_TYPE=1 $INCLUDE_FLAGS examples/fk_computation_test.cpp $LINK_FLAGS $OPENMP_FLAGS -o fk_test

**Makefile Example:**

.. code-block:: makefile

   CXX = g++
   CXXFLAGS = -std=c++17 -O3 -DNDEBUG -Wall -Wextra
   INCLUDES = -I./include
   LIBS = -lflint -lgomp
   SOURCES = $(wildcard examples/*.cpp)
   TARGETS = $(SOURCES:.cpp=)

   all: $(TARGETS)

   %: %.cpp
   	$(CXX) $(CXXFLAGS) $(INCLUDES) $< $(LIBS) -o $@

   clean:
   	rm -f $(TARGETS)

   .PHONY: all clean

Building Documentation
----------------------

The documentation system uses Sphinx with Doxygen integration for comprehensive C++ API documentation.

Documentation Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Install Python documentation tools
   pip3 install -r docs/requirements.txt

   # Or install manually
   pip3 install sphinx>=4.0
   pip3 install breathe>=4.30
   pip3 install exhale>=0.3
   pip3 install sphinx-rtd-theme>=1.0

**docs/requirements.txt:**

.. code-block:: text

   sphinx>=4.0.0
   breathe>=4.30.0
   exhale>=0.3.0
   sphinx-rtd-theme>=1.0.0

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Navigate to docs directory
   cd docs

   # Clean previous builds
   make clean

   # Build HTML documentation
   make html

   # Open documentation
   open _build/html/index.html  # macOS
   xdg-open _build/html/index.html  # Linux

**Alternative build methods:**

.. code-block:: bash

   # Build with specific builder
   sphinx-build -b html . _build/html

   # Build with warnings as errors
   sphinx-build -W -b html . _build/html

   # Build PDF documentation (requires LaTeX)
   make latexpdf

   # Build EPUB documentation
   make epub

Documentation Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The documentation is configured in ``docs/conf.py``. Key settings:

.. code-block:: python

   # Project information
   project = 'FK Computation C++ Library'
   release = '1.0.0'

   # Doxygen integration
   breathe_projects = {"fk-compute": "_build/doxygen/xml/"}
   breathe_default_project = "fk-compute"

   # Exhale for automatic API generation
   exhale_args = {
       "containmentFolder": "./api",
       "rootFileName": "library_root.rst",
       "exhaleExecutesDoxygen": True,
       "exhaleDoxygenStdin": "INPUT = ../include\\n"
   }

Troubleshooting Build Issues
----------------------------

Common Problems and Solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**FLINT Not Found:**

.. code-block:: bash

   # Check FLINT installation
   pkg-config --cflags --libs flint

   # If not found, install or set PKG_CONFIG_PATH
   export PKG_CONFIG_PATH=/path/to/flint/lib/pkgconfig:$PKG_CONFIG_PATH

   # Or specify FLINT location in CMake
   cmake -DFLINT_ROOT=/path/to/flint ..

**OpenMP Issues:**

.. code-block:: bash

   # GCC: Use -fopenmp
   g++ -fopenmp ...

   # Clang on macOS: Install libomp and use specific flags
   brew install libomp
   clang++ -Xpreprocessor -fopenmp -lomp ...

   # Check OpenMP availability
   echo '#include <omp.h>' | g++ -x c++ -fopenmp -E - > /dev/null && echo "OpenMP available"

**Template Compilation Errors:**

.. code-block:: bash

   # Ensure C++17 support
   g++ -std=c++17 ...

   # Check compiler version
   g++ --version

   # Update compiler if necessary

**Memory Issues During Build:**

.. code-block:: bash

   # Reduce parallel jobs
   make -j2  # Instead of -j$(nproc)

   # Use Ninja for more efficient builds
   cmake -GNinja ..
   ninja

**Documentation Build Issues:**

.. code-block:: bash

   # Check Sphinx installation
   sphinx-build --version

   # Check Doxygen installation
   doxygen --version

   # Install missing Python packages
   pip3 install --upgrade sphinx breathe exhale

   # Clean and rebuild
   cd docs
   make clean
   make html

**Windows-Specific Issues:**

.. code-block:: cmd

   REM Ensure vcpkg integration
   vcpkg integrate install

   REM Use correct architecture
   cmake -A x64 ..

   REM Check Visual Studio version
   cmake --version

Performance Optimization
------------------------

Compiler Optimizations
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Maximum optimization
   cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" ..

   # Profile-guided optimization (advanced)
   # Step 1: Build with profiling
   cmake -DCMAKE_CXX_FLAGS="-O3 -fprofile-generate" ..
   make && ./run_benchmarks

   # Step 2: Rebuild with profile data
   cmake -DCMAKE_CXX_FLAGS="-O3 -fprofile-use" ..
   make

Link-Time Optimization
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Enable LTO for better optimization
   cmake -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE ..

   # Or manually
   cmake -DCMAKE_CXX_FLAGS="-flto" ..

Testing the Build
-----------------

Running Tests
^^^^^^^^^^^^^

.. code-block:: bash

   # Build and run basic example
   cd build/examples  # or wherever examples are built
   ./example

   # Run FK computation test
   ./fk_computation_test

   # Test parallel processing
   ./parallel_demo

   # Performance benchmark
   time ./trefoil_complete_example

Validation
^^^^^^^^^^

.. code-block:: bash

   # Check library symbols
   nm libfk_compute.a | grep -E "(FKComputation|MultivariablePolynomial)"

   # Verify OpenMP integration
   ./parallel_demo  # Should show multiple threads

   # Test with different polynomial backends
   for backend in 0 1 2; do
       echo "Testing backend $backend"
       g++ -DPOLYNOMIAL_TYPE=$backend -std=c++17 -I../include -lflint examples/example.cpp -o test_$backend
       ./test_$backend
   done

Installation
------------

System Installation
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Install library and headers
   cd build
   sudo make install

   # Or with custom prefix
   cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ..
   make install

   # Update library cache (Linux)
   sudo ldconfig

Package Creation
^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Create packages with CPack
   cd build
   cpack

   # Create specific package types
   cpack -G DEB     # Debian package
   cpack -G RPM     # RPM package
   cpack -G TGZ     # Tarball

Next Steps
----------

After successful build:

1. Run the examples to verify functionality
2. Read the :doc:`user_guide` for usage patterns
3. Explore the :doc:`api/library_root` for detailed API reference
4. Consider contributing improvements or reporting issues

For ongoing development:

1. Set up your IDE with the include paths
2. Configure debugging with appropriate symbols
3. Set up continuous integration for your project
4. Consider using package managers like Conan or vcpkg for dependencies