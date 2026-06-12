Development Guide
==================

This page describes how to set up a development environment, run tests, and
contribute to **fkcompute**.

Development Installation
-------------------------

Clone the repository and install in editable mode so that changes to Python
source files are immediately visible without reinstalling:

.. code-block:: bash

   git clone <repository-url>
   cd fk-compute

   # Install all optional extras in editable mode
   pip install -e ".[full]"

   # Build the C++ backend
   cmake -B build -S .
   cmake --build build

   # Copy the binary into the package location
   cp build/fk_main src/fkcompute/_bin/fk_main

After any change to C++ source files, re-run ``cmake --build build`` and copy
the binary again.

Repository Layout
-----------------

.. code-block:: text

   fk-compute/
   ├── CMakeLists.txt          # CMake build configuration
   ├── pyproject.toml          # Python package metadata (scikit-build-core)
   ├── README.md
   ├── src/
   │   └── fkcompute/          # Python source
   │       ├── __init__.py     # Public API
   │       ├── api/            # High-level entry points
   │       ├── cli/            # Typer CLI commands
   │       ├── domain/         # Pure domain logic (braid, constraints)
   │       ├── inversion/      # Phase 1: sign assignment search
   │       ├── solver/         # Phase 2: ILP generation (Gurobi)
   │       ├── output/         # Symbolic formatting (SymPy)
   │       ├── infra/          # Binary execution, config I/O
   │       ├── interactive/    # Rich-based interactive UI
   │       └── _bin/           # Bundled fk_main binary (populated at build)
   ├── cpp/                    # C++ backend source
   │   ├── main.cpp
   │   ├── src/
   │   ├── include/
   │   ├── tests/
   │   └── Makefile            # Legacy build (prefer CMake)
   ├── tests/                  # Python test suite
   │   ├── test_fk_baseline.py
   │   ├── all_knots.yaml
   │   └── all_links.yaml
   ├── mathematica/            # Wolfram Language paclet wrapper
   └── docs/                   # This documentation

Running the Test Suite
-----------------------

.. code-block:: bash

   # Run all tests
   pytest

   # Run with verbose output
   pytest -v

   # Run a specific test file
   pytest tests/test_fk_baseline.py

   # Run tests matching a name pattern
   pytest -k "trefoil"

The baseline tests compare FK invariants against a set of pre-computed
reference values for knots up to 11 crossings.

Coding Conventions
------------------

Python Style
~~~~~~~~~~~~~

* Type hints on all public functions and methods.
* NumPy-style docstrings (``Parameters``, ``Returns``, ``Examples``).
* No global mutable state in domain modules.
* No I/O in domain or constraint modules (only in ``infra/``).
* Use ``logging.getLogger("fk_logger")`` for log output; never ``print``
  except in CLI/output modules.

Module Architecture Invariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following rules maintain the layered architecture:

* ``domain/`` — pure functions and data classes only; no subprocess calls,
  no file I/O, no Gurobi imports.
* ``solver/`` — may import Gurobi; no CLI or output code.
* ``infra/`` — all file and process I/O; no domain logic.
* ``api/`` — orchestrates domain and infra; no direct CLI code.
* ``cli/`` — translates CLI arguments to API calls; minimal business logic.

C++ Style
~~~~~~~~~

* Use ``#pragma omp parallel for`` with ``schedule(dynamic)`` for the outer
  state loop.
* Allocate per-thread polynomial accumulators; reduce with a critical section.
* Use ``fmpz`` (FLINT) for all coefficient arithmetic to avoid overflow.

Building the Documentation
---------------------------

Install the documentation dependencies:

.. code-block:: bash

   pip install sphinx sphinx-rtd-theme sphinx-copybutton sphinx-design

Build the HTML documentation:

.. code-block:: bash

   cd docs
   sphinx-build -b html . _build/html

Open ``docs/_build/html/index.html`` in a browser.

Build with automatic rebuilding on changes (requires ``sphinx-autobuild``):

.. code-block:: bash

   pip install sphinx-autobuild
   sphinx-autobuild docs docs/_build/html
