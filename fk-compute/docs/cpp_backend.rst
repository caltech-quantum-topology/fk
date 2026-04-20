C++ Backend
===========

The C++ backend, compiled to the ``fk_main`` binary, performs the
computationally intensive state-sum evaluation of the FK polynomial.  It is
invoked automatically during Phase 3 of the computation pipeline.

Architecture
------------

The backend is organized around four main concerns:

1. **Input parsing** — reads the ILP ``.csv`` file produced by the Python layer.
2. **State enumeration** — iterates over all valid state assignments compatible
   with the constraint system.
3. **Polynomial arithmetic** — accumulates the FK polynomial using
   arbitrary-precision coefficients.
4. **Output serialization** — writes the sparse polynomial to a JSON file.

Source Layout
-------------

.. code-block:: text

   cpp/
   ├── main.cpp                    # Binary entry point
   ├── src/
   │   ├── fk_computation.cpp      # Core state-sum driver
   │   ├── fmpoly.cpp              # FLINT modular polynomial wrapper
   │   ├── fmpoly_class.cpp        # FLINT polynomial class
   │   ├── qalg_links.cpp          # Link-specific computations
   │   ├── linalg.cpp              # Linear algebra (BLAS/OpenMP)
   │   └── inequality_solver.cpp   # Constraint checking
   └── include/
       └── *.h                     # Corresponding header files

Dependencies
------------

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Library
     - Role
     - Notes
   * - **FLINT**
     - Arbitrary-precision polynomial arithmetic
     - Also brings in **GMP** as a dependency
   * - **OpenMP**
     - Shared-memory parallelism
     - Controlled via ``--threads N``
   * - **BLAS** (OpenBLAS)
     - Dense linear algebra
     - Used in ``linalg.cpp``

Polynomial Arithmetic Strategy
--------------------------------

FK polynomial coefficients can grow very large for high-degree computations
(thousands of decimal digits).  The backend uses FLINT to manage large integer
computations.

Parallelism
-----------

The outer loop over state assignments is parallelised with OpenMP:

.. code-block:: cpp

   #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
   for (auto& assignment : all_assignments) {
       // evaluate R-matrices and accumulate contribution
   }

Each thread accumulates its own partial polynomial; the master thread
reduces them with a lock-protected sum at the end.

Set ``threads`` in your configuration file or ``--threads N`` on the CLI
to control the thread count.

Building the Backend
--------------------

The preferred build system is **CMake** via scikit-build-core (invoked
automatically by ``pip install``):

.. code-block:: bash

   # Build during pip install (automatic)
   pip install .

   # Manual CMake build for development
   cmake -B build -S .
   cmake --build build

   # Copy the binary into the Python package location
   cp build/fk_main src/fkcompute/_bin/fk_main

A legacy ``Makefile`` in ``cpp/`` is retained for reference but CMake is
preferred for reproducibility.

Binary Resolution
-----------------

When ``run_fk_binary`` is called, the binary is located in this order:

1. An ``fk_main`` binary on the system ``PATH``.
2. The bundled binary at ``fkcompute/_bin/fk_main`` (packaged inside the
   wheel at install time).

You can override the bundled binary by placing your own ``fk_main`` earlier
on ``PATH``.

Input / Output Format
----------------------

Input (``.csv``)
~~~~~~~~~~~~~~~~

The ``.csv`` file is a custom constraint format produced by
:func:`fkcompute.solver.ilp.ilp`.  It encodes:

* The degree bound.
* The number and ranges of the free integer variables.
* The constraint matrix (inequalities and conservation laws).
* R-matrix type assignments for each crossing.

Output (``.json``)
~~~~~~~~~~~~~~~~~~

The backend writes a JSON file with this structure:

.. code-block:: json

   {
     "terms": [
       {"x": [0], "q_terms": [{"q": -2, "c": "1"}, {"q": 0, "c": "-1"}]},
       {"x": [2], "q_terms": [{"q": 4, "c": "1"}]}
     ],
     "metadata": {
       "num_x_variables": 1,
       "max_x_degrees": [2],
       "storage_type": "flint"
     }
   }

See :doc:`output` for the full description of this format.
