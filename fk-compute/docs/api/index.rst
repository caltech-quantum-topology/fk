API Reference
=============

This section provides auto-generated documentation for every public module
in **fkcompute**.  The documentation is extracted from the source docstrings
using Sphinx autodoc.

.. toctree::
   :maxdepth: 2

   fkcompute
   domain
   inversion
   solver
   output
   infra
   interactive

Module Summary
--------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Module
     - Description
   * - :mod:`fkcompute`
     - Top-level public API: :func:`~fkcompute.fk`, presets, and re-exports.
   * - :mod:`fkcompute.api`
     - High-level computation functions and batch processing.
   * - :mod:`fkcompute.domain`
     - Pure domain logic: braid topology, state systems, and constraints.
   * - :mod:`fkcompute.inversion`
     - Sign assignment search (Phase 1 of the pipeline).
   * - :mod:`fkcompute.solver`
     - ILP formulation and Gurobi integration (Phase 2).
   * - :mod:`fkcompute.output`
     - Symbolic polynomial formatting using SymPy.
   * - :mod:`fkcompute.infra`
     - Infrastructure: binary execution, config I/O.
   * - :mod:`fkcompute.interactive`
     - Rich-based interactive wizard (optional dependency).
