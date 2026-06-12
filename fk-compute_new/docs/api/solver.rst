``fkcompute.solver`` — ILP Generation
======================================

The ``solver`` subpackage implements Phase 2 of the computation pipeline:
translating the reduced constraint system into a linear programming problem
and writing the ILP ``.csv`` file consumed by the C++ backend.

ILP Generation
--------------

.. automodule:: fkcompute.solver.ilp
   :members:
   :undoc-members:
   :show-inheritance:

Usage Example
~~~~~~~~~~~~~

.. code-block:: python

   from fkcompute.domain.braid.states import BraidStates
   from fkcompute.domain.constraints.reduction import full_reduce
   from fkcompute.solver.ilp import ilp

   bs = BraidStates([1, 1, 1])
   relations = bs.get_state_relations()
   reduced = full_reduce(relations)

   ilp_path = ilp(degree=2, relations=reduced, braid_states=bs,
                  outfile="trefoil_ilp.csv")
   print("ILP written to:", ilp_path)

Gurobi Integration
------------------

Gurobi (``gurobipy``) is used in two places:

1. :func:`~fkcompute.solver.ilp.integral_bounded` — feasibility check.
2. :func:`~fkcompute.solver.ilp._check_sign_assignment` — validates a
   candidate sign assignment during the inversion search.

Gurobi must be installed and licensed separately.  See :doc:`../installation`
for details.

ILP File Format
---------------

The ``.csv`` file is a custom format (not standard LP or MPS).  It encodes:

* The degree bound and variable count.
* The integer domain of each free variable (lower and upper bounds).
* Inequality and conservation constraints as linear coefficient rows.
* R-matrix type labels for each crossing.

This format is consumed exclusively by the ``fk_main`` binary and is
considered an internal implementation detail.  Do not rely on it being
stable across package versions.
