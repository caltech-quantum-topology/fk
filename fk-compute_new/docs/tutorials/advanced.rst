Advanced Usage
==============

This tutorial covers less common but powerful features: supplying precomputed
inversion data, constraining the sign search, tuning performance, and using
the low-level Python API.

Supplying Precomputed Inversion Data
--------------------------------------

If you have already found a valid sign assignment (e.g. from a previous run
or from theoretical knowledge), you can skip Phase 1 entirely.

Via Python:

.. code-block:: python

   from fkcompute import fk

   # Previously computed sign assignment for the figure-eight knot
   inversion_data = {0: [1, -1, 1, -1, 1, -1, 1, -1]}

   result = fk(
       [1, -2, 1, -2],
       degree=5,
       inversion={"inversion_data": inversion_data, "braid": [1, -2, 1, -2], "degree": 5},
   )

Via config file:

.. code-block:: yaml

   braid: [1, -2, 1, -2]
   degree: 5
   inversion:
     0: [1, -1, 1, -1, 1, -1, 1, -1]

Via a saved inversion file:

.. code-block:: yaml

   braid: [1, -2, 1, -2]
   degree: 5
   inversion_file: data/figure_eight_inversion.json

Saving and reusing inversion data:

.. code-block:: python

   import json
   from fkcompute import fk
   from fkcompute.inversion.api import find_sign_assignment

   braid = [1, -2, 1, -2]
   inv_result = find_sign_assignment(braid, degree=5, max_workers=4)
   inv_dict = inv_result.to_result_dict()

   # Save for later
   with open("figure_eight_inversion.json", "w") as f:
       json.dump(inv_dict, f)

   # Reuse
   result = fk(braid, degree=8, inversion=inv_dict)


Partial Sign Constraints
-------------------------

For large braids where the full inversion search is slow, you can fix the
signs of some strand segments and let the search fill in the rest.

.. code-block:: python

   from fkcompute import fk

   # Fix the first sign of component 0, leave the rest free (0 = unknown)
   result = fk(
       [1, -2, 3, -4],
       degree=3,
       partial_signs={0: [1, 0, 0, 0], 1: [0, 0, 0, 0]},
   )

Via config file:

.. code-block:: yaml

   braid: [1, -2, 3, -4]
   degree: 3
   partial_signs:
     0: [1, 0, 0, 0]
     1: [0, 0, 0, 0]

A value of ``1`` or ``-1`` pins that strand segment; ``0`` means the
searcher may choose either sign.  This reduces the search space from
:math:`2^n` to :math:`2^{n - \text{fixed}}` candidates.


Braid Variants in the Inversion Search
----------------------------------------

By default, the inversion search tries cyclic shifts of the braid word.
You can also include the horizontally flipped braid:

.. code-block:: python

   result = fk(
       [1, -2, 1, -2],
       degree=4,
       include_flip=True,
       max_shifts=10,      # limit number of cyclic shifts
   )

``include_flip=True`` roughly doubles the search space but may find a valid
assignment faster for certain knot classes.

Low-Level API — BraidStates
-----------------------------

:class:`fkcompute.domain.braid.states.BraidStates` exposes the full
topological structure of the braid:

.. code-block:: python

   from fkcompute.domain.braid.states import BraidStates

   bs = BraidStates([1, -2, 1, -2])

   print("Strands:", bs.n_strands)             # 3
   print("Crossings:", bs.n_crossings)          # 4
   print("Components:", bs.n_components)        # 1  (figure-eight is a knot)
   print("Writhe:", bs.writhe)                  # 0
   print("Crossing signs:", bs.crossing_signs)  # [1, -1, 1, -1]
   print("States:", bs.states)

Checking Braid Properties
--------------------------

.. code-block:: python

   from fkcompute.domain.braid.word import is_positive_braid, is_homogeneous_braid

   print(is_positive_braid([1, 1, 1]))          # True
   print(is_positive_braid([1, -2, 1, -2]))     # False
   print(is_homogeneous_braid([1, 1, 1]))       # True
   print(is_homogeneous_braid([1, -2, 1, -2]))  # False

Low-Level API — Inversion
--------------------------

Run the inversion phase standalone:

.. code-block:: python

   from fkcompute.inversion.api import find_sign_assignment

   result = find_sign_assignment(
       [1, -2, 1, -2],
       degree=3,
       max_workers=4,
       chunk_size=1 << 12,
       include_flip=False,
       verbose=True,
   )

   if result.success:
       print("Braid type:", result.braid_type)
       print("Sign assignment:", result.sign_assignment)
   else:
       print("No valid sign assignment found at this degree")

Low-Level API — Constraint System
------------------------------------

Inspect the constraint system generated from a braid:

.. code-block:: python

   from fkcompute.domain.braid.states import BraidStates
   from fkcompute.domain.constraints.reduction import full_reduce

   bs = BraidStates([1, 1, 1])
   relations = bs.get_state_relations()
   reduced = full_reduce(relations)

   print(f"Relations before reduction: {len(relations)}")
   print(f"Relations after reduction:  {len(reduced)}")

   for r in reduced[:5]:
       print(" ", r)

Low-Level API — Symbolic Output
---------------------------------

Convert a result to SymPy and manipulate it:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import (
       matrix_to_polynomial,
       collect_by_variables,
       format_symbolic_output,
   )
   import sympy as sp

   result = fk([1, 1, 1], 2)
   poly = matrix_to_polynomial(result)

   # Collect by powers of x
   collected = collect_by_variables(poly)
   print(collected)

   # Extract coefficient of x^2
   x = sp.Symbol("x")
   coeff_x2 = sp.Poly(collected, x).nth(2)
   print("Coefficient of x^2:", coeff_x2)

   # Evaluate at q=1
   q = sp.Symbol("q")
   print("F_K(x, 1):", collected.subs(q, 1))

Performance Tuning
-------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Recommendation
   * - ``threads``
     - Set to the number of physical CPU cores for the C++ backend.
   * - ``max_workers``
     - Set to ``cpu_count - 1`` for the Python inversion search.
   * - ``chunk_size``
     - Default (``16384``) works well.  Reduce to ``4096`` for fast
       interactive feedback; increase to ``65536`` for deep searches.
   * - ``include_flip``
     - Disable (default) unless the standard search fails.
   * - ``max_shifts``
     - Leave ``null`` unless you want to bound search time strictly.

For high-degree computations (``degree >= 10``), the C++ backend dominates
runtime.  For non-homogeneous braids, inversion can dominate at low degrees.
Profile both phases by enabling ``verbose=True``.
