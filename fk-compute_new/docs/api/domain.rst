``fkcompute.domain`` — Braid and Constraint Domain
====================================================

The domain subpackage contains all pure mathematical logic with no I/O
or side effects.  It is organised into two main areas:

* **Braid** — data structures for braid words, topological data, and sign
  assignments.
* **Constraints** — classes and reduction algorithms for the constraint system
  derived from R-matrix equations.

``fkcompute.domain.braid``
--------------------------

Braid Topology
~~~~~~~~~~~~~~

.. autoclass:: fkcompute.domain.braid.topology.BraidTopology
   :members:
   :undoc-members:
   :show-inheritance:

   **Attributes set by** ``__init__``

   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Attribute
        - Description
      * - ``braid``
        - The braid word as ``list[int]``.
      * - ``max_strand``
        - Largest generator index present.
      * - ``n_strands``
        - Total number of strands (``max_strand + 1``).
      * - ``n_crossings``
        - Number of crossings (``len(braid)``).
      * - ``crossing_signs``
        - ``list[int]`` of +1 / −1 per crossing.
      * - ``writhe``
        - Sum of crossing signs.
      * - ``n_components``
        - Number of components in the braid closure.
      * - ``states``
        - Canonical state representatives (strand, column) tuples.
      * - ``state_equivalence_classes``
        - Groups of locations that map to the same state.
      * - ``strand_endpoints``
        - Dict mapping component index → list of endpoint location pairs.

Signed Braid
~~~~~~~~~~~~

.. autoclass:: fkcompute.domain.braid.signed.SignedBraid
   :members:
   :undoc-members:
   :show-inheritance:

Braid States (High-Level Wrapper)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: fkcompute.domain.braid.states.BraidStates
   :members:
   :undoc-members:
   :show-inheritance:

   ``BraidStates`` is the primary entry point for braid computations.  It
   combines :class:`~fkcompute.domain.braid.topology.BraidTopology` and
   :class:`~fkcompute.domain.braid.signed.SignedBraid` and provides the
   :meth:`get_state_relations` method used by the pipeline.

   Example:

   .. code-block:: python

      from fkcompute.domain.braid.states import BraidStates

      bs = BraidStates([1, 1, 1])
      print(bs.n_components)  # 1
      print(bs.writhe)        # 3

      relations = bs.get_state_relations()
      print(f"{len(relations)} constraints generated")

Braid Word Utilities
~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.domain.braid.word
   :members:
   :undoc-members:

State Types
~~~~~~~~~~~

.. automodule:: fkcompute.domain.braid.types
   :members:
   :undoc-members:


``fkcompute.domain.constraints``
---------------------------------

Constraint Relation Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.domain.constraints.relations
   :members:
   :undoc-members:
   :show-inheritance:

   The constraint classes form a hierarchy:

   .. code-block:: text

      Constraint (ABC)
      ├── Leq          — first ≤ second
      ├── Less         — first < second
      ├── Zero         — state = [0]
      ├── NegOne       — state = [-1]
      ├── Alias        — state ≡ alias
      └── Conservation — Σ inputs = Σ outputs

   **Leq** and **Less** encode ordering between states.

   **Zero** and **NegOne** pin a state to the special values
   :data:`~fkcompute.ZERO_STATE` and :data:`~fkcompute.NEG_ONE_STATE`.

   **Alias** identifies two state variables as equal; the constraint reducer
   uses this to eliminate variables.

   **Conservation** encodes the R-matrix conservation law: the sum of
   incoming state values equals the sum of outgoing values at each crossing.

Constraint Reduction
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.domain.constraints.reduction
   :members:
   :undoc-members:

   ``full_reduce`` runs a fixed-point loop applying:

   1. Alias propagation — substitutes alias pairs throughout all other
      constraints.
   2. Constant propagation — resolves Zero/NegOne constraints by substituting
      their values everywhere.
   3. Variable elimination — removes constraints that become trivially
      satisfied after substitution.

   The loop terminates when no more reductions can be applied.

Symbolic Constraint System
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.domain.constraints.symbols
   :members:
   :undoc-members:

.. automodule:: fkcompute.domain.constraints.system
   :members:
   :undoc-members:
   :no-index:

.. automodule:: fkcompute.domain.constraints.pipeline
   :members:
   :undoc-members:


``fkcompute.domain.solver``
-----------------------------

.. automodule:: fkcompute.domain.solver.assignment
   :members:
   :undoc-members:

.. automodule:: fkcompute.domain.solver.symbolic_constraints
   :members:
   :undoc-members:
