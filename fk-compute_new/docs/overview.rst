Mathematical Overview
=====================

This page provides background on the FK invariant and explains at a
high level how **fkcompute** computes it.  Knowledge of knot theory at
the level of braid groups and quantum invariants is helpful but not
required for using the package.

The FK Invariant
-----------------

The **FK invariant** :math:`F_K(x, q)` is a two-variable power series
introduced by Gukov and Manolescu [GukovManolescu2021]_.
It is an invariant of knot and link complements: two equivalent knots
have the same :math:`F_K`.

For a knot :math:`K`, the invariant takes the form

.. math::

   F_K(x, q) = \sum_{n, k} a_{n,k}\, x^n\, q^k \;\in\; \mathbb{Z}[\![x, q^{\pm 1}]\!],

where :math:`n` ranges over non-negative integers and :math:`k` ranges over
integers. 

For a link :math:`L` with :math:`\ell` components, each component contributes
its own :math:`x` variable, so the invariant lives in

.. math::

   F_L(x_1, \ldots, x_\ell, q) \;\in\; \mathbb{Z}[\![x_1, \ldots, x_\ell, q^{\pm 1}]\!].

Relationship to Other Invariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`F_K` extends several classical quantum invariants:

* :math:`F_K` agrees with the WRT invariant at certain roots of unity.
* :math:`F_K` can be regarded as an extension of the :math:`\widehat{Z}` invariant to manifolds with boundary.

State-Sum / R-Matrix Computation
----------------------------------

**fkcompute** implements the **large-colour R-matrix / inverted state-sum**
approach developed by Park [Park2020a]_ [Park2020b]_.

Intuitively, the braid closure is cut into elementary crossings.  At each
crossing one of two *R-matrices* and hteir inverses (:math:`R_1, R_2=R_1^{-1}, R_3, R_4=R_3^{-1}`) acts,
depending on the crossing sign and the *sign assignment* (see below).
The FK polynomial is then assembled as a weighted sum over all *states*
— assignments of integers to the strands at each intermediate slice of
the braid.

Formally, the state sum reads

.. math::

   F_K(x,q) = \sum_{\text{states } s} w(s)\, \prod_{\text{crossings}} R_{\text{type}}(s),

where the weight :math:`w(s)` encodes the x-variables and the
degree truncation cuts the sum to a finite computation.

Sign Assignment (Inversion Phase)
-----------------------------------

A crucial ingredient in the Park approach is a **sign assignment**: a
choice of :math:`\pm 1` for each strand segment of the braid.  Not
every assignment is consistent with the R-matrix equations; finding a
valid one is referred to internally as *inversion*.

For **homogeneous braids** (every generator appears with a single sign),
the sign assignment is determined uniquely and computed automatically.

For general braids, fkcompute searches for a valid sign assignment by:

1. Generating *braid variants* — cyclic shifts and reflections of the
   braid word.
2. For each variant, enumerating sign-assignment *permutation candidates*
   (derived from a combinatorial structure on the sign diagrams).
3. Testing each candidate against the constraint system using Gurobi.

See :doc:`pipeline` for implementation details.

Integer Linear Programming
---------------------------

Once a sign assignment is found, the state-sum constraints are cast as
an **integer linear program (ILP)**.  The constraint variables are the
integer values assigned to each braid state, and the constraints encode:

* R-matrix consistency (equality, inequality, and conservation laws).
* Degree bounds derived from the truncation parameter.

The ILP feasibility check (using Gurobi) validates the constraint system,
and the reduced constraint matrix is written to a ``.csv`` file consumed
by the C++ backend.

C++ Backend
------------

The heavy computation — evaluating the state sum at each term — is
performed by the compiled ``fk_main`` binary.  The backend uses:

* **FLINT** for arbitrary-precision polynomial arithmetic (coefficients
  can grow very large for high-degree computations).
* **OpenMP** for shared-memory parallelism across state assignments.
* **BLAS/OpenBLAS** for dense linear algebra sub-routines.

The binary reads the ILP ``.csv`` file and writes a JSON result containing
sparse polynomial terms.

Degree Truncation
------------------

Computing :math:`F_K` exactly requires an infinite state sum.  The
``degree`` parameter limits the maximum power of :math:`x` (and
equivalently of :math:`x_1, \ldots, x_\ell` for links) that appears in
the output.  Higher degrees are more expensive but recover more terms.

References
-----------

.. [GukovManolescu2021] S. Gukov and C. Manolescu,
   "A two-variable series for knot complements,"
   *Quantum Topology* **12** (2021), no. 1, 1–109.
   `arXiv:1904.06057 <https://arxiv.org/abs/1904.06057>`_.

.. [Park2020a] S. Park,
   "Large color R-matrix for knot complements and strange identities,"
   *Journal of Knot Theory and Its Ramifications* **29** (2020).
   `arXiv:2004.02087 <https://arxiv.org/abs/2004.02087>`_.

.. [Park2020b] S. Park,
   "Inverted state sums, inverted Habiro series, and indefinite theta functions,"
   2021. `arXiv:2106.03495 <https://arxiv.org/abs/2106.03495>`_.
