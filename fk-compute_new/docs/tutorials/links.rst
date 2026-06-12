Computing FK for Links
======================

The FK invariant generalises to links with multiple components.  For a link
:math:`L` with :math:`\ell` components, the result is a polynomial in
:math:`\ell` x-variables :math:`x_1, \ldots, x_\ell` and the
quantum parameter :math:`q`.

Recognising Multi-Component Links
-----------------------------------

A braid closure is a knot (single component) if and only if its permutation
representation is a single cycle.  The ``BraidStates`` class computes this
automatically:

.. code-block:: python

   from fkcompute.domain.braid.states import BraidStates

   # Hopf link: 2-component link
   bs = BraidStates([1, 1])
   print(bs.n_components)  # 2

   # Trefoil: knot
   bs = BraidStates([1, 1, 1])
   print(bs.n_components)  # 1

After computing, the number of components is also available in
``result["metadata"]["components"]``.

The Hopf Link
--------------

The Hopf link :math:`L2a1` is the closure of :math:`\sigma_1^2`, encoded
as ``[1, 1]``.

.. code-block:: python

   from fkcompute import fk

   result = fk([1, 1], degree=2)

   print("Components:", result["metadata"]["components"])  # 2
   print("x-variables:", result["metadata"]["num_x_variables"])  # 2

The FK polynomial now has two x-variables:

.. code-block:: python

   from fkcompute.output.symbolic import matrix_to_polynomial

   poly = matrix_to_polynomial(result)
   print(poly)   # polynomial in x, y, q  (or x1, x2, q for more components)

Variable naming:

* 1 component → ``x``
* 2 components → ``x``, ``y``
* 3+ components → ``a``, ``b``, ``c``, ... (skipping ``q``)

CLI:

.. code-block:: bash

   fk simple "[1,1]" 2 --symbolic

The Whitehead Link
-------------------

The Whitehead link is a 2-component link with linking number 0:

.. code-block:: python

   # Whitehead link braid representation
   result = fk([1, -2, 1, -2, 1], degree=3)
   print("Components:", result["metadata"]["components"])  # 2

The Borromean Rings
--------------------

The Borromean rings form a 3-component link:

.. code-block:: python

   # One standard braid representation of the Borromean rings
   result = fk([1, -2, 1, -2, 3, -2, 3, -2], degree=2)
   print("Components:", result["metadata"]["components"])  # 3

For 3-component links, the FK polynomial has three variables (``a``, ``b``,
``c`` in symbolic output):

.. code-block:: python

   from fkcompute.output.symbolic import format_symbolic_output

   print(format_symbolic_output(result, "pretty"))

Iterating Over Terms
---------------------

Each term in ``result["terms"]`` has an ``x`` list of length equal to the
number of components.  The :math:`i`-th element is the exponent of the
:math:`i`-th component's variable:

.. code-block:: python

   result = fk([1, 1], degree=2)   # Hopf link

   for term in result["terms"]:
       x1_power = term["x"][0]
       x2_power = term["x"][1]
       for qt in term["q_terms"]:
           coeff = int(qt["c"])
           q_pow = qt["q"]
           if coeff != 0:
               print(f"  {coeff} * x^{x1_power} * y^{x2_power} * q^{q_pow}")

Batch Computation for Links
----------------------------

.. code-block:: yaml

   # links.yaml
   computations:
     - name: hopf
       braid: [1, 1]
       degree: 2
       symbolic: true

     - name: whitehead
       braid: [1, -2, 1, -2, 1]
       degree: 3
       preset: parallel

.. code-block:: python

   from fkcompute import fk
   results = fk("links.yaml")

   for name, res in results.items():
       n = res["metadata"]["components"]
       print(f"{name}: {n} component(s), {len(res['terms'])} terms")
