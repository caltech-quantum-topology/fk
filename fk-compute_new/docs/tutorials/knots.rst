Computing FK for Knots
======================

This tutorial computes the FK invariant for several classical knots
and demonstrates the main output options.

The Trefoil Knot
-----------------

The trefoil knot is the closure of the 3-braid
:math:`\sigma_1^3 = \sigma_1 \sigma_1 \sigma_1`, encoded as ``[1, 1, 1]``.

.. code-block:: python

   from fkcompute import fk

   result = fk([1, 1, 1], degree=2)

   print("Components:", result["metadata"]["components"])  # 1
   print("Number of terms:", len(result["terms"]))

   for term in result["terms"]:
       x_pow = term["x"]
       for qt in term["q_terms"]:
           print(f"  {qt['c']} * x^{x_pow} * q^{qt['q']}")

CLI equivalent:

.. code-block:: bash

   fk simple "[1,1,1]" 2

Symbolic Output
~~~~~~~~~~~~~~~

Install SymPy (``pip install "fkcompute[symbolic]"``) to obtain a
human-readable polynomial:

.. code-block:: python

   result = fk([1, 1, 1], degree=2, symbolic=True)
   print(result["metadata"]["symbolic"])

.. code-block:: bash

   fk simple "[1,1,1]" 2 --symbolic
   fk simple "[1,1,1]" 2 --format latex

Or convert manually:

.. code-block:: python

   from fkcompute.output.symbolic import matrix_to_polynomial, format_symbolic_output

   # Get a SymPy expression
   poly = matrix_to_polynomial(result)
   print(poly)

   # Format as LaTeX
   print(format_symbolic_output(result, "latex"))

Saving Results
~~~~~~~~~~~~~~

Use ``save_data=True`` to write intermediate and result files:

.. code-block:: python

   result = fk([1, 1, 1], degree=2, save_data=True, name="trefoil_d2")

This writes:

* ``data/trefoil_d2_inversion.json``
* ``data/trefoil_d2_ilp.csv``
* ``data/trefoil_d2.json``

Reformat the saved result at a later time without recomputing:

.. code-block:: bash

   fk print-as data/trefoil_d2.json --format mathematica


The Figure-Eight Knot
----------------------

The figure-eight knot :math:`4_1` is the closure of
:math:`\sigma_1 \sigma_2^{-1} \sigma_1 \sigma_2^{-1}`,
encoded as ``[1, -2, 1, -2]``.

.. code-block:: python

   result = fk([1, -2, 1, -2], degree=3)
   print("Components:", result["metadata"]["components"])  # 1


Speeding Up with Parallel Workers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For non-homogeneous braids at higher degrees, use multiple workers:

.. code-block:: python

   result = fk(
       [1, -2, 1, -2],
       degree=5,
       max_workers=4,   # Python inversion workers
       threads=4,       # C++ compute threads
   )

Or via a config file:

.. code-block:: yaml

   braid: [1, -2, 1, -2]
   degree: 5
   preset: parallel


Torus Knots
-----------

The torus knot :math:`T(2, n)` is the closure of :math:`\sigma_1^n`.

.. code-block:: python

   # T(2,3) = trefoil
   result_t23 = fk([1, 1, 1], degree=2)

   # T(2,5) = cinquefoil (5_1)
   result_t25 = fk([1, 1, 1, 1, 1], degree=3)

   # T(2,7)
   result_t27 = fk([1, 1, 1, 1, 1, 1, 1], degree=4)

These are all positive homogeneous braids; Phase 1 completes instantly.

The torus knot :math:`T(3, 4)` uses a 3-strand braid:

.. code-block:: python

   # T(3,4) = 8_19 in Rolfsen table
   result = fk([1, 2, 1, 2, 1, 2, 1, 2], degree=3)

Comparing Multiple Degrees
---------------------------

Compare results at several truncation degrees to see convergence:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import matrix_to_polynomial, collect_by_variables

   braid = [1, 1, 1]

   for deg in [1, 2, 3]:
       result = fk(braid, degree=deg)
       poly = matrix_to_polynomial(result)
       print(f"degree={deg}: {collect_by_variables(poly)}")

Using a Configuration File
---------------------------

The same computations can be run from a YAML file:

.. code-block:: yaml

   # knots.yaml
   computations:
     - name: trefoil_d2
       braid: [1, 1, 1]
       degree: 2
       symbolic: true

     - name: figure_eight_d3
       braid: [1, -2, 1, -2]
       degree: 3
       preset: parallel

     - name: torus_t25_d4
       braid: [1, 1, 1, 1, 1]
       degree: 4

.. code-block:: bash

   fk config knots.yaml

.. code-block:: python

   from fkcompute import fk
   results = fk("knots.yaml")
   print(results.keys())  # dict_keys(['trefoil_d2', 'figure_eight_d3', 'torus_t25_d4'])
