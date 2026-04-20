Quickstart
==========

This page shows the fastest path from installation to a working FK
computation.  For a deeper explanation of what the results mean, see
:doc:`overview`.

Braid Word Conventions
-----------------------

Braids are written as a list of **signed generator indices**:

.. math::

   \sigma_1^{+1},\; \sigma_2^{-1},\; \sigma_1^{+1} \;\longmapsto\; [1, -2, 1]

- Each integer :math:`k` represents the standard braid generator
  :math:`\sigma_{|k|}`.
- A **positive** integer means the over-strand goes left-to-right.
- A **negative** integer means the inverse crossing.
- Strands are 0-indexed; generators are 1-indexed up to
  :math:`n_\text{strands} - 1`.

The FK invariant is computed for the **closure** of the braid.

Command-Line Interface
-----------------------

The simplest invocation uses the ``simple`` sub-command:

.. code-block:: bash

   fk simple "[1,1,1]" 2

This computes :math:`F_K(x, q)` for the **trefoil knot** (braid closure of
:math:`\sigma_1^3`) truncated at degree 2.

The braid argument accepts several equivalent formats:

.. code-block:: bash

   fk simple "[1,1,1]" 2        # JSON-array style
   fk simple "1,1,1" 2          # comma-separated
   fk simple "1 1 1" 2          # space-separated

.. warning::

   Always quote the braid argument so that negative integers are not
   misinterpreted as CLI flags:

   .. code-block:: bash

      fk simple "[1,-2,1,-2]" 3   # correct
      fk simple  1,-2,1,-2  3     # may fail – leading dash

Symbolic output (requires ``fkcompute[symbolic]``):

.. code-block:: bash

   fk simple "[1,1,1]" 2 --symbolic          # pretty Unicode
   fk simple "[1,1,1]" 2 --format latex      # LaTeX string
   fk simple "[1,1,1]" 2 --format mathematica

Interactive wizard (requires ``fkcompute[interactive]``):

.. code-block:: bash

   fk                 # full wizard
   fk interactive     # same

Python API
-----------

.. code-block:: python

   from fkcompute import fk

   # Simplest call: braid + degree
   result = fk([1, 1, 1], 2)

The returned dictionary has two top-level keys:

.. code-block:: python

   result["terms"]     # sparse polynomial terms
   result["metadata"]  # auxiliary information

Print the number of link components:

.. code-block:: python

   print(result["metadata"]["components"])   # 1  (trefoil is a knot)

Get all terms:

.. code-block:: python

   for term in result["terms"]:
       x_powers = term["x"]
       for qt in term["q_terms"]:
           print(f"  x^{x_powers} * {qt['c']} * q^{qt['q']}")

Symbolic polynomial (requires SymPy):

.. code-block:: python

   result = fk([1, 1, 1], 2, symbolic=True)
   print(result["metadata"]["symbolic"])

Configuration Files
--------------------

For reproducible computations or batch jobs, use a YAML (or JSON) config file.

Create a template:

.. code-block:: bash

   fk template create my_run.yaml

Minimal ``my_run.yaml``:

.. code-block:: yaml

   braid: [1, 1, 1]
   degree: 2

Run it:

.. code-block:: bash

   fk config my_run.yaml

Or from Python:

.. code-block:: python

   result = fk("my_run.yaml")

See :doc:`configuration` for all available options.

Next Steps
----------

* :doc:`overview` — The mathematics behind the FK invariant.
* :doc:`tutorials/index` — Worked examples for knots, links, and batch jobs.
* :doc:`cli` — Complete CLI reference.
* :doc:`api/index` — Full Python API documentation.
