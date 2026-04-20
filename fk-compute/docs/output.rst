Output Format
=============

Every FK computation returns a single Python dictionary (or a dictionary of
dictionaries in batch mode) with a uniform structure.

Top-Level Result
-----------------

Single computation:

.. code-block:: python

   {
     "terms": [...],       # sparse polynomial terms
     "metadata": {...},    # auxiliary information
   }

Batch computation:

.. code-block:: python

   {
     "trefoil":      {"terms": [...], "metadata": {...}},
     "figure_eight": {"terms": [...], "metadata": {...}},
   }

``terms``
---------

The ``terms`` list encodes the FK polynomial as a sparse collection of
monomials.  Each monomial is a product of a x-variable term
:math:`x_1^{a_1} \cdots x_\ell^{a_\ell}` and a polynomial in the quantum
parameter :math:`q`:

.. code-block:: python

   {
     "x": [a1, a2, ...],        # exponents of the x-variables
     "q_terms": [
       {"q": k1, "c": "n1"},    # n1 * q^k1
       {"q": k2, "c": "n2"},    # n2 * q^k2
       ...
     ]
   }

Field descriptions:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Field
     - Type
     - Description
   * - ``x``
     - list[int]
     - Exponents for each x-variable.  Length equals
       ``metadata["num_x_variables"]``.
   * - ``q_terms``
     - list
     - Sparse list of :math:`(k, c)` pairs where :math:`c \cdot q^k` is a
       term of the q-polynomial coefficient.
   * - ``q_terms[i]["q"]``
     - int
     - Exponent of :math:`q` (may be negative).
   * - ``q_terms[i]["c"]``
     - string
     - Integer coefficient as a decimal string (arbitrary precision).

.. note::

   Coefficients are stored as strings by the C++ backend to preserve
   arbitrary precision.  Convert to Python integers with ``int(c)`` or use
   the :func:`~fkcompute.output.symbolic.matrix_to_polynomial` helper.

Example — Trefoil Knot
~~~~~~~~~~~~~~~~~~~~~~

A typical result for the trefoil at degree 2 looks like:

.. code-block:: python

   {
     "terms": [
       {
         "x": [0],
         "q_terms": [{"q": -2, "c": "1"}, {"q": 0, "c": "-1"}]
       },
       {
         "x": [2],
         "q_terms": [{"q": 4, "c": "1"}]
       }
     ],
     "metadata": { ... }
   }

This encodes the polynomial

.. math::

   F_K(x, q) = (q^{-2} - 1) + x^2 q^4 + \cdots

``metadata``
-------------

The ``metadata`` dictionary contains fields written by both the C++ backend
and the Python layer.

Fields Written by the C++ Backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``num_x_variables``
     - int
     - Number of x-variables (:math:`= \ell`, the number of link
       components).
   * - ``max_x_degrees``
     - list[int]
     - Maximum :math:`x`-variable exponent actually appearing, one per
       variable.
   * - ``storage_type``
     - string
     - Internal coefficient storage format, e.g. ``"flint"``.

Fields Added by the Python Layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``braid``
     - list[int]
     - The braid word actually used for the computation (may differ from the
       input if a cyclic shift was applied during inversion).
   * - ``inversion``
     - dict or null
     - The inversion data (component → sign list), or ``null`` when an ILP
       file was supplied directly.
   * - ``components``
     - int
     - Number of link components, computed from the braid topology.
   * - ``symbolic``
     - string
     - Pretty-printed polynomial string (present only when ``symbolic=True``
       and SymPy is installed).


Symbolic Output
---------------

When ``symbolic=True`` (or ``--symbolic`` on the CLI), fkcompute uses
SymPy to convert the sparse term list into a human-readable polynomial and
attaches it as ``metadata["symbolic"]``.

Output Formats
~~~~~~~~~~~~~~

The symbolic output supports five format types:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Format
     - Description
   * - ``pretty``
     - Multi-line Unicode pretty-print (SymPy's ``sp.pretty``).
   * - ``inline``
     - One-line Python-style string (SymPy's ``str``).
   * - ``latex``
     - LaTeX string (e.g. ``x^{2} q^{4} - 1``).
   * - ``mathematica``
     - Mathematica/Wolfram Language syntax.
   * - ``str``
     - Alias for ``inline``.

CLI:

.. code-block:: bash

   fk simple "[1,1,1]" 2 --symbolic
   fk simple "[1,1,1]" 2 --format latex
   fk simple "[1,1,1]" 2 --format mathematica

Python:

.. code-block:: python

   from fkcompute.output.symbolic import format_symbolic_output

   # result is a previously computed result dict
   print(format_symbolic_output(result, "latex"))

Variable Naming
~~~~~~~~~~~~~~~

Variable names for x-variables are chosen automatically:

* **1 component:** ``x``
* **2 components:** ``x``, ``y``
* **3+ components:** ``a``, ``b``, ``c``, ... (skipping ``q``)

Converting to SymPy
~~~~~~~~~~~~~~~~~~~~

Use :func:`~fkcompute.output.symbolic.matrix_to_polynomial` to obtain a
full SymPy expression object:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import matrix_to_polynomial

   result = fk([1, 1, 1], 2)
   poly = matrix_to_polynomial(result)
   print(poly)              # SymPy expression
   print(poly.free_symbols) # {x, q}


Saved Files
-----------

When ``save_data=True`` (or ``--save`` on the CLI), three files are written
to the save directory:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Contents
   * - ``<name>_inversion.json``
     - Inversion data (sign assignment + braid metadata).
   * - ``<name>_ilp.csv``
     - ILP constraint matrix consumed by the C++ backend.
   * - ``<name>.json``
     - Final FK result (terms + metadata), augmented by Python.

Reformat Without Recomputing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you saved the result, you can apply a different symbolic format at any
time without running the computation again:

.. code-block:: bash

   fk print-as data/trefoil.json --format mathematica

Or from Python:

.. code-block:: python

   import json
   from fkcompute.output.symbolic import format_symbolic_output

   with open("data/trefoil.json") as f:
       result = json.load(f)

   print(format_symbolic_output(result, "latex"))
