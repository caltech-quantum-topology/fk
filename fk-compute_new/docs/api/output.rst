``fkcompute.output`` — Symbolic Formatting
==========================================

The ``output`` subpackage converts the raw numeric result dictionary into
human-readable symbolic polynomial expressions using SymPy.

.. note::

   All functions in this module require SymPy.  Install it with
   ``pip install "fkcompute[symbolic]"``.  The module-level constant
   :data:`~fkcompute.SYMPY_AVAILABLE` indicates whether SymPy is installed.

``fkcompute.output.symbolic``
------------------------------

.. automodule:: fkcompute.output.symbolic
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
~~~~~~~~~~~~~~

High-level formatting:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import format_symbolic_output

   result = fk([1, 1, 1], 2)

   print(format_symbolic_output(result, "pretty"))
   print(format_symbolic_output(result, "latex"))
   print(format_symbolic_output(result, "mathematica"))

Full SymPy expression:

.. code-block:: python

   from fkcompute.output.symbolic import matrix_to_polynomial
   import sympy as sp

   poly = matrix_to_polynomial(result)
   print(type(poly))        # <class 'sympy.core.add.Add'>
   print(poly.free_symbols) # {x, q}

   # Evaluate at x=1
   x = sp.Symbol("x")
   print(poly.subs(x, 1))

Q-polynomial from coefficient list:

.. code-block:: python

   from fkcompute.output.symbolic import list_to_q_polynomial

   poly_q = list_to_q_polynomial([[0, -1], [2, 1]])
   print(poly_q)   # q**2 - 1

Format Reference
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``format_type``
     - Description
   * - ``"pretty"``
     - Multi-line Unicode (uses ``sympy.pretty``).
   * - ``"latex"``
     - LaTeX string (uses ``sympy.latex``).
   * - ``"mathematica"``
     - Wolfram Language syntax (uses ``sympy.printing.mathematica``).
   * - ``"inline"``
     - One-line Python string (uses ``str``).
   * - ``"str"``
     - Alias for ``"inline"``.
