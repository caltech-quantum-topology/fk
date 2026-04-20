Mathematica / Wolfram Language Integration
==========================================

A Wolfram Language paclet wrapper is included under ``mathematica/FkCompute/``.
It provides a thin Mathematica interface that delegates all computation to the
Python ``fkcompute`` package via a Python bridge.

Prerequisites
-------------

* **fkcompute** installed and working in a Python environment (``fk --help``
  should succeed).
* **Wolfram Mathematica** 12.0 or later (or a compatible Wolfram Engine).

Installation
------------

Install the paclet from the repository:

.. code-block:: mathematica

   (* From within Mathematica *)
   PacletInstall["path/to/mathematica/FkCompute"]

Or load it directly without installing:

.. code-block:: mathematica

   Get["path/to/mathematica/FkCompute/FkCompute.wl"]

Basic Usage
-----------

.. code-block:: mathematica

   (* Load the package *)
   Needs["FkCompute`"]

   (* Compute FK for the trefoil knot *)
   result = FKCompute[{1, 1, 1}, 2]

   (* result is an Association with "terms" and "metadata" keys *)
   result["metadata"]["components"]    (* 1 *)

   (* Pretty polynomial display *)
   FKPolynomial[{1, 1, 1}, 2]

Symbolic Output from the Python Layer
---------------------------------------

The Python backend can produce LaTeX or Mathematica-syntax strings directly:

From Python:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import format_symbolic_output

   result = fk([1, 1, 1], 2)
   mma_str = format_symbolic_output(result, "mathematica")
   print(mma_str)

From the CLI:

.. code-block:: bash

   fk simple "[1,1,1]" 2 --format mathematica

The resulting string can be pasted directly into a Mathematica notebook.

Reformatting a Saved Result
----------------------------

If you saved an FK result to JSON, you can reformat it for Mathematica at
any time:

.. code-block:: bash

   fk print-as data/trefoil.json --format mathematica

Or from Python:

.. code-block:: python

   import json
   from fkcompute.output.symbolic import format_symbolic_output

   with open("data/trefoil.json") as f:
       result = json.load(f)
   print(format_symbolic_output(result, "mathematica"))

Architecture of the Paclet
----------------------------

The paclet calls:

.. code-block:: bash

   python -m fkcompute.mathematica_bridge <braid_json> <degree>

and parses the JSON output.  The bridge module simply wraps
:func:`fkcompute.fk` and serialises the result.

Troubleshooting
---------------

Python Not Found by Mathematica
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the Python executable path explicitly before loading the paclet:

.. code-block:: mathematica

   FkCompute`PythonPath = "/usr/bin/python3";
   Needs["FkCompute`"]

Wrong Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~

If ``fkcompute`` is installed in a virtual environment, activate it or point
``FkCompute`PythonPath`` to the environment's Python binary:

.. code-block:: mathematica

   FkCompute`PythonPath = "/home/user/myenv/bin/python3";
