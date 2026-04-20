Batch Computations
==================

This tutorial shows how to compute FK invariants for many braids in a single
invocation, collect the results programmatically, and export them for further
analysis.

Batch Config Files
-------------------

A config file with a top-level ``computations`` list runs multiple braids in
sequence and returns a dictionary keyed by computation name.

.. code-block:: yaml

   # batch.yaml
   preset: parallel
   save_data: true
   save_dir: results

   computations:
     - name: trefoil
       braid: [1, 1, 1]
       degree: 2
       symbolic: true

     - name: figure_eight
       braid: [1, -2, 1, -2]
       degree: 3

     - name: torus_knot_t25
       braid: [1, 1, 1, 1, 1]
       degree: 4

     - name: hopf_link
       braid: [1, 1]
       degree: 2
       symbolic: true
       format: latex

Run from the CLI:

.. code-block:: bash

   fk config batch.yaml

Run from Python:

.. code-block:: python

   from fkcompute import fk

   results = fk("batch.yaml")
   print(results.keys())
   # dict_keys(['trefoil', 'figure_eight', 'torus_knot_t25', 'hopf_link'])

Processing Results
-------------------

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import matrix_to_polynomial

   results = fk("batch.yaml")

   for name, result in results.items():
       n_components = result["metadata"]["components"]
       n_terms = len(result["terms"])
       print(f"{name}: {n_components} component(s), {n_terms} terms")

       # Symbolic representation (requires SymPy)
       poly = matrix_to_polynomial(result)
       print(f"  poly: {poly}\n")

Programmatic Batch Invocation
------------------------------

You can also drive a batch loop directly in Python without a config file:

.. code-block:: python

   from fkcompute import fk
   from fkcompute.output.symbolic import format_symbolic_output

   knot_table = {
       "trefoil":       ([1, 1, 1],            2),
       "figure_eight":  ([1, -2, 1, -2],       3),
       "5_1":           ([1, 1, 1, 1, 1],       3),
       "5_2":           ([1, 1, -2, 1, -2],     3),
       "6_2":           ([1, 1, -2, 1, 1, -2],  3),
   }

   all_results = {}
   for name, (braid, degree) in knot_table.items():
       print(f"Computing {name}...", flush=True)
       all_results[name] = fk(braid, degree, threads=2)

   # Print LaTeX table
   for name, result in all_results.items():
       latex = format_symbolic_output(result, "latex")
       print(f"  {name}: ${latex}$")

Exporting to JSON
-----------------

.. code-block:: python

   import json
   from fkcompute import fk

   results = fk("batch.yaml")

   with open("all_results.json", "w") as f:
       json.dump(results, f, indent=2)

Loading and Reformatting Later
-------------------------------

.. code-block:: python

   import json
   from fkcompute.output.symbolic import format_symbolic_output

   with open("results/trefoil.json") as f:
       result = json.load(f)

   print(format_symbolic_output(result, "mathematica"))

Or with the CLI:

.. code-block:: bash

   fk print-as results/trefoil.json --format mathematica

Multiple Config Files in One CLI Call
--------------------------------------

Pass multiple config files to ``fk config``; they are executed in order:

.. code-block:: bash

   fk config knots.yaml links.yaml torus_knots.yaml

.. note::

   Each config file is independent.  Global defaults defined in one file do
   not carry over to the next.

Performance Tips for Batch Jobs
---------------------------------

* Use ``preset: parallel`` to make full use of available CPUs.
* For large tables, consider distributing files across multiple processes or
  machines using your own orchestration (e.g. ``xargs``, GNU Parallel,
  or a workflow manager like Snakemake).
* Pre-compute and cache inversion data for braids that appear in multiple
  computations at different degrees (use ``inversion_file``).
* Use ``save_data: true`` to preserve intermediate files; then resume failed
  jobs by supplying ``ilp_file`` directly.

Example: Resuming a Failed Job
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the C++ computation fails or times out, you can resume from the ILP file:

.. code-block:: yaml

   braid: [1, -2, 3, -4, 3, -2]
   degree: 6
   ilp_file: results/my_knot_ilp.csv
   threads: 8
