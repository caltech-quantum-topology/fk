Configuration Files
===================

Configuration files allow you to specify all computation parameters in a
reproducible, shareable format.  They are the recommended approach for
non-trivial computations and batch jobs.

Supported Formats
-----------------

* **YAML** — ``.yaml`` / ``.yml``  (requires ``fkcompute[yaml]``)
* **JSON** — ``.json``  (built-in, no extra dependency)

.. code-block:: bash

   fk config my_run.yaml
   fk config my_run.json

Or from Python:

.. code-block:: python

   from fkcompute import fk
   result = fk("my_run.yaml")

Two Modes
---------

Single Computation
~~~~~~~~~~~~~~~~~~

A config file with a top-level ``braid`` and ``degree`` key runs a single
computation and returns one result dictionary.

.. code-block:: yaml

   braid: [1, 1, 1]
   degree: 2

Batch Mode
~~~~~~~~~~

A config file with a top-level ``computations`` key runs multiple
computations and returns a dictionary keyed by computation name.

.. code-block:: yaml

   computations:
     - name: trefoil
       braid: [1, 1, 1]
       degree: 2
     - name: figure_eight
       braid: [1, -2, 1, -2]
       degree: 3

Global defaults placed at the top level of a batch config apply to every
entry unless overridden:

.. code-block:: yaml

   # Global defaults
   preset: parallel
   save_data: true
   save_dir: results

   computations:
     - name: trefoil
       braid: [1, 1, 1]
       degree: 2
       symbolic: true            # local override
     - name: figure_eight
       braid: [1, -2, 1, -2]
       degree: 3
       preset: single thread     # local override of global preset

All Configuration Keys
-----------------------

Required Keys
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Key
     - Type
     - Description
   * - ``braid``
     - list[int]
     - Braid word as a list of signed generator indices, e.g. ``[1,-2,3]``.
   * - ``degree``
     - int
     - Degree truncation for the FK computation.

Optional Keys — Computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Default
     - Description
   * - ``preset``
     - ``null``
     - Named preset (``"single thread"`` or ``"parallel"``).  Preset values
       are applied first; explicit keys override them.
   * - ``threads``
     - ``1``
     - Number of C++ compute threads (``fk_main --threads``).
   * - ``max_workers``
     - ``1``
     - Python multiprocessing workers for the inversion search.
   * - ``chunk_size``
     - ``16384``
     - Permutation candidates per worker task.
   * - ``include_flip``
     - ``false``
     - Include horizontally flipped braid variants in the inversion search.
   * - ``max_shifts``
     - ``null``
     - Maximum number of cyclic shifts to try during inversion (``null``
       means unlimited).
   * - ``partial_signs``
     - ``null``
     - Partial sign constraints for the inversion search (see below).

Optional Keys — Output
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Default
     - Description
   * - ``symbolic``
     - ``false``
     - Compute and attach a symbolic polynomial representation using SymPy.
   * - ``format``
     - ``pretty``
     - Symbolic output format.  Choices: ``pretty``, ``inline``, ``latex``,
       ``mathematica``, ``str``.
   * - ``save_data``
     - ``false``
     - Write intermediate files and the final JSON result to disk.
   * - ``save_dir``
     - ``"data"``
     - Directory for saved files.
   * - ``name``
     - auto
     - Base filename for saved files.  Defaults to a timestamp string.

Optional Keys — Precomputed Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can skip expensive phases by supplying precomputed data:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Default
     - Description
   * - ``inversion``
     - ``null``
     - Inline inversion data (component-index → sign list mapping).  Skips
       Phase 1 entirely.
   * - ``inversion_file``
     - ``null``
     - Path to a JSON file containing saved inversion data.
   * - ``ilp_file``
     - ``null``
     - Path to a pre-generated ILP ``.csv`` file.  Skips both Phase 1 and
       Phase 2.

Optional Keys — Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Default
     - Description
   * - ``verbose``
     - ``false``
     - Enable verbose logging from the Python layer.
   * - ``name``
     - auto
     - Human-readable label shown in log output and used as the base
       filename when ``save_data`` is enabled.


Presets
-------

Two built-in presets provide sensible defaults for different scenarios:

``single thread``
~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   preset: single thread

Equivalent to:

.. code-block:: yaml

   max_workers: 1
   threads: 1
   chunk_size: 4096
   include_flip: false
   verbose: false
   save_data: false

Use this for lightweight computations or when resources are shared.

``parallel``
~~~~~~~~~~~~~

.. code-block:: yaml

   preset: parallel

Equivalent to (values depend on the machine's CPU count):

.. code-block:: yaml

   max_workers: <cpu_count>
   threads: <cpu_count>
   chunk_size: 16384
   include_flip: false
   verbose: true
   save_data: false

Use this to make full use of your hardware.


Precomputed Inputs
-------------------

Inversion Data (Inline)
~~~~~~~~~~~~~~~~~~~~~~~

Supply the sign assignment directly to skip the Phase 1 search:

.. code-block:: yaml

   braid: [1, -2, 1, -2]
   degree: 11
   inversion:
     0: [-1, 1, -1, 1, -1, 1, -1, 1]

The key is a component index (zero-based); the value is the list of
:math:`\pm 1` signs for each strand segment of that component.

Inversion Data (from File)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Load inversion data saved from a previous run:

.. code-block:: yaml

   braid: [1, -2, 1, -2]
   degree: 11
   inversion_file: data/figure_eight_inversion.json

The JSON file format is:

.. code-block:: json

   {
     "inversion_data": {"0": [1, -1, 1, -1]},
     "braid": [1, -2, 1, -2],
     "degree": 11
   }

ILP File
~~~~~~~~

Skip both inversion and ILP generation by supplying a pre-built ILP file:

.. code-block:: yaml

   braid: [1, 1, 1]
   degree: 10
   ilp_file: data/trefoil_ilp.csv
   threads: 8

.. note::

   The C++ backend expects the file **base path** (without the ``.csv``
   extension).  The Python layer derives this automatically from the
   ``ilp_file`` path.


Partial Signs (Advanced)
--------------------------

You can constrain the inversion search by fixing some strand signs while
leaving others free:

.. code-block:: yaml

   braid: [1, -2, 3, -4]
   degree: 3
   partial_signs:
     0: [1, 0, -1, 0]    # 1/-1 fixed, 0 means unknown
     1: [0, 0, 0, 0]     # all unknown for component 1

This reduces the search from :math:`2^n` to :math:`2^{n - \text{fixed}}`
candidates.


Complete Example
-----------------

The following config file demonstrates all major options in a batch job:

.. code-block:: yaml

   # Batch run: compute FK for several knots
   preset: parallel
   save_data: true
   save_dir: results
   symbolic: false

   computations:

     - name: trefoil_d2
       braid: [1, 1, 1]
       degree: 2
       symbolic: true

     - name: figure_eight_d3
       braid: [1, -2, 1, -2]
       degree: 3
       preset: single thread      # override for this entry

     - name: torus_knot_d5
       braid: [1, 1, 1, 1, 1]
       degree: 5
       threads: 8

     - name: hopf_link
       braid: [1, 1]
       degree: 2
       symbolic: true
       format: latex

Run this with:

.. code-block:: bash

   fk config batch.yaml
