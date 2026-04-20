``fkcompute.infra`` — Infrastructure
=====================================

The ``infra`` subpackage handles all I/O and external-process concerns: locating
and executing the C++ binary, parsing configuration files, and managing
intermediate files.

``fkcompute.infra.binary``
---------------------------

Locate and execute the ``fk_main`` C++ binary.

**Key functions:**

* :func:`~fkcompute.infra.binary.binary_path` — return the path to the binary
  (searches ``PATH`` first, then the packaged ``_bin/fk_main``).
* :func:`~fkcompute.infra.binary.run_fk_binary` — execute the backend as a
  subprocess.
* :func:`~fkcompute.infra.binary.safe_unlink` — remove a file if it exists.

.. code-block:: python

   from fkcompute.infra.binary import run_fk_binary

   # ilp_file_base = path without .csv extension
   # output_base   = path without .json extension
   run_fk_binary("data/trefoil_ilp", "data/trefoil", threads=4)

.. automodule:: fkcompute.infra.binary
   :members:
   :undoc-members:
   :show-inheritance:

``fkcompute.infra.config``
---------------------------

Parse configuration files and braid input strings.

**Key functions:**

* :func:`~fkcompute.infra.config.load_config_file` — parse a JSON or YAML
  configuration file into a Python dict.
* :func:`~fkcompute.infra.config.load_inversion_file` — load saved inversion
  data from a JSON file.
* :func:`~fkcompute.infra.config.load_ilp_file` — load an ILP specification
  from a CSV file.
* :func:`~fkcompute.infra.config.parse_int_list` — parse a braid string in
  any supported format.

.. code-block:: python

   from fkcompute.infra.config import parse_int_list

   parse_int_list("[1,-2,3]")  # [1, -2, 3]
   parse_int_list("1,-2,3")    # [1, -2, 3]
   parse_int_list("1 -2 3")    # [1, -2, 3]

.. automodule:: fkcompute.infra.config
   :members:
   :undoc-members:
   :show-inheritance:

``fkcompute.infra.io``
-----------------------

.. automodule:: fkcompute.infra.io
   :members:
   :undoc-members:
