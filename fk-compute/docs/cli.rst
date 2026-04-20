Command-Line Interface
======================

The **fkcompute** CLI is invoked as ``fk``.  It is built with
`Typer <https://typer.tiangolo.com/>`_ and provides several sub-commands.

.. code-block:: bash

   fk --help

Global Help
-----------

.. code-block:: text

   Usage: fk [OPTIONS] COMMAND [ARGS]...

   Options:
     --help   Show this message and exit.

   Commands:
     simple       Quick FK computation for a single braid.
     interactive  Interactive computation wizard.
     config       Run one or more configuration files.
     print-as     Reformat a saved JSON result without recomputing.
     template     Create or show configuration file templates.

.. note::

   Invoking ``fk`` with no arguments launches the **interactive wizard**
   (requires ``fkcompute[interactive]``).

   The legacy shorthand ``fk "[braid]" degree`` (without a sub-command) is
   still accepted and is equivalent to ``fk simple "[braid]" degree``.


``fk simple``
--------------

Compute the FK invariant for a single braid from the command line.

.. code-block:: bash

   fk simple [OPTIONS] BRAID DEGREE

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Argument
     - Type
     - Description
   * - ``BRAID``
     - string
     - Braid word. Accepted formats: ``"[1,-2,3]"``, ``"1,-2,3"``,
       ``"1 -2 3"``.  Always quote this argument.
   * - ``DEGREE``
     - integer
     - Degree truncation for the computation.

**Options**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Default
     - Description
   * - ``--symbolic``
     - off
     - Print the result as a symbolic polynomial (requires SymPy).
   * - ``--format FORMAT``
     - ``pretty``
     - Symbolic output format.  Choices: ``pretty``, ``inline``,
       ``latex``, ``mathematica``, ``str``.
   * - ``--threads N``
     - 1
     - Number of threads for the C++ backend.
   * - ``--workers N``
     - 1
     - Number of Python multiprocessing workers for the inversion search.
   * - ``--save``
     - off
     - Save intermediate and result files to ``data/``.
   * - ``--save-dir DIR``
     - ``data``
     - Directory for saved files (used with ``--save``).
   * - ``--name NAME``
     - auto
     - Base filename when saving (used with ``--save``).
   * - ``--verbose``
     - off
     - Enable verbose logging.
   * - ``--help``
     - —
     - Show help and exit.

**Examples**

.. code-block:: bash

   # Trefoil knot at degree 2
   fk simple "[1,1,1]" 2

   # Figure-eight knot with symbolic LaTeX output
   fk simple "[1,-2,1,-2]" 3 --format latex

   # Use 4 threads and save result files
   fk simple "[1,-2,3,-4,3,-2]" 4 --threads 4 --save --name my_knot

   # Hopf link (2-component link)
   fk simple "[1,1]" 2


``fk interactive``
-------------------

Launch the interactive computation wizard.

.. code-block:: bash

   fk interactive [OPTIONS]

**Options**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Option
     - Default
     - Description
   * - ``--quick``
     - off
     - Skip non-essential prompts and use sensible defaults.

The wizard guides you through:

1. Entering a braid word (with syntax hints).
2. Choosing a degree.
3. Selecting a preset (``single thread`` or ``parallel``).
4. Optionally enabling symbolic output and file saving.
5. Displaying a live progress view during computation.

Requires ``fkcompute[interactive]`` (``pip install "fkcompute[interactive]"``).


``fk config``
--------------

Run one or more configuration files.

.. code-block:: bash

   fk config [OPTIONS] CONFIG_FILES...

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Argument
     - Type
     - Description
   * - ``CONFIG_FILES``
     - paths (one or more)
     - YAML (``.yaml``, ``.yml``) or JSON (``.json``) config files.

**Options**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Option
     - Default
     - Description
   * - ``--help``
     - —
     - Show help and exit.

**Examples**

.. code-block:: bash

   # Single config file
   fk config my_run.yaml

   # Multiple files in one invocation (runs sequentially)
   fk config trefoil.yaml figure_eight.yaml torus_link.json

   # Batch config (top-level 'computations:' list)
   fk config batch_run.yaml

See :doc:`configuration` for the full config file reference.


``fk print-as``
----------------

Reformat a previously saved JSON result file without recomputing.

.. code-block:: bash

   fk print-as [OPTIONS] RESULT_FILE

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Argument
     - Type
     - Description
   * - ``RESULT_FILE``
     - path
     - Path to a saved ``.json`` result file.

**Options**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Default
     - Description
   * - ``--format FORMAT``
     - ``pretty``
     - Output format. Choices: ``pretty``, ``inline``, ``latex``,
       ``mathematica``, ``str``.
   * - ``--help``
     - —
     - Show help and exit.

**Examples**

.. code-block:: bash

   fk print-as data/trefoil.json --format inline
   fk print-as data/trefoil.json --format latex
   fk print-as data/trefoil.json --format mathematica

Requires ``fkcompute[symbolic]``.


``fk template``
----------------

Create or display a configuration file template.

**Sub-commands**

.. code-block:: bash

   fk template create [OPTIONS] OUTPUT_FILE
   fk template show

``fk template create``
~~~~~~~~~~~~~~~~~~~~~~~

Write a fully annotated template config file to ``OUTPUT_FILE``.

.. code-block:: bash

   fk template create my_run.yaml

``fk template show``
~~~~~~~~~~~~~~~~~~~~~

Print the template to stdout.

.. code-block:: bash

   fk template show


Braid Input Formats
--------------------

The ``BRAID`` argument to ``fk simple`` (and the ``braid:`` key in config
files) accepts three equivalent string formats:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Format
     - Example
   * - JSON array
     - ``"[1,-2,3]"``
   * - Comma-separated
     - ``"1,-2,3"``
   * - Space-separated
     - ``"1 -2 3"``

Always wrap the argument in quotes so that shell expansion does not
swallow the leading ``-`` of negative generators.


Man Page
---------

A Unix man page for the ``fk`` command is bundled with the package.
Install and view it with:

.. code-block:: bash

   fk-install-man
   man fk
