``fkcompute.interactive`` — Interactive Wizard
===============================================

The ``interactive`` subpackage provides a Rich-based terminal user interface
for guided FK computation.

.. note::

   This subpackage requires the ``rich`` library.  Install it with
   ``pip install "fkcompute[interactive]"``.

Invoking the Wizard
--------------------

The wizard is launched by running ``fk`` (with no sub-command) or
``fk interactive``:

.. code-block:: bash

   fk                       # full interactive session
   fk interactive --quick   # skip optional prompts

From Python (not normally needed):

.. code-block:: python

   from fkcompute.interactive.wizard import FKWizard

   wizard = FKWizard()
   wizard.run()

Module Reference
-----------------

``fkcompute.interactive.wizard``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.interactive.wizard
   :members:
   :undoc-members:

``fkcompute.interactive.progress``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.interactive.progress
   :members:
   :undoc-members:

   **FKProgressTracker** is passed to the internal ``_fk_compute`` function
   as the ``_progress_callback`` argument.  It receives phase events and
   renders Rich progress panels:

   * ``start_inversion`` / ``complete_inversion``
   * ``start_ilp_generation`` / ``complete_ilp_generation``
   * ``start_fk_computation`` / ``update_fk_progress`` / ``complete_fk_computation``

``fkcompute.interactive.ui``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.interactive.ui
   :members:
   :undoc-members:

``fkcompute.interactive.prompts``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.interactive.prompts
   :members:
   :undoc-members:

``fkcompute.interactive.history``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: fkcompute.interactive.history
   :members:
   :undoc-members:
