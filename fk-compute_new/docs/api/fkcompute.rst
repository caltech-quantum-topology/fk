``fkcompute`` — Top-Level API
=============================

The ``fkcompute`` package exposes a minimal, stable public API at the top
level.  Import everything you need from here.

.. code-block:: python

   from fkcompute import fk, PRESETS, BraidStates

Main Computation Function
--------------------------

.. autofunction:: fkcompute.fk

Batch and Config-File Functions
---------------------------------

.. autofunction:: fkcompute.fk_from_config

.. autofunction:: fkcompute.fk_batch_from_config

Presets
-------

.. data:: fkcompute.PRESETS

   Dictionary of named preset configurations.  Keys are ``"single thread"``
   and ``"parallel"``; values are parameter dicts.

   .. code-block:: python

      from fkcompute import PRESETS
      print(PRESETS.keys())
      # dict_keys(['single thread', 'parallel'])

      print(PRESETS["parallel"])
      # {'max_workers': 8, 'threads': 8, 'chunk_size': 16384, ...}

Domain Types
------------

These are re-exported from the domain submodules for convenience.

.. autoclass:: fkcompute.BraidStates
   :members:
   :undoc-members:
   :show-inheritance:

.. autodata:: fkcompute.ZERO_STATE

.. autodata:: fkcompute.NEG_ONE_STATE

Braid Word Utilities
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: fkcompute.is_positive_braid

.. autofunction:: fkcompute.is_homogeneous_braid

Constraint Classes
------------------

.. autoclass:: fkcompute.Leq
   :members:

.. autoclass:: fkcompute.Less
   :members:

.. autoclass:: fkcompute.Zero
   :members:

.. autoclass:: fkcompute.NegOne
   :members:

.. autoclass:: fkcompute.Alias
   :members:

.. autoclass:: fkcompute.Conservation
   :members:

Symbolic Output
---------------

.. autofunction:: fkcompute.format_symbolic_output

.. autodata:: fkcompute.SYMPY_AVAILABLE

   ``True`` if SymPy is installed; ``False`` otherwise.  Check this before
   calling any symbolic function if you are not sure whether SymPy is
   available.

   .. code-block:: python

      from fkcompute import SYMPY_AVAILABLE, format_symbolic_output

      if SYMPY_AVAILABLE:
          print(format_symbolic_output(result, "latex"))

``fkcompute.api``
-----------------

.. automodule:: fkcompute.api.compute
   :members:
   :undoc-members:

.. automodule:: fkcompute.api.batch
   :members:
   :undoc-members:

.. automodule:: fkcompute.api.presets
   :members:
   :undoc-members:
