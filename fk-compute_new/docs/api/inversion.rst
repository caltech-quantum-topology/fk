``fkcompute.inversion`` — Sign Assignment Search
=================================================

The ``inversion`` subpackage implements Phase 1 of the computation pipeline:
finding a valid :math:`\pm 1` sign assignment for the strands of a braid.

Entry Point
-----------

.. autofunction:: fkcompute.inversion.api.find_sign_assignment

Result Type
-----------

.. autoclass:: fkcompute.inversion.api.InversionResult
   :members:
   :undoc-members:
   :no-index:

   Example:

   .. code-block:: python

      from fkcompute.inversion.api import find_sign_assignment

      result = find_sign_assignment([1, -2, 1, -2], degree=3, max_workers=2)

      if result.success:
          print("Braid type:", result.braid_type)         # "fibered"
          print("Sign assignment:", result.sign_assignment)
          print("Braid used:", result.braid)

          # Convert to the legacy dict format used by _fk_compute
          d = result.to_result_dict()
      else:
          print("Search failed — try a higher degree or more workers")

Braid Type Detection
---------------------

.. automodule:: fkcompute.inversion.search
   :members: braid_type, BraidType
   :undoc-members:

Parallel Search
---------------

The full parallel search is exposed through
:func:`fkcompute.inversion.search.parallel_try_sign_assignments`.
It returns ``False`` if no valid assignment was found, or a
``(braid, sign_assignment)`` tuple on success.

.. code-block:: python

   from fkcompute.domain.braid.states import BraidStates
   from fkcompute.inversion.search import parallel_try_sign_assignments

   bs = BraidStates([1, -2, 1, -2])
   sol = parallel_try_sign_assignments(
       degree=3,
       braid_states=bs,
       partial_signs=None,
       max_workers=4,
       chunk_size=1 << 12,
       include_flip=False,
   )
   if sol:
       braid, signs = sol
       print(signs)

.. autofunction:: fkcompute.inversion.search.parallel_try_sign_assignments
   :no-index:

Braid Variants
--------------

.. automodule:: fkcompute.inversion.variants
   :members:
   :undoc-members:

Permutation Candidates
----------------------

.. automodule:: fkcompute.inversion.permutations
   :members:
   :undoc-members:

Sign Validation
---------------

.. automodule:: fkcompute.inversion.validation
   :members:
   :undoc-members:
