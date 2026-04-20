fkcompute documentation
=======================

**fkcompute** is a Python package and command-line tool for computing the
**FK invariant** :math:`F_K(x, q)` of knots and links presented as braid closures.
It combines a Python orchestration layer with a compiled C++ backend for
high-performance polynomial computation.

.. code-block:: python

   from fkcompute import fk

   # Compute F_K for the trefoil knot at degree 2
   result = fk([1, 1, 1], 2)
   print(result["metadata"]["components"])  # 1

.. code-block:: bash

   # Equivalent CLI invocation
   fk simple "[1,1,1]" 2

----

.. grid:: 2

   .. grid-item-card:: Installation
      :link: installation
      :link-type: doc

      System dependencies, pip install, optional extras, and Gurobi setup.

   .. grid-item-card:: Quickstart
      :link: quickstart
      :link-type: doc

      Get up and running with your first FK computation in minutes.

.. grid:: 2

   .. grid-item-card:: Mathematical Overview
      :link: overview
      :link-type: doc

      What the FK invariant is, where it comes from, and how this package
      computes it.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc

      Step-by-step guides for knots, links, batch jobs, and config files.

.. grid:: 2

   .. grid-item-card:: CLI Reference
      :link: cli
      :link-type: doc

      Complete reference for the ``fk`` command-line tool and all sub-commands.

   .. grid-item-card:: Python API Reference
      :link: api/index
      :link-type: doc

      Full auto-generated API documentation for all public modules.

----

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   installation
   quickstart
   overview

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   tutorials/index
   cli
   configuration
   output

.. toctree::
   :maxdepth: 2
   :caption: Internals
   :hidden:

   pipeline
   cpp_backend

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   :hidden:

   api/index

.. toctree::
   :maxdepth: 1
   :caption: Development
   :hidden:

   development
   troubleshooting
   mathematica

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
