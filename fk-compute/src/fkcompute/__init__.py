"""
fk-compute
============

A package for computing the FK invariant of braids via inversion,
ILP reduction, and a compiled helper binary.

This module exposes:

- `fk`   – core function to compute the invariant (from cli.py)
- `main` – command-line entrypoint (from cli.py)

Example
-------
>>> from fk_invariant import fk
>>> result = fk([1, -2, -2, 3], degree=3)
>>> print(result)
{'q^2': 1, 'q^0': -1, ...}
"""

from .cli import fk, main

__all__ = ["fk", "main"]

__version__ = "0.1.0"
