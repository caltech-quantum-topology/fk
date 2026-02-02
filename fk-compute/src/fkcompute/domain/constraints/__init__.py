"""
Constraint system for FK computation.

This subpackage contains:
- relations: Constraint classes (Leq, Less, Zero, Nunity, Alias, Conservation)
- symbols: Symbol class for linear algebra
- reduction: Constraint propagation and reduction
"""

from .relations import Leq, Less, Zero, Nunity, Alias, Conservation
from .symbols import Symbol, symbols, one, zero, nunity, solve
from .reduction import full_reduce, reduce_relations

__all__ = [
    # Relations
    "Leq",
    "Less",
    "Zero",
    "Nunity",
    "Alias",
    "Conservation",
    # Symbols
    "Symbol",
    "symbols",
    "one",
    "zero",
    "nunity",
    "solve",
    # Reduction
    "full_reduce",
    "reduce_relations",
]
