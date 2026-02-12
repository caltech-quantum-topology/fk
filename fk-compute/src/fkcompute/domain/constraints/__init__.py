"""
Constraint system for FK computation.

This subpackage contains:
- relations: Constraint classes (Leq, Less, Zero, NegOne, Alias, Conservation)
- symbols: Symbol class for linear algebra
- reduction: Constraint propagation and reduction
- system: ConstraintSystem dataclass
- pipeline: Phase 2 pipeline entry point
"""

from .relations import Constraint, Leq, Less, Zero, NegOne, Alias, Conservation
from .symbols import Symbol, symbols, one, zero, neg_one, solve
from .reduction import full_reduce, reduce_relations
from .system import ConstraintSystem
from .pipeline import build_constraint_system

__all__ = [
    # Relations
    "Constraint",
    "Leq",
    "Less",
    "Zero",
    "NegOne",
    "Alias",
    "Conservation",
    # Symbols
    "Symbol",
    "symbols",
    "one",
    "zero",
    "neg_one",
    "solve",
    # Reduction
    "full_reduce",
    "reduce_relations",
    # System
    "ConstraintSystem",
    "build_constraint_system",
]
