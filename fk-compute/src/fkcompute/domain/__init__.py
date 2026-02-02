"""
Core domain logic for FK computation.

This module contains pure domain logic with no I/O dependencies:
- braid/: Braid-related domain objects
- constraints/: Constraint system and relations
- solver/: Sign assignment logic
"""

from .braid.types import StateLiteral, ZERO_STATE, NUNITY_STATE
from .braid.word import is_positive_braid, is_homogeneous_braid
from .braid.states import BraidStates
from .constraints.relations import Leq, Less, Zero, Nunity, Alias, Conservation
from .constraints.symbols import Symbol, symbols, one, zero, nunity, solve
from .constraints.reduction import full_reduce, reduce_relations

__all__ = [
    # Types
    "StateLiteral",
    "ZERO_STATE",
    "NUNITY_STATE",
    # Braid
    "BraidStates",
    "is_positive_braid",
    "is_homogeneous_braid",
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
