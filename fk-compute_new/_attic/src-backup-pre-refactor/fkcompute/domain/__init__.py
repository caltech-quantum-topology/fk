"""
Core domain logic for FK computation.

This module contains pure domain logic with no I/O dependencies:
- braid/: Braid-related domain objects
- constraints/: Constraint system and relations
- solver/: Sign assignment logic
"""

from .braid.types import StateLiteral, ZERO_STATE, NEG_ONE_STATE
from .braid.word import is_positive_braid, is_homogeneous_braid
from .braid.states import BraidStates
from .constraints.relations import Leq, Less, Zero, NegOne, Alias, Conservation
from .constraints.symbols import Symbol, symbols, one, zero, neg_one, solve
from .constraints.reduction import full_reduce, reduce_relations

__all__ = [
    # Types
    "StateLiteral",
    "ZERO_STATE",
    "NEG_ONE_STATE",
    # Braid
    "BraidStates",
    "is_positive_braid",
    "is_homogeneous_braid",
    # Relations
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
]
