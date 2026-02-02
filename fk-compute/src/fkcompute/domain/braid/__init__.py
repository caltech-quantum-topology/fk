"""
Braid-related domain objects.

This subpackage contains:
- types: StateLiteral and state constants
- word: Braid word operations
- states: BraidStates class for braid state management
"""

from .types import StateLiteral, ZERO_STATE, NUNITY_STATE
from .word import is_positive_braid, is_homogeneous_braid
from .states import BraidStates

__all__ = [
    "StateLiteral",
    "ZERO_STATE",
    "NUNITY_STATE",
    "is_positive_braid",
    "is_homogeneous_braid",
    "BraidStates",
]
