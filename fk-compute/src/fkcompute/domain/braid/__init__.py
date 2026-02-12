"""
Braid-related domain objects.

This subpackage contains:
- types: StateLiteral and state constants
- word: Braid word operations
- topology: BraidTopology class for pure topological braid data
- signed: SignedBraid class combining topology with sign assignment
- states: BraidStates backward-compatible wrapper
"""

from .types import StateLiteral, ZERO_STATE, NUNITY_STATE
from .word import is_positive_braid, is_homogeneous_braid
from .topology import BraidTopology
from .signed import SignedBraid
from .states import BraidStates

__all__ = [
    "StateLiteral",
    "ZERO_STATE",
    "NUNITY_STATE",
    "is_positive_braid",
    "is_homogeneous_braid",
    "BraidTopology",
    "SignedBraid",
    "BraidStates",
]
