"""
BraidStates class for braid state management.

Backward-compatible wrapper that combines BraidTopology and SignedBraid.
"""

from typing import List

from .topology import BraidTopology
from .signed import SignedBraid


class BraidStates(SignedBraid):
    """
    Represents a braid and its associated state information.

    This is a backward-compatible wrapper that creates a BraidTopology
    and initializes a SignedBraid from it.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    """

    def __init__(self, braid: List[int]):
        topology = BraidTopology(braid)
        super().__init__(topology)
