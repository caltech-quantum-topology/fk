"""
State type definitions for braid computations.

This module defines the StateLiteral class and special state constants
used throughout the FK computation.
"""


class StateLiteral:
    """
    A literal state value in the braid state system.

    Represents fixed states like zero (0) and negative unity (-1)
    that are used as boundary conditions in the constraint system.
    """

    def __init__(self, state: int):
        self.state = state

    def __repr__(self):
        return f'[{self.state}]'

    def __eq__(self, other):
        if isinstance(other, StateLiteral):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)


# Special state constants
ZERO_STATE = StateLiteral(0)
NUNITY_STATE = StateLiteral(-1)
