"""
Constraint relation classes for FK computation.

This module defines the constraint types used in the FK constraint system:
- Leq: Less than or equal inequality
- Less: Strict less than inequality
- Zero: Equality to zero state
- NegOne: Equality to negative one state
- Alias: Equality between two states
- Conservation: Sum conservation constraint
"""

from abc import ABC, abstractmethod
from typing import List, Any, Optional, Tuple

from ..braid.types import StateLiteral


def _sort_any(xs):
    """Sort any list by string representation."""
    return list(sorted(xs, key=lambda x: str(x)))


class Constraint(ABC):
    """Base class for all constraint relations in the FK system."""

    @abstractmethod
    def variables(self) -> List:
        """Return the list of variables involved in this constraint."""
        ...


class Leq(Constraint):
    """Less than or equal constraint: first <= second."""

    def __init__(self, first: Any, second: Any):
        self.first = first
        self.second = second

    def variables(self) -> List:
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} <= {self.second}'

    def __eq__(self, other):
        if isinstance(other, Leq):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))


class Less(Constraint):
    """Strict less than constraint: first < second."""

    def __init__(self, first: Any, second: Any):
        self.first = first
        self.second = second

    def variables(self) -> List:
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} < {self.second}'

    def __eq__(self, other):
        if isinstance(other, Less):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))


class Zero(Constraint):
    """Zero state constraint: state = [0]."""

    def __init__(self, state: Any):
        self.state = state

    def variables(self) -> List:
        return [self.state]

    def __repr__(self):
        return f'Zero {self.state} = [0]'

    def __eq__(self, other):
        if isinstance(other, Zero):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)


class NegOne(Constraint):
    """Negative one state constraint: state = [-1]."""

    def __init__(self, state: Any):
        self.state = state

    def variables(self) -> List:
        return [self.state]

    def __repr__(self):
        return f'NegOne {self.state} = [-1]'

    def __eq__(self, other):
        if isinstance(other, NegOne):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)



class Alias(Constraint):
    """Alias constraint: state and alias are equal."""

    def __init__(self, state: Any, alias: Any):
        state, alias = _sort_any([state, alias])
        self.state = state
        self.alias = alias

    def variables(self) -> List:
        return [self.state, self.alias]

    def __repr__(self):
        return f'Alias {self.alias} := {self.state}'

    def __eq__(self, other):
        if isinstance(other, Alias):
            return self.state == other.state and self.alias == other.alias
        return False

    def __hash__(self):
        return hash((self.state, self.alias))


class Conservation(Constraint):
    """Conservation constraint: sum of inputs = sum of outputs."""

    def __init__(self, inputs: List, outputs: List):
        inputs = _sort_any(inputs)
        outputs = _sort_any(outputs)
        inputs, outputs = _sort_any([inputs, outputs])
        self.inputs = inputs
        self.outputs = outputs

    def variables(self) -> List:
        return self.inputs + self.outputs

    def try_sum_alias(self) -> Optional[Tuple[Any, List]]:
        """
        Try to extract a sum alias from this conservation constraint.

        Returns None if both inputs and outputs have more than one element.
        Otherwise returns (alias, sum_vars) where alias = sum(sum_vars).
        """
        # case a + b = c + d (no variable is bound)
        if len(self.inputs) != 1 and len(self.outputs) != 1:
            return None

        # case a = b + c, return a
        elif len(self.inputs) == 1 and not isinstance(self.inputs[0], StateLiteral):
            return self.inputs[0], self.outputs

        # case a + b = c, return c
        elif len(self.outputs) == 1 and not isinstance(self.outputs[0], StateLiteral):
            return self.outputs[0], self.inputs

        return None

    def __repr__(self):
        input_sum = " + ".join(str(x) for x in self.inputs)
        output_sum = " + ".join(str(x) for x in self.outputs)
        return f'Conservation {input_sum} = {output_sum}'

    def __eq__(self, other):
        if isinstance(other, Conservation):
            return self.inputs == other.inputs and self.outputs == other.outputs
        return False

    def __hash__(self):
        return hash((tuple(self.inputs), tuple(self.outputs)))
