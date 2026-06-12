"""
ConstraintSystem dataclass for FK computation.

Represents the fully processed constraint system ready for ILP generation.
"""

from dataclasses import dataclass, field
from typing import Dict, List

from .symbols import Symbol


@dataclass
class ConstraintSystem:
    """The fully processed constraint system, ready for ILP generation.

    Attributes
    ----------
    assignment
        Mapping from state locations to their symbolic expressions.
    degree_criteria
        Mapping from component index to its degree constraint expression.
    multi_var_inequalities
        List of symbolic expressions where 0 <= expr must hold.
    single_var_signs
        Mapping from elementary symbols to their signs (+1.0 or -1.0).
    raw_relations
        Original relations before reduction (for debugging).
    reduced_relations
        Relations after fixed-point reduction (for debugging).
    """
    assignment: Dict
    degree_criteria: Dict
    multi_var_inequalities: List[Symbol]
    single_var_signs: Dict
    raw_relations: List = field(default_factory=list)
    reduced_relations: List = field(default_factory=list)
