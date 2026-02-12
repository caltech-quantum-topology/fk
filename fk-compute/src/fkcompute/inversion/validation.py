"""
Sign assignment feasibility checking for FK computation.

This module provides functions for validating sign assignments
against constraint relations during the Phase 1 inversion search.
"""

from typing import Dict, List, Optional

from ..domain.solver.assignment import symbolic_variable_assignment
from ..domain.solver.symbolic_constraints import process_assignment


def check_sign_assignment(degree: int, relations: List, braid_states) -> Optional[Dict]:
    """
    Check if a sign assignment is valid for a given degree.

    Parameters
    ----------
    degree
        Degree bound to check.
    relations
        List of reduced constraint relations.
    braid_states
        BraidStates object with sign assignment.

    Returns
    -------
    dict or None
        Dictionary with criteria, multiples, single_signs, and assignment if valid,
        None otherwise.
    """
    from ..solver.ilp import integral_bounded

    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations)
    for value in criteria.values():
        multiples.append(degree - value)
    if not integral_bounded(multiples, singlesigns):
        return None
    return {
        "criteria": criteria,
        "multiples": multiples,
        "single_signs": singlesigns,
        "assignment": assignment,
    }
