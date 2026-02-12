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
        Dictionary with criteria, multi_var_inequalities, single_var_signs, and assignment if valid,
        None otherwise.
    """
    from ..solver.ilp import integral_bounded

    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multi_var_inequalities, single_var_signs = process_assignment(assignment, braid_states, relations)
    for value in criteria.values():
        multi_var_inequalities.append(degree - value)
    if not integral_bounded(multi_var_inequalities, single_var_signs):
        return None
    return {
        "criteria": criteria,
        "multi_var_inequalities": multi_var_inequalities,
        "single_var_signs": single_var_signs,
        "assignment": assignment,
    }
