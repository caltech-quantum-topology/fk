"""
Phase 2 pipeline entry point for FK computation.

Provides a single clear function that performs the full Phase 2 flow:
raw constraints -> reduction -> symbolic assignment -> criteria extraction.
"""

from typing import List

from .reduction import full_reduce
from .system import ConstraintSystem
from ..solver.assignment import symbolic_variable_assignment
from ..solver.symbolic_constraints import process_assignment


def build_constraint_system(braid_states, raw_relations: List = None) -> ConstraintSystem:
    """
    Phase 2: Build and reduce the constraint system for a signed braid.

    Steps:
    1. Generate raw constraints (boundary, periodicity, sign bounds, crossing relations)
    2. Reduce constraints to minimal form (fixed-point simplification)
    3. Assign symbolic variables to free states
    4. Extract degree criteria and inequalities

    Parameters
    ----------
    braid_states
        A BraidStates (or SignedBraid) object with sign data loaded.
    raw_relations
        Pre-computed raw relations. If None, generated from braid_states.

    Returns
    -------
    ConstraintSystem
        The fully processed constraint system ready for ILP generation.
    """
    if raw_relations is None:
        raw_relations = braid_states.get_state_relations()

    reduced_relations = full_reduce(raw_relations)

    assignment = symbolic_variable_assignment(reduced_relations, braid_states)

    criteria, multi_var_inequalities, single_var_signs = process_assignment(
        assignment, braid_states, reduced_relations,
    )

    return ConstraintSystem(
        assignment=assignment,
        degree_criteria=criteria,
        multi_var_inequalities=multi_var_inequalities,
        single_var_signs=single_var_signs,
        raw_relations=raw_relations,
        reduced_relations=reduced_relations,
    )
