"""
Shared symbolic constraint processing for FK computation.

This module contains functions used by both Phase 1 (inversion/validation)
and Phase 2 (solver/ilp) for processing symbolic constraints:
- _minimum_degree_symbolic: compute degree constraints from crossing data
- _inequality_manager: extract single/multiple variable inequalities
- _process_assignment: combine the above into criteria, multiples, singlesigns
"""

from typing import Dict, List, Any

from ..braid.types import ZERO_STATE, NEG_ONE_STATE
from ..constraints.relations import Leq, Less
from ..constraints.symbols import Symbol, one, zero


def _expr_from_dict(dict_: Dict) -> Any:
    """Build a symbolic expression from a coefficients dictionary."""
    expression = 0
    for key, value in dict_.items():
        expression += value * key
    return expression


def minimum_degree_symbolic(assignment: Dict, braid_states, verbose: bool = False) -> Dict:
    """
    Compute minimum degree constraints symbolically.

    Parameters
    ----------
    assignment
        Symbolic variable assignment.
    braid_states
        BraidStates object.
    verbose
        Whether to print progress.

    Returns
    -------
    dict
        Dictionary mapping components to their degree constraints.
    """
    conditions = {val: zero for val in range(braid_states.n_components)}
    for index in range(0, braid_states.n_strands):
        if verbose:
            print(braid_states.closed_strand_components)
        conditions[braid_states.closed_strand_components[index]] -= 1/2
    for index in range(braid_states.n_crossings):
        crossing_type = braid_states.r_matrices[index]
        in1 = braid_states.top_input_state_locations[index]
        in2 = (in1[0] + 1, in1[1])
        out1 = (in1[0], in1[1] + 1)
        out2 = (out1[0] + 1, out1[1])
        in1 = braid_states.get_state(in1)
        in2 = braid_states.get_state(in2)
        out1 = braid_states.get_state(out1)
        out2 = braid_states.get_state(out2)
        if crossing_type == "R1" or crossing_type == "R2":
            conditions[braid_states.top_crossing_components[index]] += (assignment[out1] + assignment[in2] + 1) / 4
            conditions[braid_states.bottom_crossing_components[index]] += (3 * assignment[out2] - assignment[in1] + 1) / 4
        elif crossing_type == "R3" or crossing_type == "R4":
            conditions[braid_states.top_crossing_components[index]] -= (3 * assignment[in2] - assignment[out1] + 1) / 4
            conditions[braid_states.bottom_crossing_components[index]] -= (assignment[in1] + assignment[out2] + 1) / 4
        else:
            raise Exception("Crossing type is not one of the four acceptable values: 'R1', 'R2', 'R3', or 'R4'.")
    if verbose:
        for value in conditions.values():
            print(value.var)
    return conditions


def inequality_manager(relations: List, assignment: Dict, braid_states):
    """
    Process inequalities from relations into single and multiple variable forms.

    Parameters
    ----------
    relations
        List of constraint relations.
    assignment
        Variable assignment dictionary.
    braid_states
        BraidStates object.

    Returns
    -------
    tuple
        (singles, multiples) lists of inequality expressions.
    """
    singles = []
    multiples = []

    for inequality in relations:
        if not isinstance(inequality, (Leq, Less)):
            continue

        if inequality.first == ZERO_STATE:
            a = 0
        elif inequality.first == NEG_ONE_STATE:
            a = -1
        elif isinstance(inequality.first, tuple):
            a = assignment[braid_states.get_state(inequality.first)]
        else:
            continue

        if inequality.second == ZERO_STATE:
            b = 0
        elif inequality.second == NEG_ONE_STATE:
            b = -1
        elif isinstance(inequality.second, tuple):
            b = assignment[braid_states.get_state(inequality.second)]
        else:
            continue

        if isinstance(a, Symbol):
            a_dict = a.as_coefficients_dict()
        elif a == 0:
            a_dict = {}
        elif a == -1:
            a_dict = {one: -1}
        else:
            raise Exception(f'Expected variable "a" to be a Symbol, 0, or -1, but "a" was {a}!')

        if isinstance(b, Symbol):
            b_dict = b.as_coefficients_dict()
        elif b == 0:
            b_dict = {}
        elif b == -1:
            b_dict = {one: -1}
        else:
            raise Exception(f'Expected variable "b" to be a Symbol, 0, or -1, but "b" was {b}!')

        c_dict = {}
        for key in b_dict:
            if key in a_dict:
                c_dict[key] = b_dict[key] - a_dict[key]
            else:
                c_dict[key] = b_dict[key]
        for key in a_dict:
            if key not in b_dict:
                c_dict[key] = -a_dict[key]

        c_dict = {k: v for k, v in c_dict.items() if v != 0}

        if c_dict:
            expression = _expr_from_dict(c_dict)
            if len(set(c_dict.keys()) - {one}) == 1:
                singles.append(expression)
            else:
                multiples.append(expression)

    return list(set(singles)), list(set(multiples))


def process_assignment(assignment: Dict, braid_states, relations: List):
    """
    Process an assignment to extract criteria, multi-variable inequalities,
    and single-variable signs.

    Parameters
    ----------
    assignment
        Variable assignment dictionary.
    braid_states
        BraidStates object.
    relations
        List of constraint relations.

    Returns
    -------
    tuple
        (criteria, multi_var_inequalities, single_var_signs)
    """
    criteria = minimum_degree_symbolic(assignment, braid_states)
    singles, multi_var_inequalities = inequality_manager(relations, assignment, braid_states)
    single_var_signs = {}
    for entry in singles:
        dict_ = entry.as_coefficients_dict()
        single_var_signs[list(set(dict_.keys()) - {one})[0]] = list(dict_.values())[0]
    multi_var_inequalities = list(set(multi_var_inequalities))
    return criteria, multi_var_inequalities, single_var_signs
