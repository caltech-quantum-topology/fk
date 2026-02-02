"""
Sign assignment validation for FK computation.

This module provides functions for validating sign assignments
against constraint relations.
"""

from typing import Dict, List, Any, Optional

from ..braid.types import StateLiteral, ZERO_STATE, NUNITY_STATE
from ..constraints.relations import Leq, Less, Zero, Nunity, Alias, Conservation
from ..constraints.symbols import Symbol, one, zero, nunity
from .assignment import symbolic_variable_assignment


def violates_relation(assignment: Dict, relation: Any, verbose: bool = False) -> bool:
    """
    Check if an assignment violates a single relation.

    Parameters
    ----------
    assignment
        Variable assignment dictionary.
    relation
        A constraint relation.
    verbose
        Whether to print violation details.

    Returns
    -------
    bool
        True if the assignment violates the relation.
    """
    relation_type = type(relation)
    if relation_type == Leq:
        first = relation.first
        second = relation.second
        first_type = type(first)
        second_type = type(second)
        if first_type == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if second_type == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
        return (first_value > second_value)
    elif relation_type == Less:
        first = relation.first
        second = relation.second
        first_type = type(first)
        second_type = type(second)
        if first_type == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if second_type == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
        return (first_value >= second_value)
    elif relation_type == Alias:
        alias = relation.alias
        state = relation.state
        alias_value = assignment[alias]
        state_value = assignment[state]
        if (alias_value != state_value):
            if verbose:
                print('Alias')
        return (alias_value != state_value)
    elif relation_type == Conservation:
        sum_alias = relation.try_sum_alias()
        if sum_alias is not None:
            alias, sum_vars = sum_alias
            alias_value = assignment[alias]
            sum_value = sum([assignment[x] for x in sum_vars])
            if (alias_value != sum_value):
                if verbose:
                    print("Sum Alias")
            return (alias_value != sum_value)
        else:
            inputs = relation.inputs
            outputs = relation.outputs
            input_sum = sum([assignment[x] for x in inputs])
            output_sum = sum([assignment[x] for x in outputs])
            if (input_sum != output_sum):
                if verbose:
                    print("Conservation")
            return (input_sum != output_sum)
    elif relation_type == Zero:
        state = relation.state
        state_value = assignment[state]
        if (state_value != 0):
            if verbose:
                print("Zero")
        return (state_value != 0)
    elif relation_type == Nunity:
        state = relation.state
        state_value = assignment[state]
        if (state_value != -1):
            if verbose:
                print("Nunity")
        return (state_value != -1)
    return False


def violates_any_relation(assignment: Dict, relations: List, verbose: bool = False) -> bool:
    """
    Check if an assignment violates any relation in a list.

    Parameters
    ----------
    assignment
        Variable assignment dictionary.
    relations
        List of constraint relations.
    verbose
        Whether to print violation details.

    Returns
    -------
    bool
        True if the assignment violates any relation.
    """
    return any([violates_relation(assignment, r, verbose) for r in relations])


def _inequality_manager(relations: List, assignment: Dict, braid_states):
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
        if isinstance(inequality, Leq) or isinstance(inequality, Less):
            if inequality.first == ZERO_STATE:
                a = 0
            elif inequality.first == NUNITY_STATE:
                a = -1
            elif isinstance(inequality.first, tuple):
                a = assignment[braid_states.get_state(inequality.first)]
            if inequality.second == ZERO_STATE:
                b = 0
            elif inequality.second == NUNITY_STATE:
                b = -1
            elif isinstance(inequality.second, tuple):
                b = assignment[braid_states.get_state(inequality.second)]
            if isinstance(inequality, Leq):
                c = "<="
            elif isinstance(inequality, Less):
                c = "<"
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
            for key in b_dict.keys():
                if key in a_dict.keys():
                    c_dict[key] = b_dict[key] - a_dict[key]
                else:
                    c_dict[key] = b_dict[key]
            for key in a_dict.keys():
                if key not in b_dict.keys():
                    c_dict[key] = -a_dict[key]
            bad_keys = []
            for (key, value) in zip(c_dict.keys(), c_dict.values()):
                if value == 0:
                    bad_keys.append(key)
            for key in bad_keys:
                c_dict.pop(key)
            if c_dict:
                if len(set(c_dict.keys()) - set([one])) == 1:
                    expression = _expr_from_dict(c_dict)
                    singles.append(expression)
                else:
                    expression = _expr_from_dict(c_dict)
                    multiples.append(expression)

    return list(set(singles)), list(set(multiples))


def _expr_from_dict(dict_: Dict) -> Any:
    """Build a symbolic expression from a coefficients dictionary."""
    expression = 0
    for (key, value) in zip(dict_.keys(), dict_.values()):
        expression += value * key
    return expression


def _minimum_degree_symbolic(assignment: Dict, braid_states, verbose: bool = False) -> Dict:
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


def _process_assignment(assignment: Dict, braid_states, relations: List):
    """
    Process an assignment to extract criteria, multiples, and single signs.

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
        (criteria, multiples, singlesigns)
    """
    criteria = _minimum_degree_symbolic(assignment, braid_states)
    singles, multiples = _inequality_manager(relations, assignment, braid_states)
    singlesigns = {}
    for entry in singles:
        dict_ = entry.as_coefficients_dict()
        singlesigns[list(set(dict_.keys()) - set([one]))[0]] = list(dict_.values())[0]
    multiples = list(set(multiples))
    return criteria, multiples, singlesigns


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
    # Import here to avoid circular imports
    from ...solver.ilp import integral_bounded

    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = _process_assignment(assignment, braid_states, relations)
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


def czech_sign_assignment(degree: int, relations: List, braid_states, verbose: bool = False) -> Optional[Dict]:
    """
    Alternative sign assignment check with degree constraints transformed.

    Parameters
    ----------
    degree
        Degree bound to check.
    relations
        List of reduced constraint relations.
    braid_states
        BraidStates object with sign assignment.
    verbose
        Whether to print progress.

    Returns
    -------
    dict or None
        Dictionary with criteria, multiples, single_signs, and assignment if valid,
        None otherwise.
    """
    from ...solver.ilp import integral_bounded

    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = _process_assignment(assignment, braid_states, relations)

    if verbose:
        from .assignment import _expression_minimum
        for value in criteria.values():
            print(_expression_minimum(value, singlesigns))

    for (key, value) in criteria.items():
        criteria[key] = degree - value
    if not integral_bounded(multiples + list(criteria.values()), singlesigns):
        return None
    return {
        "criteria": criteria,
        "multiples": multiples,
        "single_signs": singlesigns,
        "assignment": assignment
    }
