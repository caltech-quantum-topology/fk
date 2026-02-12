"""
Variable assignment functions for FK computation.

This module provides functions for computing symbolic variable assignments
from constraint systems.
"""

from typing import Dict, List

from ..constraints.relations import Conservation, Alias, Zero, Nunity
from ..constraints.symbols import Symbol, symbols, solve
from ..constraints.reduction import free_variables, all_variables


def equivalence_assignment(assignment: Dict, braid_states) -> Dict:
    """
    Extend assignment to all equivalent state locations.

    Parameters
    ----------
    assignment
        Current variable assignment dictionary.
    braid_states
        BraidStates object containing state equivalence classes.

    Returns
    -------
    dict
        Extended assignment including all equivalent states.
    """
    update = {}
    for (key, value) in assignment.items():
        for state in braid_states.state_equivalence_classes[braid_states.get_state(key)]:
            update[state] = value
    assignment.update(update)
    return assignment


def find_expressions(relations: List, assignment: Dict, braid_states) -> List:
    """
    Find conservation expressions that can be used to solve for variables.

    Parameters
    ----------
    relations
        List of constraint relations.
    assignment
        Current variable assignment.
    braid_states
        BraidStates object.

    Returns
    -------
    list
        List of symbolic expressions from conservation constraints.
    """
    expressions = []
    assignment = equivalence_assignment(assignment, braid_states)
    keys = set(assignment.keys())
    for relation in relations:
        if not isinstance(relation, Conservation):
            continue
        sum_alias = relation.try_sum_alias()
        if sum_alias is not None:
            alias, sum_vars = sum_alias
            if not sum_vars:
                continue
            if not ({alias} | set(sum_vars) <= keys):
                continue
            expression = assignment[alias] - assignment[sum_vars[0]] - assignment[sum_vars[1]]
        else:
            if not relation.inputs:
                continue
            if not (set(relation.inputs + relation.outputs) <= keys):
                continue
            expression = assignment[relation.inputs[0]] + assignment[relation.inputs[1]] - assignment[relation.outputs[0]] - assignment[relation.outputs[1]]
        if expression != 0:
            expressions.append(expression)
    return expressions


def minimal_free(expressions: List, new: Dict = None, verbose: bool = False) -> Dict:
    """
    Find minimal set of free variables by eliminating through expressions.

    Parameters
    ----------
    expressions
        List of symbolic expressions to use for elimination.
    new
        Dictionary of substitutions found so far.
    verbose
        Whether to print progress.

    Returns
    -------
    dict
        Dictionary mapping symbols to their expressions in terms of remaining free variables.
    """
    if new is None:
        new = {}
    while expressions:
        considering = expressions.pop()
        syms = list(considering.free_symbols())
        if syms:
            update = {syms[0]: solve(considering, syms[0])[0]}
            for j in range(len(expressions)):
                expressions[j] = expressions[j].subs(update)
            for key in list(new.keys()):
                new[key] = new[key].subs(update)
            new.update(update)
    return new


def extend_variable_assignment(reduced_relations: List, partial_assignment: Dict, braid_states, verbose: bool = False) -> Dict:
    """
    Extend a partial variable assignment to all variables.

    Uses reduced relations to propagate assignments from assigned variables
    to unassigned ones through zeros, nunities, aliases, and sum aliases.

    Parameters
    ----------
    reduced_relations
        List of reduced constraint relations.
    partial_assignment
        Initial partial assignment of variables to symbols.
    braid_states
        BraidStates object.
    verbose
        Whether to print progress.

    Returns
    -------
    dict
        Complete assignment of all variables.
    """
    all_vars = all_variables(reduced_relations)

    if verbose:
        print("EXTENDING VARIABLE ASSIGNMENT")
        print(f"\tpartial assignment: {partial_assignment}")

    while True:
        assigned = partial_assignment.keys()
        unassigned = [x for x in all_vars if x not in assigned]

        if verbose:
            print(f"\tunassigned variables: {unassigned}")

        if not unassigned:
            return partial_assignment

        next_assignment = None
        for relation in reduced_relations:
            relation_type = type(relation)
            if relation_type == Zero:
                state = relation.state
                if verbose:
                    print(f"\tconsidering zero {state} := 0")
                if state not in assigned:
                    next_assignment = (state, 0)
                    break
            elif relation_type == Nunity:
                state = relation.state
                if verbose:
                    print(f"\tconsidering nunity {state} := -1")
                if state not in assigned:
                    next_assignment = (state, -1)
                    break
            elif relation_type == Alias:
                alias = relation.alias
                state = relation.state
                if verbose:
                    print(f"\tconsidering alias {alias} := {state}")
                if state in assigned and alias in unassigned:
                    next_assignment = (alias, partial_assignment[state])
                    break
            elif relation_type == Conservation:
                sum_alias = relation.try_sum_alias()
                if sum_alias is not None:
                    alias, sum_vars = sum_alias
                    if verbose:
                        print(f"\tconsidering sum alias {alias} := {sum_vars}")
                    if alias in unassigned and all(x in assigned for x in sum_vars):
                        if verbose:
                            print(f"restoring sumalias {alias} {sum_vars}")
                        next_assignment = (alias, sum(partial_assignment[x] for x in sum_vars))
                        break

        if next_assignment is not None:
            partial_assignment[next_assignment[0]] = next_assignment[1]
        else:
            raise Exception(f"could not find next assignment, even though not all locations have been assigned.\n\tunassigned variables: {unassigned}\n\tassignments: {partial_assignment}")


def symbolic_variable_assignment(relations: List, braid_states) -> Dict:
    """
    Create a symbolic variable assignment for a constraint system.

    Assigns fresh symbols to free variables, then extends the assignment
    to all variables using the constraint relations.

    Parameters
    ----------
    relations
        List of reduced constraint relations.
    braid_states
        BraidStates object.

    Returns
    -------
    dict
        Complete symbolic assignment of all state variables.
    """
    vars_list = free_variables(relations)
    assignment = dict(zip(vars_list, symbols(len(vars_list))))

    assignment = extend_variable_assignment(relations, assignment, braid_states)
    expressions = list(set(find_expressions(relations, assignment, braid_states)))
    minimalizer = minimal_free(expressions, {})
    for key, value in zip(list(assignment.keys()), list(assignment.values())):
        if value != 0 and value != -1:
            assignment[key] = value.subs(minimalizer)

    return assignment
