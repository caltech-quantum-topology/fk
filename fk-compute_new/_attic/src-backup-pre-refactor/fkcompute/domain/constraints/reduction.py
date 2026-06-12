"""
Constraint reduction and propagation for FK computation.

This module provides functions for reducing and simplifying constraint systems
through various propagation rules.
"""

from typing import List, Dict, Any, Type

from ..braid.types import ZERO_STATE, NEG_ONE_STATE
from .relations import _sort_any, Leq, Less, Zero, NegOne, Alias, Conservation


def get_zeros(relations: List) -> List:
    """Get all states constrained to be zero."""
    return [r.state for r in relations if isinstance(r, Zero)]


def get_neg_ones(relations: List) -> List:
    """Get all states constrained to be negative unity."""
    return [r.state for r in relations if isinstance(r, NegOne)]


def get_aliases(relations: List) -> Dict:
    """Get all alias relations as a dictionary."""
    return {r.alias: r.state for r in relations if isinstance(r, Alias)}


def get_sum_aliases(relations: List) -> List:
    """Get all sum aliases from conservation relations."""
    return [x[0] for x in [r.try_sum_alias() for r in relations if isinstance(r, Conservation)] if x is not None]


def _propagate_constant_aliases(relations: List, constant_type: Type, verbose: bool = False) -> List:
    """
    Propagate constant constraints (Zero or NegOne) through aliases.

    When an alias links two states and one is known to be a constant,
    the other must also be that constant.
    """
    if constant_type == Zero:
        constants = get_zeros(relations)
    else:
        constants = get_neg_ones(relations)

    new_relations = []
    for r in relations:
        if isinstance(r, Alias):
            if r.state in constants:
                if verbose:
                    print('reduction:', r, constant_type(r.state), '===>', constant_type(r.alias))
                if r.alias not in constants:
                    new_relations.append(constant_type(r.alias))
            elif r.alias in constants:
                if verbose:
                    print('reduction:', r, constant_type(r.alias), '===>', constant_type(r.state))
                if r.state not in constants:
                    new_relations.append(constant_type(r.state))
            else:
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations


def propagate_zero_aliases(relations: List, verbose: bool = False) -> List:
    """Propagate zero constraints through aliases."""
    return _propagate_constant_aliases(relations, Zero, verbose)


def propagate_neg_one_aliases(relations: List, verbose: bool = False) -> List:
    """Propagate negative one constraints through aliases."""
    return _propagate_constant_aliases(relations, NegOne, verbose)


def _resolve_variable(var, zeros, neg_ones, aliases):
    """Resolve a variable through zeros, neg_ones, and aliases."""
    if var in zeros:
        return ZERO_STATE
    elif var in neg_ones:
        return NEG_ONE_STATE
    elif var in aliases:
        return aliases[var]
    return var


def de_alias_inequalities(relations: List, verbose: bool = False) -> List:
    """Replace aliased variables in inequalities with their canonical forms."""
    zeros = get_zeros(relations)
    neg_ones = get_neg_ones(relations)
    aliases = get_aliases(relations)

    new_relations = []

    for r in relations:
        if isinstance(r, (Leq, Less)):
            f = _resolve_variable(r.first, zeros, neg_ones, aliases)
            s = _resolve_variable(r.second, zeros, neg_ones, aliases)
            r = type(r)(f, s)

        new_relations.append(r)

    return new_relations


def symmetric_inequality(relations: List, verbose: bool = False) -> List:
    """Detect symmetric inequalities that imply equality."""
    new_relations = [r for r in relations]
    inequalities = [r for r in relations if isinstance(r, Leq)]

    for r1 in inequalities:
        for r2 in inequalities:
            if r1.second == r2.first and r2.second == r1.first:
                if r1.first == ZERO_STATE:
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.second))
                    new_relations.append(Zero(r1.second))
                elif r1.first == NEG_ONE_STATE:
                    new_relations.append(NegOne(r1.second))
                elif r1.second == ZERO_STATE:
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.first))
                    new_relations.append(Zero(r1.first))
                elif r1.second == NEG_ONE_STATE:
                    new_relations.append(NegOne(r1.first))
                else:
                    if verbose:
                        print(r1, r2, '=====>', Alias(r1.first, r1.second))
                    new_relations.append(Alias(r1.first, r1.second))
    return new_relations


def conservation_zeros(relations: List, verbose: bool = False) -> List:
    """Remove zeros from conservation constraints."""
    zeros = get_zeros(relations)
    new_relations = []
    for r in relations:
        if isinstance(r, Conservation):
            inputs = [v for v in r.inputs if v != ZERO_STATE and v not in zeros]
            outputs = [v for v in r.outputs if v != ZERO_STATE and v not in zeros]
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('remove zeros:', r, '=====>', Conservation(inputs, outputs))
            new_relations.append(Conservation(inputs, outputs))
        else:
            new_relations.append(r)
    return new_relations


def conservation_alias(relations: List, verbose: bool = False) -> List:
    """Replace aliased variables in conservation constraints."""
    aliases = get_aliases(relations)
    new_relations = []
    for r in relations:
        if isinstance(r, Conservation):
            inputs = [aliases.get(v) or v for v in r.inputs]
            outputs = [aliases.get(v) or v for v in r.outputs]
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('conservation alias:', r, '=====>', Conservation(inputs, outputs))
            new_relations.append(Conservation(inputs, outputs))
        else:
            new_relations.append(r)
    return new_relations


def unary_conservation_is_alias(relations: List, verbose: bool = False) -> List:
    """Convert unary conservation constraints to aliases."""
    new_relations = []
    for r in relations:
        if isinstance(r, Conservation) and len(r.inputs) == 1 and len(r.outputs) == 1:
            alias = Alias(r.inputs[0], r.outputs[0])
            new_relations.append(alias)
            if verbose:
                print('unary conservation is alias', r, '=====>', alias)
        else:
            new_relations.append(r)

    return new_relations


def is_vacuous(relation: Any) -> bool:
    """Check if a relation is vacuously true."""
    if isinstance(relation, Leq):
        if relation.first == relation.second:
            return True
    if isinstance(relation, Alias):
        if relation.alias == relation.state:
            return True
    if isinstance(relation, Zero) and relation.state == ZERO_STATE:
        return True
    if isinstance(relation, NegOne) and relation.state == NEG_ONE_STATE:
        return True

    return False


def remove_vacuous(relations: List, verbose: bool = False) -> List:
    """Remove vacuously true relations."""
    new_relations = []
    for r in relations:
        if not is_vacuous(r):
            new_relations.append(r)
        elif verbose:
            print('removing vacuous:', r)
    return new_relations


def reduce_relations(relations: List, verbose: bool = False) -> List:
    """
    Apply all reduction rules once.

    Rule order matters: de-aliasing inequalities first allows symmetric_inequality
    to detect equalities. Propagation rules then simplify through the new equalities.
    Conservation rules run after propagation to benefit from resolved aliases/zeros.
    Vacuous removal runs last to clean up any trivially true results.
    """
    reduction_rules = [
        de_alias_inequalities,
        symmetric_inequality,
        propagate_zero_aliases,
        propagate_neg_one_aliases,
        conservation_alias,
        conservation_zeros,
        unary_conservation_is_alias,
        remove_vacuous
    ]
    for rule in reduction_rules:
        relations = rule(relations, verbose=verbose)

    relations = list(set(relations))
    return _sort_any(relations)


def full_reduce(relations: List, verbose: bool = False, max_depth: int = 20) -> List:
    """Fully reduce relations by repeatedly applying reduction rules until fixed point."""
    if max_depth == 0:
        raise Exception('max reduction recursion depth exceeded')

    new_relations = reduce_relations(relations, verbose=verbose)
    if set(relations) == set(new_relations):
        return new_relations
    else:
        return full_reduce(new_relations, verbose=verbose, max_depth=max_depth - 1)


def all_variables(relations: List) -> List:
    """Get all variables appearing in relations."""
    all_vars = []
    for r in relations:
        for v in r.variables():
            if v != ZERO_STATE and v != NEG_ONE_STATE:
                all_vars.append(v)
    all_vars = list(set(all_vars))
    return all_vars


def free_variables(relations: List) -> List:
    """Get all free (unbound) variables in relations."""
    all_vars = all_variables(relations)
    zeros = get_zeros(relations)
    neg_ones = get_neg_ones(relations)
    aliases = get_aliases(relations)
    sum_aliases = get_sum_aliases(relations)
    res = []
    for v in all_vars:
        if v not in zeros and v not in neg_ones and v not in aliases and v not in sum_aliases:
            res.append(v)
    return _sort_any(res)
