"""
Constraint reduction and propagation for FK computation.

This module provides functions for reducing and simplifying constraint systems
through various propagation rules.
"""

from typing import List, Dict, Any

from ..braid.types import ZERO_STATE, NUNITY_STATE
from .relations import Leq, Less, Zero, Nunity, Alias, Conservation


def _sort_any(xs):
    """Sort any list by string representation."""
    return list(sorted(xs, key=lambda x: str(x)))


def get_zeros(relations: List) -> List:
    """Get all states constrained to be zero."""
    return [r.state for r in relations if type(r) == Zero]


def get_nunities(relations: List) -> List:
    """Get all states constrained to be negative unity."""
    return [r.state for r in relations if type(r) == Nunity]


def get_aliases(relations: List) -> Dict:
    """Get all alias relations as a dictionary."""
    return {r.alias: r.state for r in relations if type(r) == Alias}


def get_sum_aliases(relations: List) -> List:
    """Get all sum aliases from conservation relations."""
    return [x[0] for x in [r.try_sum_alias() for r in relations if type(r) == Conservation] if x is not None]


def propagate_zero_aliases(relations: List, verbose: bool = False) -> List:
    """Propagate zero constraints through aliases."""
    zeros = get_zeros(relations)
    aliases = get_aliases(relations)

    new_relations = []
    for r in relations:
        if type(r) == Alias:
            if r.state in zeros:
                if verbose:
                    print('reduction:', r, Zero(r.state), '===>', Zero(r.alias))
                if r.alias not in zeros:
                    new_relations.append(Zero(r.alias))
            elif r.alias in zeros:
                if verbose:
                    print('reduction:', r, Zero(r.alias), '===>', Zero(r.state))
                if r.state not in zeros:
                    new_relations.append(Zero(r.state))
            else:
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations


def propagate_nunity_aliases(relations: List, verbose: bool = False) -> List:
    """Propagate negative unity constraints through aliases."""
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []
    for r in relations:
        if type(r) == Alias:
            if r.state in nunities:
                if r.alias not in nunities:
                    new_relations.append(Nunity(r.alias))
            elif r.alias in nunities:
                if verbose:
                    print('reduction:', r, Nunity(r.alias), '===>', Nunity(r.state))
                if r.state not in nunities:
                    new_relations.append(Nunity(r.state))
            else:
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations


def de_alias_inequalities(relations: List, verbose: bool = False) -> List:
    """Replace aliased variables in inequalities with their canonical forms."""
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []

    for r in relations:
        if type(r) == Leq:
            f = r.first
            s = r.second
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in nunities:
                f = NUNITY_STATE
            elif r.first in aliases:
                f = aliases[r.first]
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in nunities:
                s = NUNITY_STATE
            elif r.second in aliases:
                s = aliases[r.second]
            r = Leq(f, s)
        elif type(r) == Less:
            f = r.first
            s = r.second
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in aliases:
                f = aliases[r.first]
            elif r.first in nunities:
                f = NUNITY_STATE
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in aliases:
                s = aliases[r.second]
            elif r.second in nunities:
                s = NUNITY_STATE
            r = Less(f, s)

        new_relations.append(r)

    return new_relations


def symmetric_inequality(relations: List, verbose: bool = False) -> List:
    """Detect symmetric inequalities that imply equality."""
    new_relations = [r for r in relations]
    inequalities = [r for r in relations if type(r) == Leq]

    for r1 in inequalities:
        for r2 in inequalities:
            if r1.second == r2.first and r2.second == r1.first:
                if r1.first == ZERO_STATE:
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.second))
                    new_relations.append(Zero(r1.second))
                elif r1.first == NUNITY_STATE:
                    new_relations.append(Nunity(r1.second))
                elif r1.second == ZERO_STATE:
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.first))
                    new_relations.append(Zero(r1.first))
                elif r1.second == NUNITY_STATE:
                    new_relations.append(Nunity(r1.first))
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
        if type(r) == Conservation:
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
        if type(r) == Conservation:
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
        if type(r) == Conservation and len(r.inputs) == 1 and len(r.outputs) == 1:
            alias = Alias(r.inputs[0], r.outputs[0])
            new_relations.append(alias)
            if verbose:
                print('unary conservation is alias', r, '=====>', alias)
        else:
            new_relations.append(r)

    return new_relations


def is_vacuous(relation: Any) -> bool:
    """Check if a relation is vacuously true."""
    if type(relation) == Leq:
        if relation.first == relation.second:
            return True
    if type(relation) == Alias:
        if relation.alias == relation.state:
            return True
    if type(relation) == Zero and relation.state == ZERO_STATE:
        return True
    if type(relation) == Nunity and relation.state == NUNITY_STATE:
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
    """Apply all reduction rules once."""
    reduction_rules = [
        de_alias_inequalities,
        symmetric_inequality,
        propagate_zero_aliases,
        propagate_nunity_aliases,
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
            if v != ZERO_STATE and v != NUNITY_STATE:
                all_vars.append(v)
    all_vars = list(set(all_vars))
    return all_vars


def free_variables(relations: List) -> List:
    """Get all free (unbound) variables in relations."""
    all_vars = all_variables(relations)
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)
    sum_aliases = get_sum_aliases(relations)
    res = []
    for v in all_vars:
        if v not in zeros and v not in nunities and v not in aliases and v not in sum_aliases:
            res.append(v)
    return _sort_any(res)
