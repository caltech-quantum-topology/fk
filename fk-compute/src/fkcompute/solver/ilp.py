"""
ILP formulation and Gurobi integration for FK computation.

This module provides functions for formulating and solving integer linear
programs used in the FK computation pipeline.
"""

from typing import Dict, List, Optional, Any

import numpy as np
import gurobipy as gp
from gurobipy import GRB

from ..domain.constraints.symbols import Symbol, one, zero, nunity
from ..domain.constraints.reduction import full_reduce
from ..domain.solver.assignment import symbolic_variable_assignment
from ..domain.solver.symbolic_constraints import (
    minimum_degree_symbolic,
    inequality_manager,
    process_assignment,
)


def _sort_any(xs):
    """Sort any list by string representation."""
    return list(sorted(xs, key=lambda x: str(x)))


# Initialize Gurobi environment
_env = gp.Env(empty=True)
_env.setParam('OutputFlag', 0)
_env.start()


def integral_bounded(multiples: List, singlesigns: Dict) -> bool:
    """
    Check if the constraint system has a bounded integer solution.

    Uses Gurobi to check if there exists a bounded integer solution
    satisfying all the inequality constraints.

    Parameters
    ----------
    multiples
        List of symbolic inequality expressions (0 <= expr).
    singlesigns
        Dictionary mapping symbols to their signs (+1 or -1).

    Returns
    -------
    bool
        True if the system has a bounded integer solution.
    """
    n_multiples = len(multiples)
    sizes = [multiple.var.size for multiple in multiples]
    for_fill = max(sizes) if sizes else 1
    tableau = []
    for index in range(n_multiples):
        if multiples[index].is_constant() and multiples[index].constant() == -1:
            return False
        tableau.append(np.concatenate((multiples[index].var, np.zeros(for_fill - sizes[index]))))
    tableau = np.array(tableau) if tableau else np.array([]).reshape(0, for_fill)
    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                tableau[:, index] *= -1
                tableau[:, 0] += tableau[:, index]
    model = gp.Model(env=_env)
    model.setParam(gp.GRB.Param.PoolSearchMode, 1)
    x = model.addVars(singlesigns.keys(), vtype=GRB.INTEGER)
    for index in range(n_multiples):
        value = tableau[index][0]
        for index_ in range(1, for_fill):
            if tableau[index][index_] != 0:
                value += tableau[index][index_] * x[Symbol(index_)]
        try:
            model.addConstr(0 <= value)
        except:
            pass
    for key in singlesigns.keys():
        model.setObjective(x[key], sense=GRB.MAXIMIZE)
        model.optimize()
        if model.Status != 2:
            return False
    if model.SolCount == 0:
        return False
    return True


def _czech_sign_assignment(degree: int, relations: List, braid_states, verbose: bool = False) -> Optional[Dict]:
    """Check sign assignment validity for ILP generation."""
    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations)

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


def ilp(degree: int, relations: List, braid_states, write_to: Optional[str] = None, verbose: bool = False) -> Optional[str]:
    """
    Generate ILP formulation for FK computation.

    Parameters
    ----------
    degree
        Degree bound for the computation.
    relations
        List of reduced constraint relations.
    braid_states
        BraidStates object with sign assignment.
    write_to
        Optional file path to write the ILP data.
    verbose
        Whether to print progress information.

    Returns
    -------
    str or None
        ILP data as a string, or None if no valid assignment exists.
    """
    check = _czech_sign_assignment(degree, relations, braid_states)
    if check is None:
        return None

    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']

    criteria = list(criteria.values())

    n_criteria = len(criteria)
    n_multiples = len(multiples)

    criteria_sizes = [criterion.var.size for criterion in criteria]
    inequality_sizes = [multiple.var.size for multiple in multiples]
    segment_sizes = [segment.var.size for segment in assignment.values() if not isinstance(segment, int)]
    for_fill = max(criteria_sizes + inequality_sizes + segment_sizes) if (criteria_sizes + inequality_sizes + segment_sizes) else 1

    criteria_tableau = []
    for index in range(n_criteria):
        if criteria[index].is_constant():
            if criteria[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {criteria[index].constant()}")
        else:
            criteria_tableau.append(np.concatenate((criteria[index].var, np.zeros(for_fill - criteria_sizes[index]))))

    inequality_tableau = []
    for index in range(n_multiples):
        if multiples[index].is_constant():
            if multiples[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {multiples[index].constant()}")
        else:
            inequality_tableau.append(np.concatenate((multiples[index].var, np.zeros(for_fill - inequality_sizes[index]))))

    assignment_tableau = []
    keys = sorted(assignment.keys())

    for key in keys:
        value = assignment[key]
        if isinstance(value, int):
            assignment_tableau.append(np.concatenate(([value], np.zeros(for_fill - 1))))
        else:
            assignment_tableau.append(np.concatenate((value.var, np.zeros(for_fill - len(value.var)))))

    criteria_tableau = np.array(criteria_tableau) if criteria_tableau else np.array([]).reshape(0, for_fill)
    inequality_tableau = np.array(inequality_tableau) if inequality_tableau else np.array([]).reshape(0, for_fill)
    assignment_tableau = np.array(assignment_tableau) if assignment_tableau else np.array([]).reshape(0, for_fill)

    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                if criteria_tableau.shape[0] > 0:
                    criteria_tableau[:, index] *= -1
                    criteria_tableau[:, 0] += criteria_tableau[:, index]
                if inequality_tableau.shape[0] > 0:
                    inequality_tableau[:, index] *= -1
                    inequality_tableau[:, 0] += inequality_tableau[:, index]
                if assignment_tableau.shape[0] > 0:
                    assignment_tableau[:, index] *= -1
                    assignment_tableau[:, 0] += assignment_tableau[:, index]

    vars_list = _sort_any([0] + [x.index() for x in list(singlesigns.keys())])

    if criteria_tableau.shape[0] > 0:
        criteria_tableau = criteria_tableau[:, vars_list]
    if inequality_tableau.shape[0] > 0:
        inequality_tableau = inequality_tableau[:, vars_list]
    if assignment_tableau.shape[0] > 0:
        assignment_tableau = assignment_tableau[:, vars_list]

    lines = []

    # header
    lines.append(f"{degree},")
    lines.append(f"{braid_states.n_components},")
    lines.append(f"{braid_states.writhe},")

    # braid data
    braid_line = []
    for c in range(braid_states.n_crossings):
        braid_line.append(str(abs(braid_states.braid[c])))
        braid_line.append(str(braid_states.r_matrices[c])[1])
    lines.append(",".join(braid_line) + ",")

    # closed strands
    closed_line = [str(braid_states.closed_strand_components[c])
                   for c in range(1, braid_states.n_strands)]
    lines.append(",".join(closed_line) + ",")

    # top & bottom crossing components
    cross_line = []
    for c in range(braid_states.n_crossings):
        cross_line.append(str(braid_states.top_crossing_components[c]))
        cross_line.append(str(braid_states.bottom_crossing_components[c]))
    lines.append(",".join(cross_line) + ",")

    # criteria tableau
    for row in criteria_tableau:
        lines.append(",".join(str(e) for e in row) + ",")
    lines.append("/")  # separator

    # inequality tableau
    for row in inequality_tableau:
        lines.append(",".join(str(int(e)) for e in row) + ",")
    lines.append("/")  # separator

    # assignment tableau
    for row in assignment_tableau:
        lines.append(",".join(str(int(e)) for e in row) + ",")

    # Join everything with newlines
    output = "\n".join(lines)

    # Write once
    if write_to:
        with open(write_to, "w") as f:
            f.write(output)
    return output


def print_symbolic_relations(degree: int, relations: List, braid_states, write_to: Optional[str] = None, verbose: bool = False) -> Optional[Dict]:
    """
    Print the reduced relations in human-readable format.

    Parameters
    ----------
    degree
        Degree bound for the computation.
    relations
        List of constraint relations.
    braid_states
        BraidStates object.
    write_to
        Optional file path to write output.
    verbose
        Whether to print additional information.

    Returns
    -------
    dict or None
        Dictionary with criteria, multiples, single_signs, and assignment if valid.
    """
    check = _czech_sign_assignment(degree, relations, braid_states)
    if check is None:
        print("No valid sign assignment exists for this degree!")
        return None

    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']

    print(f"=== SYMBOLIC VARIABLES RELATIONS AT DEGREE {degree} ===\n")

    print(f"Braid: {braid_states.braid}")
    print(f"Number of components: {braid_states.n_components}")
    print(f"Writhe: {braid_states.writhe}")
    if hasattr(braid_states, 'strand_signs') and braid_states.strand_signs:
        print(f"Inversion data: {braid_states.strand_signs}")
    print()

    print("=== VARIABLE ASSIGNMENTS ===")
    for state, expr in sorted(assignment.items()):
        print(f"{state} = {expr}")
    print()

    print("=== VARIABLE SIGNS ===")
    for var, sign in sorted(singlesigns.items(), key=lambda x: x[0].index() if hasattr(x[0], 'index') else 0):
        sign_str = "+" if sign == 1 else "-"
        print(f"{var}: {sign_str}")
    print()

    print("=== DEGREE CONSTRAINTS ===")
    for component, constraint in sorted(criteria.items()):
        if constraint.is_constant():
            print(f"Component {component}: 0 <= {constraint.constant()}")
        else:
            print(f"Component {component}: 0 <= {constraint}")
    print()

    print("=== RELATION INEQUALITIES ===")
    if multiples:
        for i, ineq in enumerate(multiples):
            if ineq.is_constant():
                print(f"Inequality {i+1}: 0 <= {ineq.constant()}")
            else:
                print(f"Inequality {i+1}: 0 <= {ineq}")
    else:
        print("No inequality constraints from relations")
    print()

    print("=== ORIGINAL RELATIONS ===")
    for i, relation in enumerate(relations):
        print(f"{i+1}. {relation}")
    print()

    free_vars = [var for var, expr in assignment.items()
                 if isinstance(expr, Symbol) and not expr.is_constant()]
    print("=== FREE VARIABLES ===")
    if free_vars:
        print(f"Free variables: {free_vars}")
    else:
        print("All variables are determined by relations")
    print()

    return {
        "criteria": criteria,
        "multiples": multiples,
        "single_signs": singlesigns,
        "assignment": assignment
    }
