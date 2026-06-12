from utils import sort_any
import numpy as np

from symabelgroup import (
    Symbol,
    symbols,
    one,
    zero,
    expr_from_dict
)
from bounded_check import integral_bounded
from braidstates import BraidStates
from relations import (
    free_variables,
    print_inequalities,
    print_conservations,
    full_reduce,
    extend_variable_assignment,
    find_expressions,
    minimal_free,
    Leq,
    Less
)
from stateliteral import NUNITY_STATE, ZERO_STATE

'''
TO-DO: Change ILP's saved to save degree offsets, so executable can be run with different degrees without regenerating ILP
'''

def symbolic_variable_assignment(relations, braid_states):
    """
    Generate symbolic variable assignment for the constraint system.

    Args:
        relations: List of reduced relation objects
        braid_states: BraidStates object providing context

    Returns:
        Dict mapping state variables to symbolic expressions

    Process:
    1. Identify free variables (not constrained by relations)
    2. Assign symbolic variables to free variables
    3. Extend to complete assignment using constraint propagation
    4. Find constraint expressions from conservation relations
    5. Minimize the symbolic system by eliminating variables
    6. Apply minimization to get final symbolic assignment

    This creates a symbolic representation of all state variables
    in terms of a minimal set of free parameters.
    """
    # Step 1: Get free variables and assign symbols
    vars = free_variables(relations)
    assignment = dict(zip(vars, symbols(len(vars))))

    # Step 2: Extend to complete assignment using relations
    assignment = extend_variable_assignment(relations, assignment, braid_states)

    # Step 3: Extract constraint expressions from conservation relations
    expressions = list(set(find_expressions(relations, assignment, braid_states)))

    # Step 4: Minimize system by eliminating variables
    minimalizer = minimal_free(expressions, {})

    # Step 5: Apply minimization to get final symbolic forms
    for key, value in list(assignment.items()):
        if value != 0 and value != -1:
            assignment[key] = value.subs(minimalizer)

    return assignment

def minimum_degree_symbolic(assignment, braid_states, weight):
    conditions = {val:zero for val in range(braid_states.n_components)}
    for index in range(0, braid_states.n_strands): # starting at 0 seems to be correct for links with >= 2 components
        conditions[braid_states.closed_strand_components[index]] -= 1/2
    for index in range(braid_states.n_components):
        conditions[index] -= 1/2
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
    if weight != -1:
        conditions.update({-1:total_weight(assignment, braid_states)})
    return conditions

def inequality_manager(relations, assignment, braid_states):

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
                raise Exception(f'Expected variable "a" to be a Symbol, 0, or -1 (the latter two being integers), but "a" was {a}!')
            if isinstance(b, Symbol):
                b_dict = b.as_coefficients_dict()
            elif b == 0:
                b_dict = {}
            elif b == -1:
                b_dict = {one: -1}
            else:
                raise Exception(f'Expected variable "b" to be a Symbol, 0, or -1 (the latter two being integers), but "b" was {b}!')
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
                    expression = expr_from_dict(c_dict)
                    singles.append(expression)
                else:
                    expression = expr_from_dict(c_dict)
                    multiples.append(expression)

    return list(set(singles)), list(set(multiples))

def process_assignment(assignment, braid_states, relations, weight):
    criteria = minimum_degree_symbolic(assignment, braid_states, weight)
    singles, multiples = inequality_manager(relations, assignment, braid_states)
    singlesigns = {}
    for entry in singles:
        dict_ = entry.as_coefficients_dict()
        singlesigns[list(set(dict_.keys()) - set([one]))[0]] = list(dict_.values())[0]
    multiples = list(set(multiples))
    return criteria, multiples, singlesigns

def expression_minimum(bounding_expression, singlesigns):
    summand = 0
    for key, value in bounding_expression.as_coefficients_dict().items():
        if key == one:
            summand += value
        else:
            if singlesigns[key] == -1:
                summand -= value
    return summand

def check_sign_assignment(degree, relations, braid_states, weight):
    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations, weight)
    for (key, value) in criteria.items():
         if key >= 0:
             multiples.append(degree - value)
         else:
             multiples.append(weight - value)
    if not integral_bounded(multiples, singlesigns):
        return None
    if not integral_bounded(multiples, singlesigns):
        return None
    return {
            "criteria" : criteria,          # inequalities from degree bounding
            "multiples" : multiples,
            "single_signs" : singlesigns,
            "assignment" : assignment,       # symbolic assignment
        }

def czech_sign_assignment(degree, relations, braid_states, weight):
    assignment = symbolic_variable_assignment(relations, braid_states)
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations, weight)
    for (key, value) in criteria.items():
         if key >= 0:
             criteria[key] = degree - value
         else:
             criteria[key] = weight - value
    if not integral_bounded(multiples + list(criteria.values()), singlesigns):
        return None
    return {
            "criteria" : criteria,          # inequalities from degree bounding
            "multiples" : multiples,
            "single_signs" : singlesigns,
            "assignment" : assignment       # symbolic assignment
        }

def unreduced_variable_assignment(braid_states):
    return dict(zip(braid_states.state_locations, symbols(len(braid_states.state_locations))))

def print_unreduced(relations, braid_states, weight):
    assignment = unreduced_variable_assignment(braid_states)
    print("\nBraid:\n")
    print(braid_states.braid)
    print("\nInversion Data:\n")
    print(braid_states.sign_assignment)
    print("\nSymbolic Assignment:\n")
    print(assignment)
    print("\nInequalities:\n")
    print_inequalities(relations, assignment)
    print("\nConservation Relations:\n")
    print_conservations(relations, assignment)
    criteria, _, singlesigns = process_assignment(assignment, braid_states, relations, weight)
    print("\nMinimal x-Degree:\n")
    print(criteria[0].as_coefficients_dict())
    print("\nSigns of Each Variable:\n")
    print(singlesigns)

def total_weight(assignment, braid_states):
    acc = 0
    for index in range(braid_states.n_strands):
        if (braid_states.sign_assignment[(index, 0)] > 0):
            acc += assignment[(index, 0)]
        else:
            acc -= assignment[(index, 0)]
    return acc

def save(braid_states : BraidStates, degree, save_to, weight=-1):
    """
    Generate Integer Linear Program (ILP) file for Gukov-Manolescu invariant computation.

    Args:
        braid: List of braid generators representing the knot/link
        degree: Target degree for the invariant computation
        write_to: Output filename for the ILP file
        inversion_data: Sign assignment data (required for non-homogeneous braids)

    Returns:
        None if computation is infeasible, otherwise writes ILP to file

    Process:
    1. Create BraidStates object and generate relations
    2. Reduce relations to minimal form
    3. Generate symbolic variable assignment
    4. Convert to matrix form suitable for C++ ILP solver
    5. Write CSV format file with all constraint matrices

    File format:
    - Line 1: degree, n_components, writhe
    - Line 2: crossing data (generator, R-matrix type)
    - Line 3: component assignments for strands
    - Line 4: component assignments for crossings
    - Matrix 1: Degree criteria constraints (0 ≤ expressions)
    - Matrix 2: Inequality constraints from relations
    - Matrix 3: Variable assignment expressions
    """
    # Setup braid state analysis
    all_relations = braid_states.get_state_relations()
    relations = full_reduce(all_relations)

    # Check feasibility and generate constraint system
    check = czech_sign_assignment(degree, relations, braid_states, weight)
    if check is None:
        return None  # Infeasible

    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']

    criteria = list(criteria.values())
    n_criteria = len(criteria)
    n_multiples = len(multiples)

    # Determine matrix dimensions by finding maximum coefficient array size
    criteria_sizes = [criterion.var.size for criterion in criteria]
    inequality_sizes = [multiple.var.size for multiple in multiples]
    segment_sizes = [segment.var.size for segment in assignment.values() if not isinstance(segment, int)]
    for_fill = max(criteria_sizes + inequality_sizes + segment_sizes)

    # Build constraint matrices with zero-padding
    criteria_tableau = []
    for index in range(n_criteria):
        if criteria[index].is_constant():
            if criteria[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {criteria[index].constant()}")
        else:
            criteria_tableau.append(np.concat((criteria[index].var, np.zeros(for_fill - criteria_sizes[index]))))

    inequality_tableau = []
    for index in range(n_multiples):
        if multiples[index].is_constant():
            if multiples[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {multiples[index].constant()}")
        else:
            inequality_tableau.append(np.concat((multiples[index].var, np.zeros(for_fill - inequality_sizes[index]))))

    # Build assignment matrix
    assignment_tableau = []
    keys = sorted(assignment.keys())
    for key in keys:
        value = assignment[key]
        if isinstance(value, int):
            assignment_tableau.append(np.concat(([value], np.zeros(for_fill - 1))))
        else:
            assignment_tableau.append(np.concat((value.var, np.zeros(for_fill - len(value.var)))))

    # Convert to numpy arrays
    criteria_tableau = np.array(criteria_tableau)
    inequality_tableau = np.array(inequality_tableau)
    assignment_tableau = np.array(assignment_tableau)

    # Apply sign transformations for negative variables
    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                # Transform negative variable: replace x with -x + offset
                if criteria_tableau.shape[0] > 0:
                    criteria_tableau[:, index] *= -1
                    criteria_tableau[:, 0] += criteria_tableau[:, index]
                if inequality_tableau.shape[0] > 0:
                    inequality_tableau[:, index] *= -1
                    inequality_tableau[:, 0] += inequality_tableau[:, index]
                if assignment_tableau.shape[0] > 0:
                    assignment_tableau[:, index] *= -1
                    assignment_tableau[:, 0] += assignment_tableau[:, index]

    # Select relevant variable columns
    vars = sort_any([0] + [x.index() for x in list(singlesigns.keys())])
    if criteria_tableau.shape[0] > 0:
        criteria_tableau = criteria_tableau[:, vars]
    if inequality_tableau.shape[0] > 0:
        inequality_tableau = inequality_tableau[:, vars]
    if assignment_tableau.shape[0] > 0:
        assignment_tableau = assignment_tableau[:, vars]

    # Write ILP file in CSV format
    with open(save_to, 'w') as f:
        print(save_to, "degree", degree)
        # Header: basic parameters
        f.write(str(degree) + ",\n")
        f.write(str(braid_states.n_components) + ",\n")
        f.write(str(braid_states.writhe) + ",\n")

        # Crossing data: generator and R-matrix type
        for c in range(braid_states.n_crossings):
            f.write(str(abs(braid_states.braid[c])) + ",")
            f.write(str(braid_states.r_matrices[c])[1] + ",")  # Extract R-matrix number
        f.write("\n")

        # Component assignments
        for c in range(1, braid_states.n_strands):
            f.write(str(braid_states.closed_strand_components[c]) + ",")
        f.write("\n")
        for c in range(braid_states.n_crossings):
            f.write(str(braid_states.top_crossing_components[c]) + ",")
            f.write(str(braid_states.bottom_crossing_components[c]) + ",")
        f.write("\n")

        # Constraint matrices
        for row in criteria_tableau:
            for e in row:
                f.write(str(e) + ",")
            f.write("\n")
        f.write("/\n")

        for row in inequality_tableau:
            for e in row:
                f.write(str(int(e)) + ",")
            f.write("\n")
        f.write("/\n")

        for row in assignment_tableau:
            for e in row:
                f.write(str(int(e)) + ",")
            f.write("\n")
