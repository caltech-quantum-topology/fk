"""
Relation management system for algebraic constraint solving in Gukov-Manolescu invariants.

This module defines various types of algebraic relations and provides algorithms for:
- Relation reduction and simplification
- Variable assignment and constraint propagation
- Integer Linear Programming (ILP) formulation
- Symbolic computation and constraint solving

Key relation types:
- Inequality relations: Leq (≤), Less (<)
- Equality relations: Zero (= 0), Nunity (= -1), Alias (variable equivalence)
- Conservation relations: sum preservation across braid crossings

The system supports automatic relation reduction, variable elimination, and
conversion to ILP format for the C++ solver component.

TO-DO: Fix the bad coding practice (Cmd + F "bad coding practice")
"""

from utils import sort_any
import numpy as np

from symabelgroup import solve, nunity, zero
from stateliteral import StateLiteral, NUNITY_STATE, ZERO_STATE

class Leq:
    """
    Less than or equal inequality relation: first ≤ second.
    
    Used to represent ordering constraints between state variables
    in the Gukov-Manolescu invariant computation.
    """
    
    def __init__(self, first, second):
        """Initialize inequality relation first ≤ second."""
        self.first = first
        self.second = second

    def variables(self):
        """Return list of variables involved in this relation."""
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} <= {self.second}'

    def __eq__(self, other):
        if isinstance(other, Leq):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))
    
class Less:
    """
    Strict less than inequality relation: first < second.
    
    Used for strict ordering constraints between state variables.
    """
    
    def __init__(self, first, second):
        """Initialize strict inequality relation first < second."""
        self.first = first
        self.second = second

    def variables(self):
        """Return list of variables involved in this relation."""
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} < {self.second}'

    def __eq__(self, other):
        if isinstance(other, Less):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))

class Zero:
    """
    Zero equality relation: state = 0.
    
    Represents state variables that must equal zero in the
    constraint system, often arising from boundary conditions.
    """
    
    def __init__(self, state):
        """Initialize zero relation state = 0."""
        self.state = state

    def variables(self):
        """Return list of variables involved in this relation."""
        return [self.state]

    def __repr__(self):
        return f'Zero {self.state} = [0]'

    def __eq__(self, other):
        if isinstance(other, Zero):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)
    
class Nunity:
    """
    Negative unity equality relation: state = -1.
    
    Represents state variables that must equal negative one,
    typically arising from specific boundary conditions in knot theory.
    """
    
    def __init__(self, state):
        """Initialize negative unity relation state = -1."""
        self.state = state

    def variables(self):
        """Return list of variables involved in this relation."""
        return [self.state]

    def __repr__(self):
        return f'Nunity {self.state} = [-1]'  # Fixed: was incorrectly "Zero"

    def __eq__(self, other):
        if isinstance(other, Nunity):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)

class Alias:
    """
    Variable alias relation: alias = state.
    
    Represents equality between two state variables, creating
    an equivalence relation for variable elimination.
    """
    
    def __init__(self, state, alias):
        """Initialize alias relation with canonical ordering."""
        state, alias = sort_any([state, alias])  # Ensure consistent ordering
        self.state = state
        self.alias = alias

    def variables(self):
        """Return list of variables involved in this relation."""
        return [self.state, self.alias]

    def __repr__(self):
        return f'Alias {self.alias} := {self.state}'

    def __eq__(self, other):
        if isinstance(other, Alias):
            return self.state == other.state and self.alias == other.alias
        return False

    def __hash__(self):
        return hash((self.state, self.alias))

class Conservation:
    """
    Conservation relation: sum(inputs) = sum(outputs).
    
    Represents conservation laws at braid crossings where the sum
    of input state values equals the sum of output state values.
    This encodes the fundamental constraint that crossing operations
    preserve the total "charge" or state content.
    """
    
    def __init__(self, inputs, outputs):
        """
        Initialize conservation relation with canonical ordering.
        
        Args:
            inputs: List of input state variables
            outputs: List of output state variables
        """
        # Canonicalize ordering for consistent comparison and hashing
        inputs = sort_any(inputs)
        outputs = sort_any(outputs)
        inputs, outputs = sort_any([inputs, outputs])
        self.inputs = inputs
        self.outputs = outputs

    def variables(self):
        """Return list of all variables involved in this relation."""
        return self.inputs + self.outputs

    def try_sum_alias(self):
        """
        Attempt to extract a sum alias from this conservation relation.
        
        Returns:
            Tuple (variable, sum_terms) if one side is a single variable
            and the other side is a sum, None otherwise.
            
        Patterns recognized:
        - a = b + c → (a, [b, c])
        - a + b = c → (c, [a, b])
        """
        # Case: both sides multi-term, no alias possible
        if len(self.inputs) != 1 and len(self.outputs) != 1:
            return None
        
        # Case: single input = multiple outputs → alias input to sum of outputs
        elif len(self.inputs) == 1 and type(self.inputs[0]) != StateLiteral:
            return self.inputs[0], self.outputs

        # Case: multiple inputs = single output → alias output to sum of inputs  
        elif len(self.outputs) == 1 and type(self.outputs[0]) != StateLiteral:
            return self.outputs[0], self.inputs
            
        return None

    def __repr__(self):
        input_sum = " + ".join(str(x) for x in self.inputs)
        output_sum = " + ".join(str(x) for x in self.outputs)
        return f'Conservation {input_sum} = {output_sum}'

    def __eq__(self, other):
        if isinstance(other, Conservation):
            return self.inputs == other.inputs and self.outputs == other.outputs
        return False

    def __hash__(self):
        return hash((tuple(self.inputs), tuple(self.outputs)))

def get_zeros(relations):
    """Extract all state variables that must equal zero."""
    return [r.state for r in relations if type(r) == Zero]

def get_nunities(relations):
    """Extract all state variables that must equal -1."""
    return [r.state for r in relations if type(r) == Nunity]

def get_aliases(relations):
    """Extract alias mapping: {alias_variable: canonical_variable}."""
    return {r.alias: r.state for r in relations if type(r) == Alias}

def get_sum_aliases(relations):
    """Extract variables that can be expressed as sums of other variables."""
    return [x[0] for x in [r.try_sum_alias() for r in relations if type(r) == Conservation] if x is not None]

def propagate_zero_aliases(relations, verbose=False):
    """
    Propagate zero values through alias relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with zero values propagated through aliases
        
    Logic:
    - If alias A := B and B = 0, then A = 0
    - If alias A := B and A = 0, then B = 0
    
    This ensures that aliased variables maintain consistent zero values.
    """
    zeros = get_zeros(relations)
    aliases = get_aliases(relations)
    new_relations = []
    
    for r in relations:
        if type(r) == Alias:
            if r.state in zeros:
                # Canonical variable is zero → alias must be zero
                if verbose:
                    print('reduction:', r, Zero(r.state), '===>', Zero(r.alias))
                if r.alias not in zeros:
                    new_relations.append(Zero(r.alias))
            elif r.alias in zeros:
                # Alias is zero → canonical variable must be zero
                if verbose:
                    print('reduction:', r, Zero(r.alias), '===>', Zero(r.state))
                if r.state not in zeros:
                    new_relations.append(Zero(r.state))
            else:
                # Neither is zero, keep the alias
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations

def propagate_nunity_aliases(relations, verbose=False):
    """
    Propagate negative unity (-1) values through alias relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with -1 values propagated through aliases
        
    Logic:
    - If alias A := B and B = -1, then A = -1
    - If alias A := B and A = -1, then B = -1
    
    This ensures that aliased variables maintain consistent -1 values.
    """
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []
    for r in relations:
        if type(r) == Alias:
            if r.state in nunities:
                # Canonical variable is -1 → alias must be -1
                if r.alias not in nunities:
                    new_relations.append(Nunity(r.alias))
            elif r.alias in nunities:
                # Alias is -1 → canonical variable must be -1
                if verbose:
                    print('reduction:', r, Nunity(r.alias), '===>', Nunity(r.state))
                if r.state not in nunities:
                    new_relations.append(Nunity(r.state))
            else:
                # Neither is -1, keep the alias
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations

def de_alias_inequalities(relations, verbose=False):
    """
    Replace aliased variables in inequality relations with their canonical forms.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with variables replaced by their canonical forms
        
    Replacement rules:
    - Variables with known values (0, -1) → corresponding state literals
    - Aliased variables → canonical variable from alias relation
    - Other variables → unchanged
    
    This simplifies inequalities by removing intermediate alias variables.
    """
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []
    
    for r in relations:
        if type(r) == Leq:
            # Process ≤ inequalities
            f = r.first
            s = r.second
            
            # Replace first variable
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in nunities:
                f = NUNITY_STATE
            elif r.first in aliases:
                f = aliases[r.first]
                
            # Replace second variable
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in nunities:
                s = NUNITY_STATE
            elif r.second in aliases:
                s = aliases[r.second]
                
            r = Leq(f, s)
            
        elif type(r) == Less:
            # Process < inequalities
            f = r.first
            s = r.second
            
            # Replace first variable
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in nunities:
                f = NUNITY_STATE
            elif r.first in aliases:
                f = aliases[r.first]
                
            # Replace second variable  
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in nunities:
                s = NUNITY_STATE
            elif r.second in aliases:
                s = aliases[r.second]
                
            r = Less(f, s)

        new_relations.append(r)
        
    return new_relations

def symmetric_inequality(relations, verbose=False):
    """
    Detect symmetric inequalities and convert them to equality relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with symmetric inequalities converted to equalities
        
    Logic:
    If we have both A ≤ B and B ≤ A, then A = B, which means:
    - If A = 0, then B = 0
    - If A = -1, then B = -1  
    - If B = 0, then A = 0
    - If B = -1, then A = -1
    - Otherwise, A and B are aliases (A := B)
    
    This is a key inference rule for constraint satisfaction.
    """
    new_relations = [r for r in relations]
    inequalities = [r for r in relations if type(r) == Leq]

    for r1 in inequalities:
        for r2 in inequalities:
            # Check for symmetric pair: A ≤ B and B ≤ A
            if r1.second == r2.first and r2.second == r1.first:
                if r1.first == ZERO_STATE:
                    # 0 ≤ B and B ≤ 0 → B = 0
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.second))
                    new_relations.append(Zero(r1.second))
                elif r1.first == NUNITY_STATE:
                    # -1 ≤ B and B ≤ -1 → B = -1
                    new_relations.append(Nunity(r1.second))
                elif r1.second == ZERO_STATE:
                    # A ≤ 0 and 0 ≤ A → A = 0
                    if verbose:
                        print(r1, r2, '=====>', Zero(r1.first))
                    new_relations.append(Zero(r1.first))
                elif r1.second == NUNITY_STATE:
                    # A ≤ -1 and -1 ≤ A → A = -1
                    new_relations.append(Nunity(r1.first))
                else:
                    # A ≤ B and B ≤ A → A = B (alias)
                    if verbose:
                        print(r1, r2, '=====>', Alias(r1.first, r1.second))
                    new_relations.append(Alias(r1.first, r1.second))
                    
    return new_relations

def conservation_zeros(relations, verbose=False):
    """
    Remove zero variables from conservation relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with zero variables removed from conservation equations
        
    Logic:
    Since zero variables contribute nothing to sums, they can be removed:
    A + 0 + B = C + 0 + D → A + B = C + D
    
    This simplifies conservation relations and may enable further reductions.
    """
    zeros = get_zeros(relations)
    new_relations = []
    
    for r in relations:
        if type(r) == Conservation:
            # Remove zeros and ZERO_STATE literals from both sides
            inputs = [v for v in r.inputs if v != ZERO_STATE and v not in zeros]
            outputs = [v for v in r.outputs if v != ZERO_STATE and v not in zeros]
            
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('remove zeros:', r, '=====>', Conservation(inputs, outputs))
                
            new_relations.append(Conservation(inputs, outputs))
        else:
            new_relations.append(r)
    return new_relations
            
def conservation_alias(relations, verbose=False):
    """
    Replace aliased variables in conservation relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with aliased variables replaced by canonical forms
        
    This simplifies conservation relations by using canonical variable names.
    """
    aliases = get_aliases(relations)
    new_relations = []
    
    for r in relations:
        if type(r) == Conservation:
            # Replace variables with their canonical forms (or keep if not aliased)
            inputs = [aliases.get(v, v) for v in r.inputs]  # More concise than aliases.get(v) or v
            outputs = [aliases.get(v, v) for v in r.outputs]
            
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('conservation alias:', r, '=====>', Conservation(inputs, outputs))
                
            new_relations.append(Conservation(inputs, outputs))
        else:
            new_relations.append(r)
    return new_relations

def unary_conservation_is_alias(relations, verbose=False):
    """
    Convert single-variable conservation relations to alias relations.
    
    Args:
        relations: List of relation objects
        verbose: If True, print reduction steps
        
    Returns:
        Updated relations with unary conservation converted to aliases
        
    Logic:
    If A = B (single variable on each side), this is equivalent to alias A := B.
    This conversion enables alias propagation algorithms to work on these relations.
    """
    new_relations = []
    
    for r in relations:
        if type(r) == Conservation and len(r.inputs) == 1 and len(r.outputs) == 1:
            # Single variable on each side → convert to alias
            alias = Alias(r.inputs[0], r.outputs[0])
            new_relations.append(alias)
            if verbose:
                print('unary conservation is alias', r, '=====>', alias)
        else:
            new_relations.append(r)

    return new_relations
                
def is_vacuous(relation):
    """
    Check if a relation is vacuously true (always satisfied).
    
    Args:
        relation: Relation object to check
        
    Returns:
        True if relation is vacuous, False otherwise
        
    Vacuous relations:
    - A ≤ A (variable ≤ itself)
    - A := A (variable aliased to itself)  
    - 0 = 0 (ZERO_STATE equals zero)
    - -1 = -1 (NUNITY_STATE equals -1)
    
    These relations provide no constraint information and can be removed.
    """
    if type(relation) == Leq:
        if relation.first == relation.second:
            return True  # A ≤ A
    elif type(relation) == Alias:
        if relation.alias == relation.state:
            return True  # A := A
    elif type(relation) == Zero and relation.state == ZERO_STATE:
        return True  # 0 = 0
    elif type(relation) == Nunity and relation.state == NUNITY_STATE:
        return True  # -1 = -1

    return False

def remove_vacuous(relations, verbose=False):
    """
    Remove vacuously true relations from the relation set.
    
    Args:
        relations: List of relation objects
        verbose: If True, print removed relations
        
    Returns:
        Filtered relations with vacuous ones removed
        
    This cleanup step removes relations that are always true and
    provide no constraint information.
    """
    new_relations = []
    for r in relations:
        if not is_vacuous(r):
            new_relations.append(r)
        elif verbose:
            print('removing vacuous:', r)
    return new_relations

def reduce_relations(relations, verbose=False):
    """
    Apply one round of relation reduction rules.
    
    Args:
        relations: List of relation objects to reduce
        verbose: If True, print reduction steps
        
    Returns:
        Reduced list of relations after applying all reduction rules
        
    Reduction rules applied in order:
    1. de_alias_inequalities: Replace aliased variables in inequalities
    2. symmetric_inequality: Detect symmetric inequalities → equality/alias
    3. propagate_zero_aliases: Propagate zero values through aliases  
    4. propagate_nunity_aliases: Propagate -1 values through aliases
    5. conservation_alias: Replace aliased variables in conservation relations
    6. conservation_zeros: Remove zero variables from conservation relations
    7. unary_conservation_is_alias: Convert single-variable conservation to alias
    8. remove_vacuous: Remove trivially satisfied relations
    """
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
    
    # Remove duplicates and canonicalize ordering
    relations = list(set(relations))
    return sort_any(relations)

def full_reduce(relations, verbose=False, max_depth=20):
    """
    Repeatedly apply relation reduction until no further changes occur.
    
    Args:
        relations: List of relation objects to reduce
        verbose: If True, print reduction steps
        max_depth: Maximum recursion depth to prevent infinite loops
        
    Returns:
        Fully reduced list of relations
        
    Raises:
        Exception: If max recursion depth is exceeded
        
    This implements a fixed-point algorithm that repeatedly applies
    reduction rules until the relation set stabilizes.
    """
    if max_depth == 0:
        raise Exception('max reduction recursion depth exceeded')
        
    new_relations = reduce_relations(relations, verbose=verbose)
    if set(relations) == set(new_relations):
        # Fixed point reached
        return new_relations
    else:
        # Continue reducing
        return full_reduce(new_relations, verbose=verbose, max_depth=max_depth-1)

def all_variables(relations):
    """
    Extract all state variables from a list of relations.
    
    Args:
        relations: List of relation objects
        
    Returns:
        List of all variables referenced in relations, excluding literals
    """
    all_variables = []
    for r in relations:
        for v in r.variables():
            if v != ZERO_STATE and v != NUNITY_STATE:
                all_variables.append(v)
    all_variables = list(set(all_variables))  # Remove duplicates
    return all_variables

def free_variables(relations):
    """
    Find variables that are not constrained to specific values or aliases.
    
    Args:
        relations: List of relation objects
        
    Returns:
        List of variables that remain free (not bound by constraints)
        
    A variable is considered bound if it:
    - Must equal zero (Zero relation)
    - Must equal -1 (Nunity relation)  
    - Is aliased to another variable (Alias relation)
    - Can be expressed as sum of other variables (Conservation relation)
    """
    all_vars = all_variables(relations)
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)
    sum_aliases = get_sum_aliases(relations)
    
    res = []
    for v in all_vars:
        if v not in zeros and v not in nunities and v not in aliases and v not in sum_aliases:
            res.append(v)
    return sort_any(res)

def equivalence_assignment (assignment, braid_states):
    update = {}
    for (key, value) in assignment.items():
        for state in braid_states.state_equivalence_classes[braid_states.get_state(key)]:
            update[state] = value
    assignment.update(update)
    return assignment

def find_expressions(relations, assignment, braid_states):
    expressions = []
    assignment = equivalence_assignment(assignment, braid_states)
    for relation in relations:
        if isinstance(relation, Conservation):
            sum_alias = relation.try_sum_alias()
            expression = None
            if sum_alias is not None:
                alias, sum_vars = sum_alias
                valid_list = [alias] + sum_vars
                if all([v in list(assignment.keys()) for v in valid_list]):
                    if len(sum_vars) > 0: # modifying: bad coding practice
                        expression = assignment[alias] - assignment[sum_vars[0]] - assignment[sum_vars[1]]
            else:
                if len(relation.inputs) > 0:# modifying: bad coding practice
                    valid_list = relation.inputs + relation.outputs
                    if all([v in list(assignment.keys()) for v in valid_list]):
                        expression = assignment[relation.inputs[0]] + assignment[relation.inputs[1]] - assignment[relation.outputs[0]] - assignment[relation.outputs[1]]
            if expression != None and expression != 0:
                expressions.append(expression)
    return expressions
        
def minimal_free(expressions, new={}, verbose=False):
    if not expressions:
        return new
    considering = expressions.pop()
    syms = list(considering.free_symbols()) 
    if syms:
        update = {syms[0]: solve(considering, syms[0])[0]}
        for j in range(len(expressions)):
            expressions[j] = expressions[j].subs(update)
        keys, values = list(new.keys()), list(new.values())
        for (key, value) in zip(keys, values):
            new[key] = value.subs(update)
        new.update(update)

    return minimal_free(expressions, new)

def extend_variable_assignment(reduced_relations, partial_assignment, braid_states, verbose=False):
    """
    Extend partial variable assignment to complete assignment using relations.
    
    Args:
        reduced_relations: List of reduced relation objects
        partial_assignment: Dict with current variable assignments
        braid_states: BraidStates object (for debugging context)
        verbose: If True, print assignment steps
        
    Returns:
        Complete variable assignment dict
        
    Raises:
        Exception: If assignment cannot be completed (inconsistent system)
        
    Algorithm:
    1. Find unassigned variables
    2. Look for relations that can determine values for unassigned variables:
       - Zero relations: assign 0
       - Nunity relations: assign -1
       - Alias relations: copy from canonical variable
       - Sum alias relations: compute sum of assigned variables
    3. Recursively continue until all variables assigned
    
    This implements a constraint propagation algorithm.
    """
    assigned = partial_assignment.keys()
    unassigned = [x for x in all_variables(reduced_relations) if x not in assigned]
    
    if verbose:
        print("EXTENDING VARIABLE ASSIGNMENT")
        print(f"\tpartial assignment: {partial_assignment}")
        print(f"\tunassigned variables: {unassigned}")
        
    match unassigned:
        case []:
            # Base case: all variables assigned
            return partial_assignment
        case _:
            # Find next variable to assign
            next_assignment = None
            
            for relation in reduced_relations:
                relation_type = type(relation)
                
                if relation_type == Zero:
                    state = relation.state
                    if verbose:
                        print(f"\tconsidering zero {state} := 0")
                    if state not in assigned:
                        next_assignment = [state, 0]
                        break
                        
                elif relation_type == Nunity:
                    state = relation.state
                    if verbose:
                        print(f"\tconsidering nunity {state} := -1")
                    if state not in assigned:
                        next_assignment = [state, -1]
                        break
                        
                elif relation_type == Alias:
                    alias = relation.alias
                    state = relation.state
                    if verbose:
                        print(f"\tconsidering alias {alias} := {state}")
                    if state in assigned and alias in unassigned:
                        next_assignment = [alias, partial_assignment[state]]
                        break
                        
                elif relation_type == Conservation:
                    sum_alias = relation.try_sum_alias()
                    if sum_alias is not None:
                        alias, sum_vars = sum_alias
                        if verbose:
                            print(f"\tconsidering sum alias {alias} := {sum_vars}")
                        if alias in unassigned and all([x in assigned for x in sum_vars]):
                            if verbose:
                                print(f"restoring sum alias {alias} {sum_vars}")
                            next_assignment = [alias, sum([partial_assignment[x] for x in sum_vars])]
                            break
                            
            # Apply assignment and recurse
            if next_assignment is not None:
                partial_assignment[next_assignment[0]] = next_assignment[1]
                return extend_variable_assignment(reduced_relations, partial_assignment, braid_states, verbose)
            else:
                raise Exception(f"Could not find next assignment, system may be inconsistent.\n"
                              f"\tunassigned variables: {unassigned}\n"
                              f"\tassignments: {partial_assignment}")

def violates_relation(assignment, relation, verbose=False):
    """
    Check if a variable assignment violates a specific relation.
    
    Args:
        assignment: Dict mapping variables to their values
        relation: Relation object to check
        verbose: If True, print violation details
        
    Returns:
        True if assignment violates the relation, False otherwise
        
    This function validates that a concrete assignment satisfies
    all the constraints encoded in the relation.
    """
    relation_type = type(relation)
    
    if relation_type == Leq:
        # Check A ≤ B constraint
        first = relation.first
        second = relation.second
        
        # Get values (literals have .state, variables use assignment)
        if type(first) == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if type(second) == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
            
        return (first_value > second_value)  # Violation: first > second
        
    elif relation_type == Less:
        # Check A < B constraint
        first = relation.first
        second = relation.second
        
        if type(first) == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if type(second) == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
            
        return (first_value >= second_value)  # Violation: first ≥ second
        
    elif relation_type == Alias:
        # Check A := B constraint (A = B)
        alias = relation.alias
        state = relation.state
        alias_value = assignment[alias]
        state_value = assignment[state]
        return (alias_value != state_value)  # Violation: A ≠ B
        
    elif relation_type == Conservation:
        # Check sum(inputs) = sum(outputs) constraint
        sum_alias = relation.try_sum_alias()
        if sum_alias is not None:
            # Single variable = sum case
            alias, sum_vars = sum_alias
            alias_value = assignment[alias]
            sum_value = sum([assignment[x] for x in sum_vars])
            return (alias_value != sum_value)  # Violation: alias ≠ sum
        else: 
            # General conservation case
            inputs = relation.inputs
            outputs = relation.outputs
            input_sum = sum([assignment[x] for x in inputs])
            output_sum = sum([assignment[x] for x in outputs])
            return (input_sum != output_sum)  # Violation: input_sum ≠ output_sum
            
    elif relation_type == Zero:
        # Check variable = 0 constraint
        state = relation.state
        state_value = assignment[state]
        return (state_value != 0)  # Violation: state ≠ 0
        
    elif relation_type == Nunity:
        # Check variable = -1 constraint
        state = relation.state
        state_value = assignment[state]
        return (state_value != -1)  # Violation: state ≠ -1
        
    return False  # Unknown relation type, assume satisfied

def violates_any_relation(assignment, relations, verbose=False):
    """
    Check if assignment violates any relation in the set.
    
    Args:
        assignment: Dict mapping variables to their values
        relations: List of relation objects to check
        verbose: If True, print violation details
        
    Returns:
        True if any relation is violated, False if all are satisfied
    """
    return any([violates_relation(assignment, r, verbose) for r in relations])


def print_conservations(relations, assignment=dict()):
    """
    Print conservation relations in symbolic form.
    
    Args:
        relations: List of relation objects
        assignment: Dict mapping variables to symbolic expressions
        
    Prints conservation equations as "0 = expression" for debugging.
    """
    printed = []
    for r in relations:
        if isinstance(r, Conservation):
            sum_alias = r.try_sum_alias()
            if sum_alias is not None:
                alias, sum_vars = sum_alias
                valid_list = [alias] + sum_vars  # Fixed: valid_list was undefined
                if all([v in list(assignment.keys()) for v in valid_list]):
                    expression = assignment[alias] - assignment[sum_vars[0]] - assignment[sum_vars[1]]
            else:
                valid_list = r.inputs + r.outputs
                if all([v in list(assignment.keys()) for v in valid_list]):
                    expression = assignment[r.inputs[0]] + assignment[r.inputs[1]] - assignment[r.outputs[0]] - assignment[r.outputs[1]]
            if expression not in printed:
                print("0", "=", expression)
                printed.append(expression)

def print_inequalities(relations, assignment=dict()):
    printed = []
    for r in relations:
        if isinstance(r, Less):
            if r.first == ZERO_STATE:
                first = zero
            elif r.first == NUNITY_STATE:
                first = nunity
            else:
                if r.first in assignment.keys():
                    first = assignment[r.first]
            if r.second == ZERO_STATE:
                second = zero
            elif r.second == NUNITY_STATE:
                second = nunity
            else:
                if r.second in assignment.keys():
                    second = assignment[r.second]
            if second - first not in printed:
                print(0, "<", second - first)
                printed.append(second - first)
        elif isinstance(r, Leq):
            if r.first == ZERO_STATE:
                first = zero
            elif r.first == NUNITY_STATE:
                first = nunity
            else:
                if r.first in assignment.keys():
                    first = assignment[r.first]
            if r.second == ZERO_STATE:
                second = zero
            elif r.second == NUNITY_STATE:
                second = nunity
            else:
                if r.second in assignment.keys():
                    second = assignment[r.second]
            if second - first not in printed:
                print(0, "<=", second - first)
                printed.append(second - first)

def print_angle_inequalities(criteria, relations, mapping, angle_assignment, n_crossings, assignment):
    acc = np.zeros(1 + n_crossings)
    dict_ = criteria[0].as_coefficients_dict()
    for key in dict_:
        acc += angle_assignment[tuple(mapping[key])] * dict_[key]
    print("\nAngle Criterion:")
    acc[1:] *= -1
    print(acc)
    print("\n")
    printed = []
    for r in relations:
        if isinstance(r, Less) or isinstance(r, Leq):
            if r.second == NUNITY_STATE:
                second = nunity
            elif r.second == ZERO_STATE:
                second = zero
            else:
                second = assignment[r.second]
            if r.first == NUNITY_STATE:
                first = nunity
            elif r.first == ZERO_STATE:
                first = zero
            else:
                first = assignment[r.first]
            dict_ = (second - first).as_coefficients_dict()
            acc = np.zeros(1 + n_crossings)
            for key in dict_:
                acc += angle_assignment[tuple(mapping[key])] * dict_[key]
            if list(acc) not in printed:
                print(acc)
                printed.append(list(acc[:]))

def print_inequality_in_angles (inequality, angle_assignment, n_crossings, mapping):
        dict_ = inequality.as_coefficients_dict()
        acc = np.zeros(1 + n_crossings)
        for key in dict_:
            acc += angle_assignment[tuple(mapping[key])] * dict_[key]
        print(acc)