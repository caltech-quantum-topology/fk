# Solver Subpackage (Domain Layer)

This subpackage contains the domain logic for sign assignment and validation. It provides functions for computing symbolic variable assignments and validating sign configurations.

> **Note:** This is the domain-layer solver logic. For the ILP solver and parallel inversion search, see `fkcompute.solver/`.

## Modules

### `assignment.py`
Variable assignment functions for computing symbolic representations of state variables.

**Key Functions:**

| Function | Description |
|----------|-------------|
| `symbolic_variable_assignment(relations, braid_states)` | Create complete symbolic assignment |
| `extend_variable_assignment(relations, partial, braid_states)` | Extend partial assignment to all variables |
| `equivalence_assignment(assignment, braid_states)` | Extend to equivalent state locations |
| `find_expressions(relations, assignment, braid_states)` | Find conservation expressions |
| `minimal_free(expressions, new)` | Minimize free variables through elimination |

```python
from fkcompute.domain.solver.assignment import symbolic_variable_assignment

# Get reduced relations
relations = braid_states.reduced_relations()

# Create symbolic assignment
assignment = symbolic_variable_assignment(relations, braid_states)

# Each state location maps to a Symbol or integer
for state, value in assignment.items():
    print(f"{state}: {value}")
```

**How it works:**

1. Identify free variables from reduced relations
2. Assign fresh symbols to each free variable
3. Extend assignment through:
   - Zero constraints (state = 0)
   - Nunity constraints (state = -1)
   - Alias constraints (state_a = state_b)
   - Sum aliases from conservation (a = b + c)
4. Minimize free variables using conservation expressions

### `validation.py`
Sign assignment validation against constraint relations.

**Key Functions:**

| Function | Description |
|----------|-------------|
| `violates_relation(assignment, relation)` | Check if assignment violates a single relation |
| `violates_any_relation(assignment, relations)` | Check if any relation is violated |
| `check_sign_assignment(degree, relations, braid_states)` | Full validation for a degree |
| `czech_sign_assignment(degree, relations, braid_states)` | Alternative validation with transformed criteria |

```python
from fkcompute.domain.solver.validation import (
    check_sign_assignment, violates_any_relation
)

# Validate a sign assignment
result = check_sign_assignment(degree=2, relations=relations, braid_states=bs)

if result is None:
    print("No valid assignment exists for this degree")
else:
    print("Valid assignment found!")
    print(f"Criteria: {result['criteria']}")
    print(f"Assignment: {result['assignment']}")
```

**Validation checks:**
- Inequality constraints (Leq, Less)
- Equality constraints (Zero, Nunity, Alias)
- Conservation constraints (sum preservation)
- Degree bounds (component-wise)

## Usage Example

```python
from fkcompute.domain.braid.states import BraidStates
from fkcompute.domain.constraints.reduction import full_reduce
from fkcompute.domain.solver.assignment import symbolic_variable_assignment
from fkcompute.domain.solver.validation import check_sign_assignment

# Create braid states
bs = BraidStates([1, -2, 1, -2])  # Figure-eight knot

# For non-homogeneous braids, you need to set strand_signs first
# (normally done by the inversion search)

# Get and reduce relations
all_relations = bs.get_state_relations()
relations = full_reduce(all_relations)

# Create symbolic assignment
assignment = symbolic_variable_assignment(relations, bs)

# Check if valid for degree 3
result = check_sign_assignment(3, relations, bs)
```

## Integration with ILP Solver

The domain solver provides the symbolic representation that the ILP solver (`fkcompute.solver.ilp`) uses to:

1. Generate the ILP tableau from symbolic constraints
2. Check feasibility using Gurobi
3. Find bounded integer solutions

The validation functions are also used by the parallel inversion search (`fkcompute.solver.inversion`) to quickly filter out invalid sign configurations.
