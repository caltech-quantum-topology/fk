# Constraints Subpackage

This subpackage implements the constraint system used in FK computation. It provides classes for representing different types of constraints and functions for manipulating and reducing constraint systems.

## Modules

### `relations.py`
Constraint relation classes representing different types of constraints in the FK system.

**Classes:**

| Class | Description | Example |
|-------|-------------|---------|
| `Leq` | Less than or equal: `first <= second` | `Leq(a, b)` means a ≤ b |
| `Less` | Strict less than: `first < second` | `Less(a, b)` means a < b |
| `Zero` | State equals zero: `state = [0]` | `Zero(a)` means a = 0 |
| `Nunity` | State equals -1: `state = [-1]` | `Nunity(a)` means a = -1 |
| `Alias` | Two states are equal | `Alias(a, b)` means a = b |
| `Conservation` | Sum conservation: `Σinputs = Σoutputs` | Conservation at crossings |

```python
from fkcompute.domain.constraints.relations import (
    Leq, Less, Zero, Nunity, Alias, Conservation
)

# Create constraints
c1 = Leq(state_a, state_b)      # a <= b
c2 = Zero(state_c)               # c = 0
c3 = Conservation([a, b], [c, d])  # a + b = c + d
```

### `symbols.py`
Symbolic linear algebra for constraint manipulation.

**Class: `Symbol`**

Represents a linear expression with integer coefficients, like `a + 2b - 3c + 5`.

```python
from fkcompute.domain.constraints.symbols import Symbol, symbols, one, zero, solve

# Create symbols
a, b, c = symbols(3)

# Arithmetic operations
expr = a + 2*b - c + 5
expr2 = 3*a - b

# Check properties
print(expr.is_constant())     # False
print(expr.free_symbols())    # [Symbol(1), Symbol(2), Symbol(3)]
print(expr.constant())        # 5.0

# Substitution
result = expr.subs({a: 2})    # 2 + 2b - c + 5 = 7 + 2b - c

# Solve for a variable
# Given expr = 0, solve for a
solution = solve(expr, a)     # Returns [expression for a]
```

**Constants:**
- `one`: The constant 1
- `zero`: The constant 0
- `nunity`: The constant -1

### `reduction.py`
Constraint propagation and reduction algorithms.

**Key Functions:**

| Function | Description |
|----------|-------------|
| `full_reduce(relations)` | Fully reduce constraints until fixed point |
| `reduce_relations(relations)` | Apply all reduction rules once |
| `free_variables(relations)` | Get unbound variables |
| `all_variables(relations)` | Get all variables in constraints |

**Reduction Rules Applied:**
1. `de_alias_inequalities`: Replace aliased variables in inequalities
2. `symmetric_inequality`: Detect `a ≤ b` and `b ≤ a` implies `a = b`
3. `propagate_zero_aliases`: If `a = 0` and `a = b`, then `b = 0`
4. `propagate_nunity_aliases`: Same for `-1`
5. `conservation_alias`: Replace aliases in conservation constraints
6. `conservation_zeros`: Remove zeros from conservation sums
7. `unary_conservation_is_alias`: `a = b` conservation becomes alias
8. `remove_vacuous`: Remove tautologies like `a ≤ a`

```python
from fkcompute.domain.constraints.reduction import (
    full_reduce, free_variables, all_variables
)

# Get relations from a braid
relations = braid_states.get_state_relations()

# Fully reduce
reduced = full_reduce(relations)

# Find free variables
free_vars = free_variables(reduced)
print(f"Free variables: {len(free_vars)}")
```

## Constraint System Overview

The FK computation uses a constraint system to track relationships between state values at different positions on the braid diagram:

1. **Boundary conditions**: First strand starts at 0 (or -1 for negative braids)
2. **Periodicity**: States wrap around the braid closure
3. **Crossing relations**: Conservation and ordering constraints at each crossing
4. **Sign constraints**: States are bounded by 0 or -1 based on strand signs

The reduction process simplifies these constraints to find the minimal set of free variables needed to describe all possible state configurations.
