# Domain Layer

The domain layer contains the core business logic for FK computation. This layer is **pure** - it has no I/O dependencies and only relies on numpy for numerical operations.

## Structure

```
domain/
├── braid/           # Braid-related domain objects
├── constraints/     # Constraint system and relations
└── solver/          # Sign assignment logic
```

## Design Principles

1. **No I/O**: This layer never reads files, makes network calls, or interacts with external systems
2. **Pure Functions**: Functions are deterministic and side-effect free where possible
3. **Self-Contained**: All domain logic can be tested without mocking external dependencies

## Subpackages

### `braid/`
Braid topology and state management:
- `types.py`: Core type definitions (`StateLiteral`, `ZERO_STATE`, `NUNITY_STATE`)
- `word.py`: Braid word operations and visualization
- `states.py`: `BraidStates` class for managing braid state information

### `constraints/`
Constraint system for FK computation:
- `relations.py`: Constraint classes (`Leq`, `Less`, `Zero`, `Nunity`, `Alias`, `Conservation`)
- `symbols.py`: Symbolic linear algebra (`Symbol` class)
- `reduction.py`: Constraint propagation and reduction algorithms

### `solver/`
Sign assignment logic:
- `assignment.py`: Variable assignment functions
- `validation.py`: Sign assignment validation

## Usage

```python
from fkcompute.domain import (
    # Types
    StateLiteral, ZERO_STATE, NUNITY_STATE,
    # Braid
    BraidStates, is_homogeneous_braid,
    # Constraints
    Leq, Less, Zero, Nunity, Alias, Conservation,
    # Symbols
    Symbol, symbols, solve,
    # Reduction
    full_reduce,
)

# Create a braid and get its states
braid = [1, -2, 3]
bs = BraidStates(braid)

# Get and reduce relations
relations = bs.get_state_relations()
reduced = full_reduce(relations)

# Work with symbols
a, b, c = symbols(3)
expr = a + 2*b - c
```

## Dependency Flow

```
braid/types.py (no deps)
       ↓
braid/word.py (no deps)
       ↓
constraints/relations.py → braid/types.py
       ↓
constraints/symbols.py (numpy only)
       ↓
constraints/reduction.py → relations.py, braid/types.py
       ↓
braid/states.py → constraints/*, braid/word.py
       ↓
solver/assignment.py → constraints/*, braid/types.py
       ↓
solver/validation.py → assignment.py, constraints/*
```
