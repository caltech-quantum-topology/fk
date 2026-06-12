# Braid Subpackage

This subpackage contains braid-related domain objects for FK computation.

## Modules

### `types.py`
Core type definitions used throughout the FK computation.

**Classes:**
- `StateLiteral`: Represents literal state values (0 or -1) in the constraint system

**Constants:**
- `ZERO_STATE`: The zero state literal `[0]`
- `NEG_ONE_STATE`: The negative one state literal `[-1]`

```python
from fkcompute.domain.braid.types import StateLiteral, ZERO_STATE, NEG_ONE_STATE

# Check if something is a zero state
if state == ZERO_STATE:
    print("State is zero")
```

### `word.py`
Braid word operations and utilities.

**Functions:**
- `is_positive_braid(braid)`: Check if all crossings are positive
- `is_homogeneous_braid(braid)`: Check if braid is homogeneous (no generator appears with both signs)
- `strand_path(strand, braid)`: Compute path of a strand through the braid
- `braid_svg(braid, ...)`: Generate SVG visualization
- `draw_braid(braid, ...)`: Display braid in Jupyter notebook

```python
from fkcompute.domain.braid.word import is_homogeneous_braid, is_positive_braid

braid = [1, 1, 1]  # Trefoil
print(is_positive_braid(braid))     # True
print(is_homogeneous_braid(braid))  # True

braid = [1, -2, 1, -2]  # Figure-eight
print(is_positive_braid(braid))     # False
print(is_homogeneous_braid(braid))  # True
```

### `states.py`
The main `BraidStates` class for managing braid state information.

**Class: `BraidStates`**

Represents a braid and its associated state information, including:
- State locations on the braid diagram
- Component structure (which strands form which link components)
- Sign assignments at crossings
- Relations between states

```python
from fkcompute.domain.braid.states import BraidStates

# Create braid states for a trefoil
bs = BraidStates([1, 1, 1])

# Access properties
print(bs.n_components)      # Number of link components
print(bs.n_crossings)       # Number of crossings
print(bs.writhe)            # Writhe of the braid
print(bs.strand_signs)      # Sign assignments

# Get state relations
relations = bs.get_state_relations()

# Get reduced relations
reduced = bs.reduced_relations()

# Get free variables
free_vars = bs.free_variables()
```

## Braid Conventions

- **Strand numbering**: Strand 0 is the top strand, strand 1 is second from top, etc.
- **Generator indexing**: Generators are indexed from 1 to (number of strands - 1)
- **Positive crossings**: Positive integer `n` means strand `n-1` crosses over strand `n`
- **Negative crossings**: Negative integer `-n` means strand `n-1` crosses under strand `n`

## Examples

```python
# Trefoil knot (positive braid)
trefoil = [1, 1, 1]

# Figure-eight knot (alternating)
figure_eight = [1, -2, 1, -2]

# Hopf link
hopf = [1, 1]

# Create BraidStates and analyze
bs = BraidStates(trefoil)
print(f"Components: {bs.n_components}")
print(f"Crossings: {bs.n_crossings}")
print(f"Writhe: {bs.writhe}")
```
