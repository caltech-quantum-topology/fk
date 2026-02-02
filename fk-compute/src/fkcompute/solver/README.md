# Solver Package

This package contains the ILP solver and parallel inversion search for FK computation. It depends on the domain layer and uses Gurobi for integer linear programming.

## Modules

### `ilp.py`
ILP (Integer Linear Programming) formulation and Gurobi integration.

**Key Functions:**

| Function | Description |
|----------|-------------|
| `ilp(degree, relations, braid_states, write_to)` | Generate ILP formulation |
| `integral_bounded(multiples, singlesigns)` | Check if system has bounded integer solution |
| `print_symbolic_relations(degree, relations, braid_states)` | Print relations in human-readable format |

```python
from fkcompute.solver.ilp import ilp, integral_bounded

# Generate ILP data for a braid
ilp_data = ilp(
    degree=2,
    relations=reduced_relations,
    braid_states=bs,
    write_to="output.csv"
)

# Check if constraints have bounded solution
is_bounded = integral_bounded(multiples, singlesigns)
```

**ILP Output Format:**

The generated ILP file is a CSV with:
1. Header: degree, n_components, writhe
2. Braid data: crossing indices and R-matrix types
3. Component mappings
4. Criteria tableau (degree constraints)
5. Inequality tableau (relation constraints)
6. Assignment tableau (variable mappings)

### `inversion.py`
Parallel sign assignment search for non-homogeneous braids.

**Key Functions:**

| Function | Description |
|----------|-------------|
| `get_sign_assignment_parallel(braid, ...)` | Main entry point for sign search |
| `parallel_try_sign_assignments(degree, braid_state, ...)` | Parallel search across variants |
| `braid_type(braid)` | Classify as homogeneous or fibered |

**Classes:**

| Class | Description |
|-------|-------------|
| `BraidType` | Enum: `HOMOGENEOUS` or `FIBERED` |

```python
from fkcompute.solver.inversion import (
    get_sign_assignment_parallel,
    braid_type,
    BraidType
)

# Check braid type
braid = [1, -2, 1, -2]
bt = braid_type(braid)
print(bt == BraidType.HOMOGENEOUS.value)  # False (it's fibered)

# Find valid sign assignment
result = get_sign_assignment_parallel(
    braid=braid,
    degree=3,
    verbose=True,
    max_workers=4,      # Parallel workers
    chunk_size=1 << 14, # Index chunk size
    include_flip=False, # Try flipped variants
    max_shifts=None,    # Cyclic rotation limit
)

if result["inversion_data"] == "failure":
    print("No valid assignment found")
else:
    print(f"Found: {result['inversion_data']}")
```

## Search Strategy

The parallel inversion search works as follows:

1. **Braid Classification**: Homogeneous braids have automatic sign assignments; fibered braids require search

2. **Variant Generation**: For each braid, generate variants through:
   - Cyclic rotations (moving first crossing to end)
   - Horizontal flips (optional)

3. **Index Space**: Each sign configuration is encoded as an integer index
   - Total space: `2^n_s_total` where `n_s_total` is total number of strand segments

4. **Parallel Search**:
   - Split index space into chunks
   - Process chunks in parallel using multiprocessing
   - Return first valid assignment found

5. **Validation**: Each candidate is validated by:
   - Computing crossing matrices
   - Checking R-matrix consistency
   - Verifying degree bounds via ILP feasibility

## Configuration Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_workers` | CPU count | Number of parallel processes |
| `chunk_size` | 16384 | Indices per worker task |
| `include_flip` | False | Include horizontally flipped variants |
| `max_shifts` | None | Limit cyclic rotations (None = all) |
| `verbose` | False | Print progress information |

## Performance Tips

1. **Chunk size**: Larger chunks reduce overhead but increase latency to first result
2. **Workers**: More workers help for large search spaces
3. **Flip symmetry**: Disable unless needed (doubles search space)
4. **Partial signs**: If you know some signs, provide them to reduce search space

```python
# With partial sign constraints
partial_signs = {
    0: [1, None, -1, None],  # Component 0: some signs known
    1: [None, None, None],   # Component 1: all unknown
}

result = get_sign_assignment_parallel(
    braid=braid,
    partial_signs=partial_signs,
    degree=3,
)
```
