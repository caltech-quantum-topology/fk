# ILP Generation (Gurobi + Custom CSV Format)

The ILP stage serves two purposes:

1. During inversion search (Phase 1): quickly reject sign assignments that cannot satisfy degree bounds.
2. During the main compute: produce an on-disk input file consumed by the C++ backend.

Implementation: `src/fkcompute/solver/ilp.py`.

## Gurobi Feasibility / Boundedness Check

Function:

- `integral_bounded(multiples, single_var_signs) -> bool`

Inputs:

- `multiples`: list of symbolic expressions representing inequalities of the form `0 <= expr`
- `single_var_signs`: mapping `Symbol(i) -> (+1.0|-1.0)`

Notes:

- The solver flips variables with negative sign into a positive form by transforming the constraint tableau columns.
- It uses a quiet Gurobi environment (`OutputFlag = 0`) created at import time.

## ILP File Writer

Function:

- `ilp(degree, relations, braid_states, write_to=...) -> str | None`

This generates a custom file (often named `*_ilp.csv`) that is *not* a standard LP/MPS format.
It is a structured CSV-like stream expected by the C++ backend.

If the sign assignment is invalid for the given degree, `ilp(...)` returns `None`.

## File Format (High Level)

The file is line-based and uses `/` as section separators.

Header lines:

1. `degree,`
2. `n_components,`
3. `writhe,`

Braid encoding:

4. One line with `2 * n_crossings` entries: for each crossing
   - `abs(generator_index)`
   - R-matrix type digit (the second character of strings like `R1`, so `1..4`)

Component mapping lines:

5. Closed strand components for strands `1..n_strands-1`
6. For each crossing: `(top_component, bottom_component)`

Then three tableaus:

7. Criteria tableau rows (degree constraints)
8. `/`
9. Inequality tableau rows
10. `/`
11. Assignment tableau rows

All numeric rows end with a trailing comma.

For an example file, see `data/4_1_ilp.csv`.

## Generating An ILP File Programmatically

```python
from fkcompute.domain.braid.states import BraidStates
from fkcompute.domain.constraints.reduction import full_reduce
from fkcompute.solver.ilp import ilp

bs = BraidStates([1, 1, 1])
relations = full_reduce(bs.get_state_relations())

csv_text = ilp(degree=5, relations=relations, braid_states=bs, write_to="trefoil_ilp.csv")
assert csv_text is not None
```

If you need a human-readable dump of the same symbolic objects:

- `src/fkcompute/solver/ilp.py:print_symbolic_relations`

There is a runnable example in `src/fkcompute/example_symbolic.py`.

## How The C++ Backend Uses The File

The C++ executable `fk_main` reads the file as `<base>.csv`.

Python calls it via:

- `src/fkcompute/infra/binary.py:run_fk_binary`

See `docs/cpp_backend.md`.
