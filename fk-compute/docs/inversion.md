# Inversion (Sign Assignment Search)

Phase 1 of the pipeline finds a valid sign assignment for the braid (or declares failure).

Entry point: `src/fkcompute/inversion/api.py:find_sign_assignment`.

## Homogeneous vs Non-Homogeneous

Homogeneous check:

- `src/fkcompute/domain/braid/word.py:is_homogeneous_braid`

Definition used here:

- A braid is homogeneous if no generator index appears with both signs (i.e., never both `k` and `-k`).

Behavior:

- Homogeneous braids: signs are initialized deterministically by `src/fkcompute/domain/braid/signed.py:SignedBraid._init_homogeneous_signs`.
- Non-homogeneous braids: `parallel_try_sign_assignments(...)` searches over sign configurations.

## What Is Being Searched?

The braid closure decomposes into components and each component has a number of strand segments.

The total number of signs to choose is:

- `BraidStates(...).n_s_total` (see `src/fkcompute/domain/braid/topology.py`)

Each sign assignment is encoded as an integer index in `[0, 2^n)`.

Decoder:

- `src/fkcompute/inversion/search.py:sign_assignment_from_index`

## Braid Variants: Rotate + (Optional) Flip

The search may try braid variants to find a configuration that validates:

- cyclic rotation (move the first generator to the end)
- horizontal flip (mirror generators)

Variant generator:

- `src/fkcompute/inversion/variants.py:generate_braid_variants`

Controls:

- `include_flip` (doubles the variant count)
- `max_shifts` (limit rotations; `None` means all shifts)

If you provide `partial_signs`, the code rotates/flips the partial assignment alongside the braid.

## Partial Signs

Partial signs reduce the search space.

Type used by the inversion code:

- `PartialSignsType = Dict[int, List[Optional[int]]]` in `src/fkcompute/inversion/variants.py`

Conventions:

- `1` / `-1` means fixed
- `0` or `null`/`None` means unknown (will be filled by the search)

Example YAML:

```yaml
partial_signs:
  0: [1, 0, -1, 0]
```

## Validation: What Makes An Assignment "Valid"?

For each candidate assignment:

1. Load signs into `BraidStates` and compute crossing matrices
2. Validate that matrices and R-matrix type classification are consistent
3. Generate raw relations, reduce them, build a symbolic assignment
4. Check degree-bounded feasibility via Gurobi

Core worker logic:

- `src/fkcompute/inversion/search.py:check_assignment_for_braid`
- `src/fkcompute/inversion/validation.py:check_sign_assignment`
- `src/fkcompute/solver/ilp.py:integral_bounded`

In other words: inversion search is not just "pick +/-"; it filters by the same constraint logic used later for ILP generation.

## Parallelization

Parallelism is Python multiprocessing across disjoint index ranges.

Controls:

- `max_workers`: number of worker processes
- `chunk_size`: indices per task (half-open ranges)

Tradeoffs:

- Larger `chunk_size` reduces overhead but increases "time to first hit" because workers check longer ranges before returning.
- Larger `max_workers` increases throughput but can increase memory and (depending on your setup) the number of concurrent Gurobi environments.

## Debugging Tips

### Print Variant Attempts

Pass `verbose=True` to `find_sign_assignment(...)` (or via the main `fk(...)` call).
`parallel_try_sign_assignments` will print which variant is being tried.

### Reproduce Inversion From Python

```python
from fkcompute.inversion.api import find_sign_assignment

braid = [1, -2, 1, -2]
res = find_sign_assignment(braid, degree=3, max_workers=4, verbose=True)
print(res.success)
print(res.braid)
print(res.sign_assignment)
```

### Skip Search

If you already know inversion data, pass it (or load from file) to skip this entire phase.

See `docs/config.md` for the config-file representation.
