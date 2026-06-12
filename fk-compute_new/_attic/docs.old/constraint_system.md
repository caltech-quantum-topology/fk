# Constraint System (Domain Layer)

The constraint system is the core "math logic" in the Python layer. It lives under `src/fkcompute/domain/` and is designed to be mostly I/O-free.

## Key Objects

- `BraidTopology`: pure topological data (components, state locations, etc)
  - `src/fkcompute/domain/braid/topology.py`

- `SignedBraid`: topology + sign assignment + relation generation
  - `src/fkcompute/domain/braid/signed.py`

- `BraidStates`: backwards-compatible wrapper used throughout the package
  - `src/fkcompute/domain/braid/states.py`

## Relations (Raw Constraints)

`SignedBraid.get_state_relations()` returns a list of relation objects (constraints).

Code: `src/fkcompute/domain/braid/signed.py:get_state_relations`.

It combines four sources:

1. Boundary conditions
   - The first strand endpoints are fixed to `[0]` (or `[-1]`) depending on the first component sign.

2. Periodicity
   - Closure aliases: `(i, len(braid)) := (i, 0)`.

3. Sign bounds
   - Every state is bounded above/below by `[0]` or `[-1]` depending on the sign assignment.

4. Crossing constraints
   - Inequalities determined by the computed R-matrix type (R1/R4 currently emit explicit inequalities).
   - Conservation constraints at every crossing.

Relation classes live in `src/fkcompute/domain/constraints/relations.py`:

- `Leq`, `Less`
- `Zero`, `NegOne`
- `Alias`
- `Conservation`

## Reduction

Raw constraints are often redundant.

The reducer `full_reduce(...)` applies a fixed-point simplification pass over relations:

- propagate aliases through inequalities and conservation
- detect symmetric inequalities (`a <= b` and `b <= a` implies `a := b`)
- propagate `Zero` / `NegOne` through aliases
- remove vacuous constraints

Implementation: `src/fkcompute/domain/constraints/reduction.py`.

Snapshot tests pin this behavior in `tests/test_phase2_snapshots.py`.

## Symbolic Assignment

After reduction, the system assigns symbolic expressions to states.

Implementation:

- `src/fkcompute/domain/solver/assignment.py:symbolic_variable_assignment`

This produces a dict:

```python
{
  (state_location): Symbol(...) | 0 | -1,
  ...
}
```

Free variables correspond to reduced states that are not forced to a constant.

Helper:

- `src/fkcompute/domain/constraints/reduction.py:free_variables`

## Criteria + Inequalities

Given the symbolic assignment, `process_assignment(...)` extracts:

- degree criteria per component
- multi-variable inequalities (0 <= expr)
- single-variable signs (+1 / -1)

Implementation:

- `src/fkcompute/domain/solver/symbolic_constraints.py:process_assignment`

These are the inputs to both:

- inversion validity checks (Phase 1)
- ILP file generation (Phase 2)

## One-Call Helper

If you want to run the whole Phase 2 flow in one call, use:

- `src/fkcompute/domain/constraints/pipeline.py:build_constraint_system`

Example:

```python
from fkcompute.domain.braid.states import BraidStates
from fkcompute.domain.constraints.pipeline import build_constraint_system

bs = BraidStates([1, 1, 1])
cs = build_constraint_system(bs)

print(cs.degree_criteria)
print(cs.multi_var_inequalities)
print(cs.single_var_signs)
```

## Visualizing Braids (Optional)

The domain layer has a simple braid SVG generator:

- `src/fkcompute/domain/braid/word.py:braid_svg`
- `src/fkcompute/domain/braid/word.py:draw_braid` (Jupyter/IPython)
