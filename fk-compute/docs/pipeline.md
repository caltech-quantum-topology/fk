# Computation Pipeline

This document explains what happens when you call the public entrypoint `fkcompute.fk(...)`.

Primary implementation: `src/fkcompute/api/compute.py`.

## Entry Modes

`fkcompute.fk` supports two invocation styles:

1. Config mode

```python
fk("config.yaml")
fk(braid_or_config=[...], config="config.yaml")
```

Config handling is delegated to `src/fkcompute/api/batch.py:fk_from_config`.

2. Direct mode

```python
fk([1, -2, 3], degree=2, threads=4, symbolic=True)
```

This flows into the internal `_fk_compute(...)` helper.

## Stages

### 0) Parameter Defaults + Presets

- Defaults are applied in `fk(...)` and `_fk_compute(...)`.
- Presets are merged in config mode via `src/fkcompute/api/presets.py`.

Important knobs:

- `max_workers`, `chunk_size`: Python multiprocessing for inversion search
- `threads`: passed to `fk_main --threads`
- `include_flip`, `max_shifts`: inversion variant search space
- `save_data`, `save_dir`, `name`: file outputs

### 1) Inversion (Sign Assignment)

Goal: determine a consistent +/- sign assignment for strand segments.

Code:

- `src/fkcompute/inversion/api.py:find_sign_assignment`
- `src/fkcompute/inversion/search.py:parallel_try_sign_assignments`

Behavior:

- If `ilp_file` is provided, inversion is skipped.
- If `inversion` is provided, it is used directly.
- If `inversion_file` is provided, it is loaded (JSON) and used.
- Otherwise, `find_sign_assignment(...)` runs:
  - homogeneous braids: immediate
  - non-homogeneous braids: parallel search

The resulting object is stored as a dict with at least:

```python
{
  "inversion_data": {component_index: [signs...]},
  "braid": [...],
  "degree": degree,
}
```

### 2) Constraint System + ILP File

Goal: turn braid topology + signs into a reduced constraint system and write a custom ILP-like file that the C++ backend consumes.

Key code paths:

- Relations generation: `src/fkcompute/domain/braid/signed.py:get_state_relations`
- Reduction: `src/fkcompute/domain/constraints/reduction.py:full_reduce`
- Symbolic assignment: `src/fkcompute/domain/solver/assignment.py:symbolic_variable_assignment`
- Criteria extraction: `src/fkcompute/domain/solver/symbolic_constraints.py:process_assignment`
- ILP file writer: `src/fkcompute/solver/ilp.py:ilp`

In the main compute path, ILP generation is wrapped by `_generate_ilp(...)` in `src/fkcompute/api/compute.py`.

Skip behavior:

- If `ilp_file` is provided, it is read and used.
- If `ilp` is provided (string), ILP generation is skipped, but the C++ backend still needs an on-disk `*.csv` file.
  In practice, prefer `ilp_file` if you want to skip ILP generation.

### 3) Run The C++ Backend

The backend binary is `fk_main`.

Python launches it via `src/fkcompute/infra/binary.py:run_fk_binary`:

- input: ILP base path (no extension)
- output: output base path (no extension)
- option: `--threads N`

The backend reads `<input>.csv` and writes `<output>.json`.

Binary resolution:

1. If a `fk_main` exists in your `PATH`, it is used.
2. Otherwise, Python uses the packaged binary at `fkcompute/_bin/fk_main`.

See `docs/cpp_backend.md`.

### 4) Assemble Result + Optional Symbolic Formatting

Python reads the JSON produced by the backend, then adds metadata:

- `metadata.braid`
- `metadata.inversion`
- `metadata.components`

If `symbolic=True` and SymPy is installed, Python also adds `metadata.symbolic` (pretty format by default).

Formatting code: `src/fkcompute/output/symbolic.py`.

### 5) Cleanup vs Save

If `save_data=False` (default), intermediate files are removed.

If `save_data=True`, these files are written under `save_dir`:

- `<name>_inversion.json`
- `<name>_ilp.csv`
- `<name>.json`

## Progress Hooks (Interactive UI)

`_fk_compute(...)` accepts a private `_progress_callback` argument.

The Rich interactive UI passes a `FKProgressTracker` (`src/fkcompute/interactive/progress.py`) to get live phase updates.
