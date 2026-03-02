# Python API

The public Python entrypoint is `fkcompute.fk`.

- Re-exported from `src/fkcompute/__init__.py`
- Implemented in `src/fkcompute/api/compute.py`

## Basic Usage

```python
from fkcompute import fk

result = fk([1, 1, 1], 2)
print(result["metadata"]["braid"])
print(result["metadata"]["components"])
```

## Config-File Mode

If the first argument is a string, `fk(...)` treats it as a config path and delegates to `fk_from_config`:

```python
result = fk("config.yaml")
```

Batch config files return a dict of results keyed by computation name.

See `docs/config.md`.

## Direct Mode Parameters

Direct mode is:

```python
fk(braid: list[int], degree: int, **kwargs)
```

Common keyword arguments (see `_fk_compute` in `src/fkcompute/api/compute.py`):

- `symbolic: bool` (default: False)
- `threads: int | None` (default: 1)
- `name: str | None` (used for saved files when `save_data=True`)
- `verbose: bool` (default: False)

Performance / search knobs:

- `max_workers: int` (default: 1)
- `chunk_size: int` (default: 16384)
- `include_flip: bool` (default: False)
- `max_shifts: int | None` (default: None)

I/O knobs:

- `save_data: bool` (default: False)
- `save_dir: str` (default: "data")

Precomputed inputs:

- `inversion: dict | None`
- `inversion_file: str | None`
- `ilp_file: str | None`

Advanced:

- `partial_signs: dict | None` (see `docs/inversion.md`)
- `_progress_callback`: internal hook used by the Rich interactive UI

Notes:

- Presets (`preset: "parallel"`, etc.) are supported in config files, not as a direct-mode kwarg.
- If you want to use a preset in Python, merge `fkcompute.PRESETS["parallel"]` into your kwargs explicitly.

## Return Value

Single mode returns a dict with `terms` and `metadata`.

See `docs/output.md` for schema details.

## Calling Lower-Level Pieces

The API is intentionally layered. You can call pieces directly when debugging:

- Inversion search:
  - `fkcompute.inversion.api.find_sign_assignment`

- Constraint system (Phase 2):
  - `fkcompute.domain.constraints.pipeline.build_constraint_system`

- ILP writer:
  - `fkcompute.solver.ilp.ilp`
  - `fkcompute.solver.ilp.print_symbolic_relations`

- Symbolic formatting:
  - `fkcompute.output.symbolic.format_symbolic_output`
