# Configuration Files

Config-file execution is implemented by `src/fkcompute/api/batch.py:fk_from_config`.

Supported formats:

- JSON (`.json`)
- YAML (`.yaml` / `.yml`, requires `PyYAML` / `fkcompute[yaml]`)

Two modes are supported:

1. Single computation (top-level `braid` and `degree`)
2. Batch mode (top-level `computations:` list)

## Single Computation

Minimal:

```yaml
braid: [1, -2, 3]
degree: 2
```

Common options:

```yaml
braid: [1, 1, 1]
degree: 2

name: trefoil_d2
preset: parallel

symbolic: true
threads: 8

verbose: true
save_data: true
save_dir: data
```

## Batch Mode

Batch config example:

```yaml
# Global defaults (applied to every entry unless overridden)
preset: parallel
threads: 8
symbolic: false
save_data: true
save_dir: data

computations:
  - name: trefoil_d2
    braid: [1, 1, 1]
    degree: 2
    symbolic: true

  - name: figure_eight_d3
    braid: [1, -2, 1, -2]
    degree: 3
    # override the global preset for this entry
    preset: single thread
```

Return value:

- Single mode returns a single result dict.
- Batch mode returns a dict keyed by computation name.

## Presets

Preset definitions live in `src/fkcompute/api/presets.py`:

- `single thread`
  - `max_workers: 1`
  - `threads: 1`
  - `chunk_size: 4096`
  - `verbose: false`

- `parallel`
  - `max_workers: cpu_count-1` (minimum 1)
  - `threads: cpu_count-1` (minimum 1)
  - `chunk_size: 16384`
  - `verbose: true`

You can still override individual parameters alongside a preset.

## Precomputed Inputs

### Inversion Data

You can provide inversion data to skip the Phase 1 search.

In config files, `inversion` is expected as a mapping from component index to a list of signs:

```yaml
braid: [1, -2, 1, -2]
degree: 11

inversion:
  0: [-1, 1, -1, 1, -1, 1, -1, 1]
```

`fk_from_config` wraps this into the internal structure:

```python
{
  "inversion_data": {0: [...]},
  "braid": [...],
  "degree": 11,
}
```

You can also load inversion data from a JSON file:

```yaml
inversion_file: data/4_1_inversion.json
```

The inversion JSON on disk is expected to look like:

```json
{"inversion_data": {"0": [1, 1, 1, 1]}, "braid": [1, 1], "degree": 2}
```

### ILP File

If you already have an ILP input file (the format produced by `fkcompute.solver.ilp.ilp(...)`), you can skip both inversion and ILP generation:

```yaml
braid: [1, 1, 1]
degree: 10
ilp_file: data/trefoil_ilp.csv
threads: 8
```

Note:

- The C++ backend expects the *base path without extension*; the Python code derives this from `ilp_file`.
- Keep `braid` consistent with the ILP file contents; Python uses it for metadata and component counting.

## Partial Signs (Advanced)

The inversion search can be constrained with partial sign information.

Current format (used by `src/fkcompute/inversion/search.py`) is a dict keyed by component index:

```yaml
partial_signs:
  0: [1, 0, -1, 0]   # 1/-1 fixed, 0 (or null) means unknown
  1: [0, 0, 0]
```

This reduces the search space from `2^n` to `2^(n - fixed)`.

## Performance Knobs

- `max_workers`: Python multiprocessing workers for the inversion search (Phase 1)
- `chunk_size`: number of sign-assignment indices per worker task (trade latency vs overhead)
- `threads`: number of C++ compute threads passed to `fk_main --threads`
- `include_flip`: include horizontally flipped braid variants in the inversion search
- `max_shifts`: limit cyclic shifts tried during inversion search

For more detail on inversion behavior, see `docs/inversion.md`.
