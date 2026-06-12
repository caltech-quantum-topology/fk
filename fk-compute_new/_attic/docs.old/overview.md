# Project Overview

`fkcompute` is a Python package + CLI for computing FK invariants of braid closures. It combines:

- Python orchestration (`src/fkcompute/`) for parsing inputs, generating constraints/ILP data, and formatting output.
- A compiled C++ backend (`fk_main`) for the heavy polynomial computation.

Primary entrypoints:

- CLI: `fk` (declared in `pyproject.toml` as `fk = "fkcompute.cli.app:main"`)
- Python API: `fkcompute.fk(...)` (implemented in `src/fkcompute/api/compute.py`)
- Mathematica wrapper (Paclet): `mathematica/FkCompute/` calling `python -m fkcompute.mathematica_bridge`

## Repository Layout

Key folders:

```
src/fkcompute/
  api/             # public compute entrypoint + config runner + presets
  cli/             # Typer CLI commands
  domain/          # pure braid/constraint logic
  inversion/       # Phase 1 sign assignment search
  solver/          # ILP generation + Gurobi feasibility checks
  output/          # symbolic formatting (SymPy)
  infra/           # config parsing + binary execution
  interactive/     # optional Rich-based interactive UI

cpp/               # C++ sources and Makefile
mathematica/       # Wolfram Language paclet wrapper
tests/             # pytest suite and baseline data
fk.1               # man page
```

## Braid Word Conventions

Braids are represented as `list[int]` where each integer is a generator index with sign.

- Strand positions are 0-indexed in the code.
- Generators are 1..(n_strands-1).
- A crossing `+k` means the strand at position `k-1` crosses over the strand at position `k`.
- A crossing `-k` means the inverse crossing.

See `src/fkcompute/domain/braid/word.py`.

## High-Level Computation Pipeline

The main call `fkcompute.fk(...)` orchestrates:

1. Phase 1 (Inversion / sign assignment)
   - Homogeneous braids get signs automatically.
   - Non-homogeneous cases search over sign assignments and braid variants.
   - Code: `src/fkcompute/inversion/api.py`, `src/fkcompute/inversion/search.py`, `src/fkcompute/inversion/variants.py`

2. Phase 2 (Constraint system + ILP)
   - Generate relations, reduce them, build a symbolic assignment, extract degree criteria + inequalities.
   - Use Gurobi to check boundedness/feasibility and write an ILP-like CSV consumed by the C++ backend.
   - Code: `src/fkcompute/domain/constraints/*`, `src/fkcompute/domain/solver/*`, `src/fkcompute/solver/ilp.py`

3. Phase 3 (C++ FK computation)
   - Run `fk_main` on the generated input file to produce a JSON result with `terms` and `metadata`.
   - Python attaches extra metadata (braid used, inversion signs, component count) and can add symbolic formatting.
   - Code: `src/fkcompute/infra/binary.py`, `cpp/main.cpp`

4. Output formatting
   - JSON result by default.
   - Optional SymPy formatting: pretty/inline/latex/mathematica.
   - Code: `src/fkcompute/output/symbolic.py`

## Result Schema (Quick Glance)

The C++ backend writes JSON like:

```json
{
  "terms": [
    {"x": [0], "q_terms": [{"q": 1, "c": "-1"}]},
    {"x": [2], "q_terms": [{"q": 2, "c": "1"}]}
  ],
  "metadata": {
    "num_x_variables": 1,
    "max_x_degrees": [9],
    "storage_type": "flint"
  }
}
```

Python augments `metadata` with fields like `braid`, `inversion`, `components`, and (optionally) `symbolic`.

See `docs/output.md` for the full schema.
