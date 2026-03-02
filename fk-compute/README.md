# fkcompute

Compute the FK invariant of braid closures (knot theory) via a Python orchestration layer plus a compiled C++ backend.

- Python package: `fkcompute`
- CLI: `fk`
- Backend executable: `fk_main` (bundled into the wheel)

## Documentation

Project docs live in `docs/`:

- Getting started: `docs/quickstart.md`
- CLI: `docs/cli.md`
- Config files: `docs/config.md`
- Python API: `docs/python_api.md`
- Output format: `docs/output.md`
- Pipeline internals: `docs/pipeline.md`

Optional: build a single PDF from `docs/` with `./docs/build_pdf.sh`.

## Install

### System dependencies

You must install these system libraries to build/run the C++ backend:

- FLINT (+ GMP)
- OpenMP runtime
- BLAS implementation (OpenBLAS recommended)

macOS (Homebrew):

```bash
brew install flint libomp openblas
```

Ubuntu/Debian:

```bash
sudo apt-get install libflint-dev libgomp1-dev libopenblas-dev
```

RHEL/Fedora:

```bash
sudo yum install flint-devel gcc-openmp openblas-devel
```

### Python package

```bash
pip install .
```

Optional extras:

```bash
pip install ".[symbolic]"      # SymPy symbolic output
pip install ".[interactive]"   # Rich interactive wizard + history
pip install ".[yaml]"          # YAML configs (PyYAML)
pip install ".[full]"          # all optional extras
```

### Gurobi

`fkcompute` uses Gurobi (`gurobipy`) for ILP feasibility/boundedness checks. You need a working Gurobi install and license.

## Quickstart

### CLI

Compute quickly:

```bash
fk simple "[1,1,1]" 2
```

Symbolic output (requires `fkcompute[symbolic]`):

```bash
fk simple "[1,1,1]" 2 --symbolic
fk simple "[1,1,1]" 2 --format latex
fk simple "[1,1,1]" 2 --format mathematica
```

Run from a config file (threading, saving, presets, batch, etc.):

```bash
fk config my_run.yaml
```

Interactive wizard:

```bash
fk
fk interactive --quick
```

### Python

```python
from fkcompute import fk

result = fk([1, 1, 1], 2)
print(result["metadata"]["components"])
print(len(result["terms"]))
```

## Config files

Config execution supports JSON and YAML (YAML requires `fkcompute[yaml]`).

Minimal `my_run.yaml`:

```yaml
braid: [1, 1, 1]
degree: 2
```

Typical config with performance and I/O options:

```yaml
braid: [1, -2, 1, -2]
degree: 3

name: figure_eight_d3
preset: parallel

threads: 8
max_workers: 8
chunk_size: 16384

symbolic: true
save_data: true
save_dir: data
```

Batch mode uses a top-level `computations:` list (see `docs/config.md`).

## Output

The result is JSON with:

- `terms`: sparse FK polynomial terms
- `metadata`: auxiliary information (and optional `metadata.symbolic`)

See `docs/output.md` for the exact schema.

If you saved results to disk (`save_data: true`), you can reformat a JSON result without recomputing:

```bash
fk print-as data/trefoil.json --format latex
```

## Man page

Install and view the manual (Unix-like systems):

```bash
fk-install-man
man fk
```

## Mathematica / Wolfram Language

This repo includes a Paclet wrapper under `mathematica/FkCompute/`.

Start here: `docs/mathematica.md`.

## Development

Build the C++ backend:

```bash
cmake -B build -S .
cmake --build build
```

Run tests:

```bash
pytest
```

More details: `docs/development.md`.

## License

MIT.

## Authors

Toby, Paul, Lara, Josef, Davide.
