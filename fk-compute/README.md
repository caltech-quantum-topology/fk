# fkcompute

`fkcompute` computes the FK (Gukov-Manolescu) invariant `F_K(x,q)` of knots and links presented as braid closures.
It combines a Python orchestration layer (sign/inversion search + constraint/ILP preparation + formatting) with a compiled C++ backend for the expensive state-sum/R-matrix evaluation.

- Python package: `fkcompute`
- CLI: `fk`
- Helper binary: `fk_main` (built and bundled at install time)

## What is the FK invariant?

`F_K` is a two-variable quantum invariant for knot and link complements introduced by Gukov and Manolescu (arXiv:1904.060597).
This repository implements the large-color R-matrix / inverted state-sum approach developed by Park (arXiv:2004.02087 and followups).
In practice, you provide a braid word and a truncation degree, and `fkcompute` returns a sparse polynomial/series in `q` and in one `x`-variable per link component.

## Installation

### System dependencies (C++ backend)

To build/run the backend you need these system libraries:

- FLINT (and GMP)
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

Build tooling:

- A C++ compiler (Clang or GCC)
- CMake >= 3.20

### Install the Python package

From the repo root:

```bash
pip install .
```

Optional extras:

```bash
pip install ".[symbolic]"      # SymPy symbolic formatting (pretty/latex/mathematica)
pip install ".[interactive]"   # Rich-based interactive wizard
pip install ".[yaml]"          # YAML config files (PyYAML)
pip install ".[full]"          # all optional extras
```

### Gurobi

`fkcompute` uses Gurobi (`gurobipy`) for ILP feasibility/boundedness checks.
You need a working Gurobi installation and license (often configured via `GRB_LICENSE_FILE`).

## Quickstart

### CLI (updated syntax)

Help:

```bash
fk --help
fk simple --help
```

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

Interactive wizard:

```bash
fk
fk interactive --quick
```

Run from one or more config files:

```bash
fk config my_run.yaml
fk config a.yaml b.yaml c.json
```

Create a starter config:

```bash
fk template create my_run.yaml
```

Reformat a saved JSON result (no recomputation):

```bash
fk print-as data/trefoil.json --format inline
```

Braid input strings can be JSON-style, comma-separated, or space-separated. Always quote the argument so `-2` is not parsed as a flag:

```bash
fk simple "1 -2 3" 2
```

Legacy shortcut (still accepted):

```bash
fk "[1,-2,3]" 2
```

### Python

```python
from fkcompute import fk

result = fk([1, 1, 1], 2)
print(result["metadata"]["components"])
print(len(result["terms"]))
```

Config-file mode:

```python
result = fk("my_run.yaml")
```

### Minimal config

`my_run.yaml`:

```yaml
braid: [1, 1, 1]
degree: 2
```

Batch mode uses a top-level `computations:` list; see `docs/config.md`.

## Documentation

Project docs live in `docs/` (Sphinx / ReadTheDocs):

- Getting started: `docs/quickstart.rst`
- Installation: `docs/installation.rst`
- CLI: `docs/cli.rst`
- Config files: `docs/configuration.rst`
- Output format: `docs/output.rst`
- Pipeline internals: `docs/pipeline.rst`
- Troubleshooting: `docs/troubleshooting.rst`

Build locally:

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
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

```bash
pip install -e .
cmake -B build -S .
cmake --build build
pytest
```

More details: `docs/development.md`.

## License

MIT.

## Authors

Paul Orland, Davide Passaro, Lara San Martin Suarez, Toby Saunders-A'Court, Josef Svoboda.
