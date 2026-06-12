# Installation

This project ships a Python package (`fkcompute`) that builds and bundles a compiled C++ helper binary (`fk_main`).

## System Dependencies

You must install these *system* libraries before `pip install` can build the C++ backend:

- FLINT (and its dependency GMP)
- OpenMP runtime
- BLAS implementation (OpenBLAS recommended)

Common installs:

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

## Python Installation

From the repo root:

```bash
pip install .
```

Optional extras:

```bash
pip install ".[symbolic]"      # SymPy symbolic output
pip install ".[interactive]"   # Rich-based interactive wizard + history
pip install ".[yaml]"          # YAML config files (PyYAML)
pip install ".[full]"          # all optional extras
```

Development install:

```bash
pip install -e .
```

## Gurobi

`fkcompute` uses Gurobi via `gurobipy` for ILP feasibility/boundedness checks.

You need:

- The `gurobipy` Python package (installed automatically by `pip install .`)
- A working Gurobi installation and license

Typical environment variables you may need (depends on your setup):

- `GRB_LICENSE_FILE` (path to your license file)

If Gurobi is not licensed correctly, imports/solves in `src/fkcompute/solver/ilp.py` will fail.

## Post-Install: Man Page

Install the man page (Unix-like systems):

```bash
fk-install-man
man fk
```

## Building The C++ Backend Manually (Optional)

The normal `pip install` path builds the binary automatically via scikit-build-core.

If you want to build the backend yourself:

```bash
cmake -B build -S .
cmake --build build
```

This produces `fk_main` and copies it into `src/fkcompute/_bin/fk_main` (see `CMakeLists.txt`).

There is also a standalone Makefile in `cpp/Makefile` (useful for C++ dev/debug):

```bash
make -C cpp main
```
