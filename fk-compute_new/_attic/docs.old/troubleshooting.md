# Troubleshooting

## Build Errors (C++ Backend)

Symptom: `pip install .` fails while compiling.

Common causes:

- FLINT headers/libs not found
- OpenMP not available
- BLAS / CBLAS missing
- CMake too old

Fixes:

- Install the system dependencies listed in `docs/installation.md`.
- On macOS, ensure `brew install libomp` is present.

## Runtime Error: Helper Binary Not Found

Symptom:

- Python raises `RuntimeError: Could not find required helper binary 'fk_main' ...`

Fixes:

- Ensure the package built successfully (wheel contains `fkcompute/_bin/fk_main`).
- If working from source, run `cmake -B build -S . && cmake --build build`.
- If you have a custom `fk_main`, ensure it is discoverable via `PATH`.

## Gurobi License / Import Errors

Symptom:

- Importing `gurobipy` fails
- Solves fail with licensing errors

Fixes:

- Verify your Gurobi installation and license.
- Set `GRB_LICENSE_FILE` if needed.
- Quick sanity check:

```bash
python -c "import gurobipy as gp; m=gp.Model(); print('ok')"
```

## YAML Config Fails To Load

Symptom:

- `ImportError: PyYAML is required for YAML config files`

Fix:

```bash
pip install "fkcompute[yaml]"
```

## Symbolic Output Not Available

Symptom:

- CLI prints: `SymPy is required for symbolic output`

Fix:

```bash
pip install "fkcompute[symbolic]"
```

Note: `pretty` output may include unicode. If you need ASCII-only output, use:

```bash
fk simple "[1,1,1]" 2 --format inline
```

## Enhanced Interactive / History Not Available

Symptom:

- CLI falls back to basic mode
- `fk history ...` reports missing dependencies

Fix:

```bash
pip install "fkcompute[interactive]"
```

## Symbolic Conversion Fails On Huge Coefficients

`src/fkcompute/output/symbolic.py` refuses to parse coefficients with more than 10,000 digits.

Workarounds:

- Use JSON output and post-process numerically.
- Use `--symbolic` only for smaller runs.

## Man Page Install Fails

Symptom:

- `fk-install-man` cannot write to system directories.

Notes:

- The installer tries user-local man directories first (`~/.local/share/man/man1`).
- If you want system-wide installation, you may need `sudo`.
