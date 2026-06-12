# Output Format

The core result format is JSON with two top-level keys:

- `terms`: sparse polynomial terms
- `metadata`: auxiliary information

The C++ backend produces the base structure, and Python augments `metadata`.

## Top-Level Result

Single computation returns:

```python
{
  "terms": [...],
  "metadata": {...},
}
```

Batch config mode returns:

```python
{
  "name1": {"terms": [...], "metadata": {...}},
  "name2": {"terms": [...], "metadata": {...}},
}
```

## `terms`

Each entry in `terms` is a monomial in the topological variables with a polynomial in `q` as coefficient.

Schema:

```json
{
  "x": [x1_power, x2_power, ...],
  "q_terms": [
    {"q": q_exponent, "c": "integer_as_string"},
    ...
  ]
}
```

Notes:

- `x` is a list of nonnegative integers. Its length is typically `metadata.num_x_variables`.
- `q_terms` is sparse in `q`. Exponents can be negative.
- Coefficients are written as strings by the C++ backend to preserve arbitrary precision.

## `metadata`

### Produced by the C++ backend

Fields written by the backend include:

- `num_x_variables` (int)
- `max_x_degrees` (list[int])
- `storage_type` (string, e.g. `"flint"`)

### Added by Python

Python adds or updates:

- `braid`: the braid word actually used
- `inversion`: inversion data dict (or `null`)
- `components`: number of link components (computed from `BraidStates`)
- `symbolic`: pretty-printed polynomial (only if `symbolic=True` and SymPy is available)

## Symbolic Formatting

Formatter:

- `src/fkcompute/output/symbolic.py:format_symbolic_output`

Supported `format_type` values:

- `pretty`: multi-line pretty print (may include unicode)
- `inline`: one-line Python-ish string
- `latex`: LaTeX string
- `mathematica`: Mathematica syntax via SymPy

CLI integration:

- `fk simple ... --symbolic`
- `fk simple ... --format latex`
- `fk print-as result.json --format mathematica`

### Variable Naming

Variable names are chosen from `metadata.num_x_variables`:

- 1 variable: `x`
- 2 variables: `x`, `y`
- 3+ variables: `a`, `b`, `c`, ... (skipping `q`)

## Converting Output To SymPy

Tests include an example converter in `tests/test_fk_baseline.py:data_file_to_sympy`.

If you have SymPy installed, the package also provides:

```python
from fkcompute.output.symbolic import matrix_to_polynomial

poly = matrix_to_polynomial(result)
```

## Reformat A Saved Result

If you saved the result as JSON (via `save_data: true`), you can reformat without recomputing:

```bash
fk print-as data/trefoil.json --format latex
```
