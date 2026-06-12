# Quickstart

## CLI

Compute an FK invariant quickly:

```bash
fk simple "[1,1,1]" 2
```

Symbolic output (requires `sympy` / `fkcompute[symbolic]`):

```bash
fk simple "[1,1,1]" 2 --symbolic
fk simple "[1,1,1]" 2 --format latex
fk simple "[1,1,1]" 2 --format mathematica
```

Interactive wizard:

```bash
fk
fk interactive
fk interactive --quick
```

Run from a config file:

```bash
fk config my_run.yaml
```

## Python

The public API is `fkcompute.fk`:

```python
from fkcompute import fk

result = fk([1, 1, 1], 2)
print(result["metadata"]["components"])
print(len(result["terms"]))
```

With symbolic formatting:

```python
result = fk([1, 1, 1], 2, symbolic=True)
print(result["metadata"].get("symbolic"))
```

## Save Data + Reformat Later

If you enable `save_data`, the computation writes intermediate files:

- `<name>_inversion.json`
- `<name>_ilp.csv`
- `<name>.json`

Example:

```python
from fkcompute import fk

fk([1, 1, 1], 2, name="trefoil", save_data=True, save_dir="data")
```

Reformat an existing JSON result file (no recomputation):

```bash
fk print-as data/trefoil.json --format inline
```

## Minimal Config Example

`my_run.yaml`:

```yaml
braid: [1, 1, 1]
degree: 2
symbolic: true
preset: parallel
threads: 4
save_data: true
save_dir: data
name: trefoil_d2
```

Run:

```bash
fk config my_run.yaml
```
