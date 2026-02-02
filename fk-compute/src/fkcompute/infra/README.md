# Infrastructure Package

This package handles all infrastructure concerns: file I/O, configuration parsing, and external binary execution. It isolates the domain and application layers from system-level details.

## Modules

### `binary.py`
C++ binary execution utilities for running the compiled FK computation engine.

**Functions:**

| Function | Description |
|----------|-------------|
| `binary_path(name)` | Find path to a required binary |
| `run_fk_binary(ilp_file, output_file, threads, verbose)` | Execute FK computation |
| `safe_unlink(path)` | Safely remove a file |

```python
from fkcompute.infra.binary import binary_path, run_fk_binary, safe_unlink

# Find binary location
path = binary_path("fk_main")
print(f"Binary at: {path}")

# Run FK computation
run_fk_binary(
    ilp_file="input_ilp",      # Without extension
    output_file="output",       # Output base name
    threads=4,                  # Thread count
    verbose=True                # Show output
)

# Clean up temp file
safe_unlink("temp_file.csv")
```

**Binary Resolution Order:**
1. Check system PATH
2. Fall back to `fkcompute/_bin/` package directory

### `config.py`
Configuration file parsing for JSON and YAML formats.

**Functions:**

| Function | Description |
|----------|-------------|
| `load_config_file(path)` | Load JSON or YAML config |
| `parse_int_list(string)` | Parse string to int list |
| `load_inversion_file(path)` | Load inversion data JSON |
| `load_ilp_file(path)` | Load ILP data file |

```python
from fkcompute.infra.config import load_config_file, parse_int_list

# Load configuration
config = load_config_file("config.yaml")
print(config["braid"])
print(config["degree"])

# Parse braid string (multiple formats supported)
braid = parse_int_list("[1, -2, 3]")    # JSON style
braid = parse_int_list("1,-2,3")         # Comma separated
braid = parse_int_list("1 -2 3")         # Space separated
```

**Supported Config Formats:**
- `.json` - JSON format
- `.yaml` / `.yml` - YAML format (requires `PyYAML`)

### `io.py`
General file I/O utilities.

**Functions:**

| Function | Description |
|----------|-------------|
| `sort_any(xs)` | Sort list by string representation |
| `find_where(lst, predicate)` | Find first matching element |
| `find_index(lst, predicate)` | Find index of first match |
| `csv_to_dicts(path)` | Read CSV to list of dicts |
| `tsv_to_dicts(path)` | Read TSV to list of dicts |
| `save_dicts_to_tsv(data, path)` | Save dicts to TSV |

```python
from fkcompute.infra.io import sort_any, csv_to_dicts, save_dicts_to_tsv

# Sort heterogeneous list
items = [(1, 2), "hello", 42, (0, 1)]
sorted_items = sort_any(items)

# Read CSV (supports .gz compression)
data = csv_to_dicts("data.csv")
data = csv_to_dicts("data.csv.gz")  # Gzipped

# Save to TSV
save_dicts_to_tsv(data, "output.tsv")
```

## Configuration File Format

### Single Computation

```yaml
# Required
braid: [1, -2, 3]
degree: 2

# Optional
name: my_knot
preset: parallel
max_workers: 4
threads: 2
verbose: true
save_data: true
```

### Batch Processing

```yaml
# Global defaults
preset: single thread
save_data: true

# Multiple computations
computations:
  - name: trefoil
    braid: [1, 1, 1]
    degree: 2

  - name: figure_eight
    braid: [1, -2, 1, -2]
    degree: 3
    preset: parallel  # Override global
```

### Pre-computed Data

```yaml
braid: [1, -2, 3]
degree: 2

# Use pre-computed inversion
inversion:
  0: [1, -1, 1, -1]
  1: [1, 1, -1]

# Or load from file
# inversion_file: path/to/inversion.json

# Or use pre-computed ILP
# ilp_file: path/to/precomputed.csv
```

## Error Handling

All functions in this package handle errors gracefully:

- `binary_path`: Raises `RuntimeError` if binary not found
- `load_config_file`: Raises `FileNotFoundError` or `ImportError` (for YAML)
- `safe_unlink`: Silently ignores errors (safe cleanup)
- File operations: Standard Python exceptions
