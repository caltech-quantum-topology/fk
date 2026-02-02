# API Package

This package provides the main application programming interface for FK computation. It orchestrates the domain layer, solver, and infrastructure to provide a simple, unified interface.

## Modules

### `compute.py`
The core `fk()` function - the main entry point for all FK computations.

**Functions:**

| Function | Description |
|----------|-------------|
| `fk(braid_or_config, degree, ...)` | Unified FK computation |
| `configure_logging(verbose)` | Set logging level |

```python
from fkcompute import fk

# Simple mode - braid and degree
result = fk([1, 1, 1], 2)

# With options
result = fk(
    [1, -2, 3],
    degree=2,
    symbolic=True,      # Add SymPy polynomial
    threads=4,          # C++ threads
    name="my_knot",     # For saved files
    verbose=True,       # Show progress
    save_data=True,     # Keep intermediate files
)

# Config file mode
result = fk("config.yaml")
result = fk(config="config.yaml")
```

**Return Value:**

```python
{
    "terms": [...],           # FK polynomial terms
    "metadata": {
        "braid": [1, 1, 1],
        "inversion": {...},   # Sign assignment
        "components": 1,      # Number of link components
        "num_x_variables": 1,
        "max_x_degrees": [2],
        "symbolic": "...",    # If symbolic=True
    }
}
```

### `batch.py`
Batch processing for multiple FK computations.

**Functions:**

| Function | Description |
|----------|-------------|
| `fk_from_config(path)` | Load and run from config file |
| `fk_batch_from_config(config_data, path)` | Run batch computations |

```python
from fkcompute.api.batch import fk_from_config

# Single computation from config
result = fk_from_config("single.yaml")

# Batch from config with "computations" array
results = fk_from_config("batch.yaml")
# Returns: {"trefoil": {...}, "figure_eight": {...}}
```

### `presets.py`
Predefined configuration presets.

**Available Presets:**

| Preset | Description |
|--------|-------------|
| `single thread` | Single-threaded, fast startup |
| `parallel` | Multi-threaded, auto-detected optimal settings |

```python
from fkcompute import PRESETS

print(PRESETS.keys())
# dict_keys(['single thread', 'parallel'])

# Use preset in config
config = {
    "braid": [1, 1, 1],
    "degree": 2,
    "preset": "parallel"
}
```

**Preset Details:**

```python
# "single thread" preset
{
    "max_workers": 1,
    "chunk_size": 4096,
    "include_flip": False,
    "max_shifts": None,
    "verbose": False,
    "save_data": False,
    "threads": 1,
}

# "parallel" preset (auto-detected)
{
    "max_workers": <cpu_count - 1>,
    "chunk_size": 16384,
    "include_flip": False,
    "max_shifts": None,
    "verbose": True,
    "save_data": False,
    "threads": <cpu_count - 1>,
}
```

## Computation Pipeline

The `fk()` function orchestrates this pipeline:

```
1. Input Parsing
   ├── Config file? → Load and parse
   └── Braid list? → Use directly

2. Sign Assignment (Inversion)
   ├── Homogeneous braid? → Use default signs
   └── Fibered braid? → Parallel search

3. ILP Generation
   └── Build constraint tableau

4. FK Computation
   └── Run C++ binary

5. Result Assembly
   ├── Load JSON output
   ├── Add metadata
   └── Optional: Add symbolic representation
```

## Usage Examples

### Basic Usage

```python
from fkcompute import fk

# Trefoil knot
result = fk([1, 1, 1], 2)
print(result["terms"])
```

### With Symbolic Output

```python
result = fk([1, 1, 1], 2, symbolic=True)
print(result["metadata"]["symbolic"])
# Prints polynomial in x and q
```

### Parallel Computation

```python
result = fk(
    [1, -2, 1, -2],  # Figure-eight (needs inversion search)
    degree=3,
    max_workers=8,   # Parallel inversion search
    threads=4,       # C++ parallelism
    verbose=True,
)
```

### Batch Processing

```python
# batch_config.yaml:
# computations:
#   - name: trefoil
#     braid: [1, 1, 1]
#     degree: 2
#   - name: figure_eight
#     braid: [1, -2, 1, -2]
#     degree: 3

results = fk("batch_config.yaml")
for name, result in results.items():
    print(f"{name}: {len(result['terms'])} terms")
```

### Pre-computed Data

```python
# Skip inversion search with known signs
result = fk(
    [1, -2, 3],
    degree=2,
    inversion={
        "inversion_data": {0: [1, -1, 1, -1, 1, -1]},
        "braid": [1, -2, 3],
        "degree": 2,
    }
)

# Skip to FK computation with ILP file
result = fk(
    [1, -2, 3],
    degree=2,
    ilp_file="precomputed.csv"
)
```

## Error Handling

The API provides clear error messages:

```python
# Missing degree
fk([1, 1, 1])
# ValueError: degree is required when providing a braid list

# Invalid config
fk("nonexistent.yaml")
# FileNotFoundError: Config file not found: nonexistent.yaml

# Inversion failure
result = fk([complex_braid], 2)
if result.get("inversion_data") == "failure":
    print("No valid sign assignment found")
```
