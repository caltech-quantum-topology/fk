#kcompute

Compute the **FK invariant** for braids using inversion, ILP reduction, and a compiled helper binary.

This package bundles both Python logic and C++ executables to make the computation portable and easy to use from the command line or within Python code.

---

## Features

### Core Functionality
- **FK Invariant Computation**: Calculate FK invariants for braids using advanced inversion and ILP reduction algorithms
- **Optimized Performance**: Uses compiled C++ helper binary (`fk_segments_links`) for maximum speed
- **Parallel Processing**: Multi-threaded inversion calculations with configurable worker count
- **Multiple Interfaces**: Command-line tool, Python API, and configuration file support

### New Symbolic Output
- **Human-Readable Polynomials**: Convert FK results to symbolic mathematical expressions using SymPy
- **Adaptive Variables**: Automatic variable naming based on braid topology (x for 1D, x,y for 2D, a,b,c... for 3D+)
- **Organized Format**: Terms collected by powers for improved readability

### Command-Line Interface
- **Multiple Usage Modes**:
  - `fk simple` - Quick computations with minimal options
  - `fk preset` - Predefined configurations (fast, accurate, parallel)
  - `fk config` - Configuration file support (JSON/YAML)
  - `fk advanced` - Full parameter control
- **Batch Processing**: Process multiple braids from configuration files
- **Smart Defaults**: Flip symmetry disabled by default for better performance

### Enhanced Output
- **Components Field**: Automatic detection of braid topology components
- **Clean Output**: No debug messages unless verbose mode is enabled
- **Multiple Formats**: JSON output or symbolic polynomial representation
- **Progress Tracking**: Real-time progress for batch computations

### Configuration & Presets
- **Preset Configurations**:
  - `fast` - Single-threaded, optimized for speed
  - `accurate` - Multi-core, comprehensive computation with data saving
  - `parallel` - High parallelism for multi-core systems
- **Flexible Config Files**: JSON/YAML support with global defaults and per-computation overrides
- **Advanced Parameters**: Fine control over ILP, inversion, workers, and more

---

## Installation

You'll need a working C++ toolchain (GCC/Clang on Linux/macOS, MSVC on Windows).
Then install directly:

```bash
pip install .
```

### Man Page Installation

After installation, you can install the man page for the `fk` command:

```bash
# Install man page (attempts user-local first, then system-wide)
fk-install-man
```

Or install manually:

```bash
# Manual installation (requires sudo for system-wide)
sudo cp fk.1 /usr/local/share/man/man1/
sudo mandb
```

Once installed, view the manual with:

```bash
man fk
```

---

## Quick Start

### Simple Usage
```bash
# Basic FK computation
fk simple "[1,-2,3]" 2

# With symbolic polynomial output
fk simple "[1,-2,3]" 2 --symbolic
```

### Preset Usage
```bash
# Fast computation (single-threaded)
fk preset "[1,-2,1,-2]" 3 --preset fast

# Accurate computation (multi-core, saves data)
fk preset "[1,-2,1,-2]" 3 --preset accurate --symbolic

# Parallel optimized computation
fk preset "[1,-2,1,-2]" 3 --preset parallel
```

### Configuration File Usage
Create a configuration file `example.json`:
```json
{
  "braid": [1, -2, 3],
  "degree": 2,
  "preset": "accurate",
  "max_workers": 8,
  "symbolic": true
}
```

Then run:
```bash
fk config example.json
```

### Batch Processing
Create a batch configuration `batch.yaml`:
```yaml
preset: fast
max_workers: 4
computations:
  - name: trefoil
    braid: [1, 1, 1]
    degree: 2
  - name: figure_eight
    braid: [1, -2, 1, -2]
    degree: 3
    preset: accurate
```

Run the batch:
```bash
fk config batch.yaml
```

### Advanced Usage
```bash
# Full parameter control
fk advanced "[1,-2,3]" 2 \
  --max-workers 8 \
  --verbose \
  --save-data \
  --symbolic \
  --save-dir results
```

---

## Output Formats

### Standard JSON Output
```json
{
  "braid": [1, -2, 3],
  "components": 4,
  "degree": 2,
  "fk": [
    [[-1, 1], [0, -1]],
    [[1, -1], [2, 1]]
  ],
  "inversion_data": {...}
}
```

### Symbolic Polynomial Output
When using `--symbolic` flag, you'll see human-readable mathematical expressions:
```
q - q³ + x²(-q² + q⁶)
```

### Braid Input Formats
All these formats are equivalent:
- JSON style: `"[1, -2, -3, 1]"`
- Comma-separated: `"1,-2,-3,1"`
- Space-separated: `"1 -2 -3 1"`

---

## Dependencies

### Required
- Python ≥ 3.9
- NumPy ≥ 1.20
- Gurobi ≥ 9.5 (for ILP solving)
- C++ toolchain for compilation

### Optional
- SymPy ≥ 1.10 (for symbolic output)
- PyYAML (for YAML configuration files)

Install symbolic dependencies:
```bash
pip install sympy pyyaml
```
