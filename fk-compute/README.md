# fk - FK Invariant Computation Tool

A high-performance tool for computing FK invariants for braids and links in knot theory. This package combines Python for orchestration with optimized C++ computation for efficiency.

## Features

- **Multiple Interfaces**: Interactive mode (basic, enhanced, quick), simple CLI, configuration files, and template generation
- **High Performance**: Optimized C++ computation with parallel processing support (Python workers + C++ threads)
- **Flexible Input**: Multiple braid input formats (JSON, comma-separated, space-separated)
- **Batch Processing**: Process multiple computations in a single run with progress tracking
- **Configuration Files**: YAML and JSON support for reproducible computations with presets
- **Symbolic Output**: Human-readable polynomial expressions using SymPy (multiple formats: pretty, inline, LaTeX, Mathematica)
- **Template Generation**: Create configuration file templates with comprehensive documentation
- **History Management**: Track, search, and reuse previous computations (when interactive dependencies available)
- **Mathematica Integration**: Complete Wolfram Language wrapper paclet for seamless Mathematica usage
- **Precomputed Data Support**: Load/save inversion and ILP data to skip expensive computation steps

## Installation

### System Requirements

Before installing, ensure you have the following system dependencies:

**macOS:**
```bash
brew install flint libomp openblas
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libflint-dev libgomp1-dev libopenblas-dev
```

**RHEL/Fedora:**
```bash
sudo yum install flint-devel gcc-openmp openblas-devel
```

### Python Package

Install from source:

```bash
# Basic installation
pip install .

# Installation with optional dependencies
pip install .[full]  # All optional features
pip install .[symbolic]  # Symbolic output only
pip install .[interactive]  # Enhanced interactive mode only
pip install .[yaml]  # YAML configuration files only
```

**Core Requirements:**
- Python 3.9+
- NumPy >= 1.20
- Gurobi >= 9.5 (optimization solver)
- Typer >= 0.9.0 (CLI framework)
- CMake >= 3.20 (for building C++ components)

**Optional Requirements:**
- SymPy >= 1.10 (for symbolic output)
- Rich >= 13 (for enhanced interactive mode)
- PyYAML >= 6 (for YAML configuration files)
- IPython >= 8, Jupyter (for visualization and development)

### Post-Installation

Install the manual page:

```bash
fk-install-man
```

Then access documentation with:

```bash
man fk
```

## Mathematica (Wolfram Language) Wrapper: `FkCompute`

This repository includes a comprehensive Wolfram Language wrapper (Paclet) that provides seamless integration between Mathematica and the `fkcompute` Python package. The wrapper supports all core functionality including symbolic computation, batch processing, and configuration file execution.

### Prerequisites

- Mathematica 12.0 or later
- Python 3.9+ with `fkcompute` package installed (see Python installation above)
- SymPy 1.10+ (for symbolic output functionality)

### Installation

**1. Install the Python package:**
```bash
pip install .
```

**2. Install the Mathematica paclet:**
```wl
PacletInstall[Directory["/path/to/fk-compute/mathematica/FkCompute"]]
```

Replace `/path/to/fk-compute` with the absolute path to this repository on your system.

This permanently installs the paclet, enabling `Needs["FkCompute"]` in any future Mathematica session.

### Load and Test

```wl
Needs["FkCompute`"]
FkComputeVersion[]
```

If successful, `FkComputeVersion[]` returns version information for both Python and the `fkcompute` package.

### Usage Examples

**Basic Computation:**
```wl
res = FkCompute[{1, 1, 1}, 2];
res["metadata", "symbolic"]
```

**Symbolic Computation:**
```wl
res = FkCompute[{1, 1, 1}, 2, "Symbolic" -> True, "Threads" -> 4];
res["metadata", "symbolic"]
(* Output: q - q^3 + x^2*(-q^2 + q^6) *)
```

**Configuration File Execution:**
```wl
res = FkCompute["config.yaml"];
```

**Advanced Options:**
```wl
res = FkCompute[{1, -2, 1, -2}, 3, 
  "Symbolic" -> True,
  "Threads" -> 8,
  "Verbose" -> True,
  "MaxWorkers" -> 4,
  "Preset" -> "parallel"
];
```

### Available Options

The `FkCompute` function supports all parameters available in the Python API:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `"PythonExecutable"` | String | `Automatic` | Path to Python executable |
| `"Symbolic"` | Boolean | `False` | Generate symbolic polynomial output |
| `"Threads"` | Integer | `Automatic` | Number of C++ computation threads |
| `"Verbose"` | Boolean | `False` | Enable verbose logging |
| `"MaxWorkers"` | Integer | `Automatic` | Maximum parallel workers |
| `"Preset"` | String | `Automatic` | Configuration preset (`"single thread"`, `"parallel"`) |
| `"Name"` | String | `Automatic` | Computation name for file naming |
| `"SaveData"` | Boolean | `Automatic` | Save intermediate computation data |

### Python Configuration

**Set Default Python Executable:**
```wl
$FkComputePythonExecutable = "/usr/bin/python3";
```

**Per-Call Python Specification:**
```wl
FkCompute[{1, 1, 1}, 2, "PythonExecutable" -> "/opt/conda/bin/python"]
```

**Helper Functions:**
```wl
SetFkComputePythonExecutable["/usr/bin/python3"]
GetFkComputePythonExecutable[]
```

### Development (Optional)

For development work, load the paclet directly from the source directory:

```wl
PacletDirectoryLoad["/path/to/fk-compute/mathematica"]
Needs["FkCompute`"]
```

Note: This only loads the paclet for the current Mathematica session and doesn't require installation.

### Error Handling

The wrapper provides comprehensive error messages:

- `FkCompute::pyfail`: Python execution failed
- `FkCompute::pybadjson`: Invalid JSON output from Python
- `FkCompute::pyerror`: fkcompute Python package errors

### Integration with Mathematica

The wrapper returns Mathematica `Association` objects, enabling natural integration:

```wl
res = FkCompute[{1, 1, 1}, 2, "Symbolic" -> True];
symbolicForm = res["metadata", "symbolic"];
components = res["metadata", "components"];
braid = res["metadata", "braid"];
```

## Quick Start

### Interactive Mode

Start an interactive session with guided prompts:

```bash
# Enhanced interactive mode (default, with progress tracking)
fk
fk interactive

# Quick interactive mode (minimal prompts)
fk interactive --quick

# Enhanced mode explicitly
fk interactive --enhanced
```

The tool will prompt you for:
- Braid word (supports multiple formats)
- Computation degree (positive integer)
- Number of threads (optional, auto-detected)
- Computation name (optional, for file naming)
- Symbolic output preference and format selection
- Data saving preferences

**Interactive Features:**
- Progress tracking during computation
- History management (search, reuse previous computations)
- Input validation and helpful error messages
- Format selection for symbolic output (pretty, inline, LaTeX, Mathematica)

### Simple Usage

For quick computations with sensible defaults:

```bash
# Basic computation
fk simple "[1,-2,3]" 2

# With symbolic output (auto-enables SymPy formatting)
fk simple "[1,1,1]" 2 --symbolic

# Different symbolic formats
fk simple "[1,1,1]" 2 --format inline
fk simple "[1,1,1]" 2 --format latex
fk simple "[1,1,1]" 2 --format mathematica

# Custom thread count
fk simple "[1,-2,3]" 2 --threads 4
```

### Using Templates

Generate a configuration file template with all options documented:

```bash
# Create template with default name (fk_config.yaml)
fk template create

# Create template with custom name
fk template create my_computation.yaml

# Overwrite existing file
fk template create config.yaml --overwrite
```

Edit the generated template and run:

```bash
fk config my_computation.yaml
```

**Template Features:**
- Comprehensive documentation for all parameters
- Example configurations for common use cases
- Preset usage examples
- Precomputed data format specifications

### Configuration Files

#### Single Computation

Create a YAML configuration file:

```yaml
# config.yaml
braid: [1, -2, 3]
degree: 2
name: my_knot
preset: parallel
max_workers: 4
save_data: true
verbose: true
symbolic: true
```

Run the computation:

```bash
fk config config.yaml
```

#### Multiple Configuration Files

Process multiple configurations in sequence:

```bash
fk config config1.yaml config2.yaml config3.json
```

The tool will:
- Process each file sequentially
- Show progress for each computation
- Display individual results
- Provide a summary of successful/failed computations

#### Batch Processing

Process multiple braids with shared settings:

```yaml
# batch.yaml
max_workers: 4
save_data: true
verbose: false
symbolic: true

computations:
  - name: trefoil
    braid: [1, 1, 1]
    degree: 2

  - name: figure_eight
    braid: [1, -2, 1, -2]
    degree: 3
    preset: parallel  # Override global settings

  - name: hopf_link
    braid: [1, 1]
    degree: 2
    max_workers: 8  # Per-computation override
    symbolic: false  # Override global symbolic setting
```

Run batch processing:

```bash
fk config batch.yaml
```

#### Advanced Configuration Examples

**Using Precomputed Data:**
```yaml
braid: [1, -2, 1, -2]
degree: 3
# Skip expensive inversion computation
inversion:
  0: [1, -1, 1, -1]
  1: [-1, 1, 1, -1]
# Load precomputed ILP
ilp_file: "precomputed/ilp_data.csv"
```

**Performance Optimization:**
```yaml
braid: [1, 1, 1, 1, 1]
degree: 10
preset: parallel
# Fine-tune performance
max_workers: 8
threads: 16
chunk_size: 65536
# Limit search space for faster computation
max_shifts: 1000
```

## Braid Input Formats

Braids can be specified in multiple formats:

```bash
# JSON-style (with brackets and commas)
fk simple "[1, -2, 3, 1]" 2

# Comma-separated
fk simple "1,-2,3,1" 2

# Space-separated
fk simple "1 -2 3 1" 2
```

## Configuration Options

### Presets

Two built-in presets for common use cases:

- **single thread**: Quick computations with minimal parallelism (1 worker, 1 thread, smaller chunks)
- **parallel**: High-performance parallel processing (auto-detects optimal CPU cores, larger chunks)

```yaml
preset: parallel
```

### Complete Configuration Parameters

This section documents all available configuration options for FK computation.

#### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `braid` | list[int] | Braid word as list of integers (e.g., `[1, -2, 3]`) |
| `degree` | int | Computation degree (positive integer) |

#### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| **Computation Settings** | | | |
| `name` | str | `auto-generated` | Computation name (used for output files when `save_data` is true) |
| `preset` | str | `none` | Use preset configuration (`"single thread"` or `"parallel"`) |
| `symbolic` | bool | `false` | Generate symbolic polynomial output using SymPy |
| **Performance Settings** | | | |
| `max_workers` | int | `1` | Maximum number of worker processes for parallel computation |
| `chunk_size` | int | `16384` | Chunk size for parallel processing (power of 2) |
| `threads` | int | `1` | Number of threads for C++ FK computation |
| `include_flip` | bool | `false` | Include flip symmetry in computation |
| `max_shifts` | int | `null` | Maximum number of shifts (null for unlimited) |
| **I/O Settings** | | | |
| `verbose` | bool | `false` | Enable verbose logging with debug information |
| `save_data` | bool | `false` | Save intermediate computation data to files |
| `save_dir` | str | `"data"` | Directory for saved data files |
| **Precomputed Data** | | | |
| `inversion` | dict | `null` | Precomputed inversion data (format: `{crossing_index: [signs]}`) |
| `inversion_file` | str | `null` | Path to JSON file containing precomputed inversion data |
| `ilp_file` | str | `null` | Path to precomputed ILP problem file |
| `partial_signs` | list[int] | `null` | Partial sign assignments to constrain computation |

#### Parameter Details

**Performance Optimization:**
- `max_workers`: Controls Python-level parallelism for inversion computation
- `threads`: Controls C++-level parallelism for FK binary computation  
- `chunk_size`: Larger chunks reduce overhead but use more memory
- Presets automatically configure optimal values for your system

**Precomputed Data Usage:**
```yaml
# Precomputed inversion data (skip expensive inversion step)
inversion:
  0: [1, -1, 1, -1]
  1: [-1, 1, 1, -1]

# Load from external file
inversion_file: "path/to/inversion.json"
ilp_file: "path/to/ilp.csv"
```

**Advanced Options:**
```yaml
# Constrain computation with partial signs
partial_signs: [1, -1, 0, 1]  # 0 = unconstrained

# Limit search space
max_shifts: 1000  # Reduce computation time for complex braids

# Memory vs speed tradeoff
chunk_size: 65536  # Larger = faster, more memory
```

## Output

### JSON Output

Default output is JSON format:

```json
{
  "braid": [1, -2, 3],
  "inversion_data": {
    "0": [1, -1, 1],
    "1": [-1, 1, 1]
  },
  "degree": 2,
  "components": 3,
  "fk": [
    [[-2, 1], [0, -1], [2, 1]]
  ]
}
```

### Symbolic Output

With `--symbolic` flag or `symbolic: true` in config:

```bash
# Default pretty format
fk simple "[1,1,1]" 2 --symbolic
# Output: q - q^3 + x^2*(-q^2 + q^6)

# Inline format
fk simple "[1,1,1]" 2 --format inline
# Output: q - q^3 + x^2*(-q^2 + q^6)

# LaTeX format
fk simple "[1,1,1]" 2 --format latex
# Output: $q - q^{3} + x^{2} \left(-q^{2} + q^{6}\right)$

# Mathematica format
fk simple "[1,1,1]" 2 --format mathematica
# Output: q - q^3 + x^2*(-q^2 + q^6)
```

**Variable Assignment:**
Variables are chosen based on braid topology:
- 1 component: `x`
- 2 components: `x`, `y`
- 3+ components: `a`, `b`, `c`, ... (skipping `q`)

**Requirements:**
- SymPy >= 1.10 must be installed (`pip install sympy`)
- Install with optional dependencies: `pip install fkcompute[symbolic]`

## Examples

### Classic Knots

**Trefoil Knot:**
```bash
fk simple "[1,1,1]" 2 --symbolic
```

**Figure-Eight Knot:**
```bash
fk simple "[1,-2,1,-2]" 3 --symbolic
```

**Hopf Link:**
```bash
fk simple "[1,1]" 2
```

### Workflow Example

```bash
# 1. Create a configuration template
fk template create trefoil.yaml

# 2. Edit trefoil.yaml to set your parameters
# (set braid: [1,1,1], degree: 2, etc.)

# 3. Run the computation
fk config trefoil.yaml

# 4. Process multiple configurations
fk config trefoil.yaml figure_eight.yaml hopf.yaml
```

### Batch Computation Workflow

```bash
# Create batch configuration
cat > knots.yaml << 'EOF'
max_workers: 4
save_data: true
verbose: true

computations:
  - name: trefoil_d2
    braid: [1, 1, 1]
    degree: 2
  - name: trefoil_d3
    braid: [1, 1, 1]
    degree: 3
  - name: figure_eight_d3
    braid: [1, -2, 1, -2]
    degree: 3
EOF

# Run batch
fk config knots.yaml
```

## Architecture

The tool uses a multi-stage pipeline:

1. **Inversion Assignment**: Computes sign assignments for braid crossings using Gurobi ILP solver with parallel processing
2. **ILP Reduction**: Reduces the problem to an Integer Linear Program for efficient solving
3. **C++ Computation**: High-performance computation using FLINT (arbitrary precision), OpenMP (parallelism), and OpenBLAS (linear algebra)
4. **Symbolic Processing**: Optional conversion to human-readable polynomial expressions using SymPy
5. **Output Formatting**: Multiple output formats (JSON, symbolic pretty-print, LaTeX, Mathematica)

**Performance Optimizations:**
- Parallel Python workers for inversion computation
- Multi-threaded C++ binary for FK computation
- Configurable chunk sizes for memory vs speed tradeoffs
- Precomputed data support to skip expensive steps
- Preset configurations for optimal system performance

## Development

### Building from Source

```bash
# Clone repository
git clone <repository-url>
cd fk-compute

# Install in development mode
pip install -e .

# Build C++ components
cmake -B build -S .
cmake --build build
```

### Running Tests

```bash
# Python tests
pytest

# C++ tests
cd build
ctest
```

### History Management

When interactive dependencies are available, use computation history:

```bash
# Show recent computations
fk history show

# Search history
fk history search "trefoil"

# Clear history
fk history clear --confirm

# Export/import history
fk history export my_history.json
fk history import my_history.json
```

### Python API

Direct Python usage:

```python
from fkcompute import fk

# Simple computation
result = fk([1, -2, 3], 2)

# With options
result = fk([1, -2, 3], 2, symbolic=True, threads=4, verbose=True)

# From config file
result = fk("config.yaml")

# Access results
fk_invariant = result["fk"]
inversion_data = result["inversion_data"]
```

## License

MIT License - see LICENSE file for details

## Authors

Toby, Paul, Lara, Josef, Davide

## Citations

If you use this tool in your research, please cite:

```bibtex
@software{fkcompute,
  title={fk: FK Invariant Computation Tool},
  author={Toby and Paul and Lara and Josef and Davide},
  year={2025},
  url={https://github.com/yourusername/fk-compute}
}
```

## Support

For issues, questions, or contributions:
- Manual: `man fk`
- Help: `fk --help`, `fk <command> --help`
- GitHub Issues: [Report a bug or request a feature]

## See Also

- Gurobi Optimization: https://www.gurobi.com/
- FLINT Library: https://flintlib.org/
- SymPy Documentation: https://docs.sympy.org/
