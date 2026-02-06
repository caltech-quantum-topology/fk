# fk - FK Invariant Computation Tool

A high-performance tool for computing FK invariants for braids and links in knot theory. This package combines Python for orchestration with optimized C++ computation for efficiency.

## Features

- **Multiple Interfaces**: Interactive mode, simple CLI, configuration files, and template generation
- **High Performance**: Optimized C++ computation with parallel processing support
- **Flexible Input**: Multiple braid input formats (JSON, comma-separated, space-separated)
- **Batch Processing**: Process multiple computations in a single run
- **Configuration Files**: YAML and JSON support for reproducible computations
- **Symbolic Output**: Human-readable polynomial expressions using SymPy
- **Template Generation**: Create configuration file templates with documentation

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
pip install .
```

Requirements:
- Python 3.9+
- NumPy >= 1.20
- Gurobi >= 9.5 (optimization solver)
- SymPy >= 1.10 (for symbolic output)
- PyYAML (for YAML configuration files)
- CMake >= 3.20 (for building C++ components)

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
fk
# or explicitly
fk interactive
```

The tool will prompt you for:
- Braid word
- Computation degree
- Number of threads (optional)
- Output preferences

### Simple Usage

For quick computations with defaults:

```bash
# Basic computation
fk simple "[1,-2,3]" 2

# With symbolic output
fk simple "[1,1,1]" 2 --symbolic
```

### Using Templates

Generate a configuration file template:

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

### Configuration Files

#### Single Computation

Create a YAML configuration file:

```yaml
# config.yaml
braid: [1, -2, 3]
degree: 2
name: my_knot
preset: accurate
max_workers: 4
save_data: true
verbose: true
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

computations:
  - name: trefoil
    braid: [1, 1, 1]
    degree: 2

  - name: figure_eight
    braid: [1, -2, 1, -2]
    degree: 3
    preset: accurate  # Override global settings

  - name: hopf_link
    braid: [1, 1]
    degree: 2
    max_workers: 8  # Per-computation override
```

Run batch processing:

```bash
fk config batch.yaml
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

Three built-in presets for common use cases:

- **single thread**: Quick computations with minimal parallelism
- **parallel**: High-performance parallel processing (auto-detects CPU cores)

```yaml
preset: single thread
```

### Available Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `braid` | list[int] | Braid word (required) |
| `degree` | int | Computation degree (required) |
| `name` | str | Computation name (for file naming) |
| `preset` | str | Preset configuration (single thread/parallel) |
| `max_workers` | int | Number of parallel workers |
| `chunk_size` | int | Chunk size for parallel processing |
| `threads` | int | Number of C++ threads |
| `include_flip` | bool | Include flip symmetry |
| `max_shifts` | int | Maximum number of shifts |
| `verbose` | bool | Enable verbose logging |
| `save_data` | bool | Save intermediate data files |
| `symbolic` | bool | Generate symbolic output |
| `ilp_file` | str | Path to precomputed ILP file |
| `inversion` | dict | Precomputed inversion data |

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

```
q - q^3 + x^2*(-q^2 + q^6)
```

Variables are chosen based on topology:
- 1 component: `x`
- 2 components: `x`, `y`
- 3+ components: `a`, `b`, `c`, ... (skipping `q`)

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

1. **Inversion Assignment**: Computes sign assignments for braid crossings using Gurobi ILP solver
2. **ILP Reduction**: Reduces the problem to an Integer Linear Program
3. **C++ Computation**: High-performance computation using FLINT (arbitrary precision), OpenMP (parallelism), and OpenBLAS (linear algebra)

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
