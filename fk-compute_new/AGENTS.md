# AGENTS.md

This file contains guidelines and commands for agentic coding agents working in this repository.

## Project Overview

`fkcompute` is a high-performance Python package for computing FK invariants for braids and links in knot theory. The project combines Python orchestration with optimized C++ computation using FLINT, OpenMP, and OpenBLAS.

**Key Components:**
- Python package (`src/fkcompute/`) for high-level logic and user interfaces
- C++ components (`cpp/`) for performance-critical computation
- CLI tool with interactive and batch processing modes
- Mathematica integration via Wolfram Language wrapper

## Build Commands

### Python Package Installation
```bash
# Basic installation
pip install .

# Installation with optional dependencies
pip install .[full]          # All optional features
pip install .[symbolic]       # Symbolic output only
pip install .[interactive]    # Enhanced interactive mode only
pip install .[yaml]           # YAML configuration files only

# Development installation
pip install -e .

# Build C++ components manually (if needed)
cmake -B build -S .
cmake --build build

# Install manual page
fk-install-man
```

### System Dependencies
Before installation, ensure these are installed:
- **macOS**: `brew install flint libomp openblas`
- **Ubuntu**: `sudo apt-get install libflint-dev libgomp1-dev libopenblas-dev`
- **RHEL**: `sudo yum install flint-devel gcc-openmp openblas-devel`

## Testing Commands

### Python Tests
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_fk_baseline.py

# Run single test with verbose output
pytest tests/test_fk_baseline.py::test_fk_matches_baseline[3_1] -v

# Run with coverage
pytest --cov=fkcompute

# Run tests in parallel
pytest -n auto

# Run only fast tests (add markers if available)
pytest -m "not slow"
```

### C++ Tests
```bash
# Build and run C++ tests
cd build
ctest

# Verbose C++ test output
ctest --verbose
```

### Manual Test Execution
```bash
# Run baseline comparison test directly
python tests/test_fk_baseline.py

# Run CLI tests
fk simple "[1,1,1]" 2 --symbolic
```

## Code Style Guidelines

### Python Code Style

**Type Hints:**
- Use type hints consistently (`from typing import List, Dict, Optional, Union`)
- Modern type syntax: `list[int]` instead of `List[int]` for Python 3.9+
- Import from `typing` only when needed for complex types

**Import Organization:**
```python
# 1. Standard library imports
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# 2. Third-party imports
import numpy as np
import pytest
import typer
import sympy

# 3. Local imports (relative imports preferred)
from .presets import PRESETS
from .batch import fk_from_config
from ..domain.braid.states import BraidStates
from ..inversion.api import find_sign_assignment
```

**Docstrings:**
- Use triple quotes with descriptive summary
- Follow Google/NumPy style for complex functions
- Include parameter types and return values
- Example:
```python
def compute_fk(knot_name: str, config: dict) -> dict:
    """Compute FK polynomial for a knot.
    
    Args:
        knot_name: Name identifier for the knot.
        config: Configuration dictionary containing braid and parameters.
        
    Returns:
        Dictionary containing FK computation results.
    """
```

**Naming Conventions:**
- Functions: `snake_case`
- Variables: `snake_case` (descriptive names)
- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Private members: prefix with `_`

**Error Handling:**
- Use specific exceptions (`ValueError`, `FileNotFoundError`, etc.)
- Log errors with appropriate levels
- Provide meaningful error messages
- Use `try/except` blocks for external dependencies

**Logging:**
```python
import logging
logger = logging.getLogger("fk_logger")

def configure_logging(verbose: bool) -> None:
    """Set logging level based on verbose flag."""
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)
```

### C++ Code Style

**Headers and Includes:**
- Use appropriate standard headers (`<iostream>`, `<vector>`, etc.)
- Group includes: standard, third-party, local
- Use include guards or `#pragma once`

**Naming:**
- Functions: `snake_case`
- Variables: `snake_case`
- Classes: `PascalCase`
- Constants: `kCamelCase`

## File Organization

### Directory Structure
```
src/fkcompute/
├── __init__.py           # Public API exports
├── api/                  # Core API functions
│   ├── compute.py        # Main fk() function
│   ├── batch.py          # Batch processing
│   └── presets.py        # Configuration presets
├── domain/               # Domain models and types
│   ├── braid/           # Braid-related classes
│   └── constraints/     # Constraint definitions
├── inversion/           # Sign assignment computation
├── solver/             # ILP solver interface
├── cli/                # Command-line interface
├── infra/              # Infrastructure utilities
└── output/             # Output formatting

tests/
├── test_fk_baseline.py  # Main regression tests
├── test_phase2_snapshots.py
├── data/               # Test data files
└── data_baseline/      # Baseline comparison data

cpp/                   # C++ computation components
├── src/              # C++ source files
├── include/          # C++ headers
└── Makefile          # C++ build configuration

mathematica/          # Wolfram Language wrapper
└── FkCompute/       # Mathematica paclet
```

## Development Workflow

### Making Changes
1. **Understand the impact**: Check if changes affect Python API, C++ components, or CLI
2. **Write tests**: Add appropriate test cases before implementing
3. **Implement changes**: Follow existing code patterns and style
4. **Test thoroughly**: Run both Python and C++ tests
5. **Check imports**: Ensure no unused imports remain

### Performance Considerations
- C++ components are performance-critical; profile before optimizing
- Use parallel processing (`max_workers`, `threads`) for heavy computations
- Consider memory vs speed tradeoffs with `chunk_size` parameter
- Precomputed data can skip expensive computation steps

### Dependencies Management
- **Core requirements**: numpy>=1.20, gurobipy>=9.5, typer>=0.9.0
- **Optional dependencies**: sympy>=1.10 (symbolic), rich>=13 (interactive), PyYAML>=6 (yaml)
- **System dependencies**: FLINT, OpenMP, OpenBLAS (handled by installation docs)

## CLI Interface Guidelines

### Command Structure
```bash
# Main entry point: fk
fk --help

# Interactive modes
fk                    # Enhanced interactive (default)
fk interactive        # Explicit interactive
fk interactive --quick  # Quick prompts

# Simple computation
fk simple "[1,-2,3]" 2

# Configuration file
fk config config.yaml
fk template create    # Generate template

# History management (when available)
fk history show
```

### Adding New CLI Commands
- Use Typer for CLI definitions
- Add command functions in `src/fkcompute/cli/commands.py`
- Include help text and parameter validation
- Consider both interactive and programmatic usage

## Testing Strategy

### Test Categories
1. **Unit tests**: Individual functions and classes
2. **Integration tests**: End-to-end FK computation
3. **Baseline tests**: Compare against known results
4. **CLI tests**: Command-line interface functionality
5. **Performance tests**: Benchmark critical computations

### Test Data
- Store test data in `tests/data/` and `tests/data_baseline/`
- Use descriptive naming: `knot_name_degree.json`
- Include both simple and complex braid examples
- Maintain baseline data for regression testing

## Common Patterns

### API Function Structure
```python
def fk(
    braid_or_config: Union[List[int], str],
    *args,
    **kwargs
) -> Dict[str, Any]:
    """Main FK computation function."""
    # 1. Parse and validate inputs
    # 2. Configure logging
    # 3. Execute computation pipeline
    # 4. Format and return results
```

### Configuration Processing
```python
def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML/JSON file."""
    # 1. Validate file existence
    # 2. Parse file format
    # 3. Apply presets if specified
    # 4. Validate required parameters
```

### Error Handling Pattern
```python
try:
    result = risky_operation()
except SpecificError as e:
    logger.error(f"Operation failed: {e}")
    raise ValueError(f"Cannot compute FK: {e}") from e
```

## Important Notes

1. **Gurobi Dependency**: The project requires Gurobi optimization solver. Ensure proper licensing and installation.
2. **Cross-Platform**: C++ components must work on macOS, Linux, and handle different OpenMP availability.
3. **Memory Management**: Large computations can be memory-intensive; monitor memory usage in tests.
4. **Backward Compatibility**: Maintain compatibility with existing configuration files and API calls.
5. **Mathematica Integration**: Changes to Python API may require updates to the Wolfram Language wrapper.

## Linting and Formatting

While no specific linting configuration was found, follow Python best practices:
- Use `black` or similar for consistent formatting
- Use `flake8` or `ruff` for linting
- Ensure type hints are complete and accurate
- Keep docstrings consistent with existing style

Before committing changes:
```bash
# If linting tools are configured
ruff check src/
black src/

# Run full test suite
pytest
```