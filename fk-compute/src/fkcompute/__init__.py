"""
fkcompute
=========

A package for computing the FK invariant of braids via inversion,
ILP reduction, and a compiled helper binary.

The main interface is a single intelligent `fk()` function that automatically
detects the type of call based on the arguments provided:

1. **Simple mode**: `fk([1,-2,3], 2)`
   - Just braid and degree, uses sensible defaults

2. **Config file mode**: `fk("config.yaml")` or `fk(config="config.yaml")`
   - Loads all parameters from JSON/YAML file
   - Supports batch processing of multiple braids
   - Can include presets and customized parameters in config files

CLI Interface
-------------
The command-line tool `fk` provides two modes:

- `fk simple "[1,-2,3]" 2` - Quick computation with defaults
- `fk config myconfig.yaml` - Configuration file with custom options

Examples
--------
>>> from fkcompute import fk
>>>
>>> # Simple - automatic quiet mode with defaults
>>> result = fk([1, -2, -2, 3], 2)
>>> # Returns: {"braid": [1,-2,-2,3], "inversion_data": {...}, "degree": 2, "fk": {...}}
>>>
>>> # Config file - load from file (supports presets and custom parameters)
>>> result = fk("myconfig.yaml")
>>>
>>> # Access FK coefficients
>>> fk_invariant = result["fk"]
>>> print(fk_invariant)  # {"q^0": 1, "q^2": -1}
>>>
>>> # Access inversion data
>>> inversion = result["inversion_data"]
"""

# Public API imports
from .api.compute import fk
from .api.presets import PRESETS
from .api.batch import fk_from_config, fk_batch_from_config

# Domain types (commonly used)
from .domain.braid.states import BraidStates
from .domain.braid.types import StateLiteral, ZERO_STATE, NUNITY_STATE
from .domain.braid.word import is_positive_braid, is_homogeneous_braid

# Constraint classes
from .domain.constraints.relations import Leq, Less, Zero, Nunity, Alias, Conservation

# Output formatting
from .output.symbolic import format_symbolic_output, SYMPY_AVAILABLE

# CLI entry point
from .cli.app import main

__all__ = [
    # Main API
    "fk",
    "PRESETS",
    "fk_from_config",
    "fk_batch_from_config",
    # Domain types
    "BraidStates",
    "StateLiteral",
    "ZERO_STATE",
    "NUNITY_STATE",
    "is_positive_braid",
    "is_homogeneous_braid",
    # Constraints
    "Leq",
    "Less",
    "Zero",
    "Nunity",
    "Alias",
    "Conservation",
    # Output
    "format_symbolic_output",
    "SYMPY_AVAILABLE",
    # CLI
    "main",
]

__version__ = "0.2.18"
