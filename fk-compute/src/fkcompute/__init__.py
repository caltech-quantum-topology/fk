"""
fk-compute
============

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

from .fk import fk, PRESETS
from .cli import main

__all__ = ["fk", "PRESETS", "main"]

__version__ = "0.1.0"
