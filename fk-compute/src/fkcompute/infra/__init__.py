"""
Infrastructure layer for FK computation.

This subpackage handles all I/O and subprocess operations:
- binary: C++ binary execution
- config: Config file parsing (JSON/YAML)
- io: File I/O utilities
"""

from .binary import binary_path, run_fk_binary, safe_unlink
from .config import load_config_file, parse_int_list
from .io import sort_any, csv_to_dicts, tsv_to_dicts, save_dicts_to_tsv

__all__ = [
    # Binary
    "binary_path",
    "run_fk_binary",
    "safe_unlink",
    # Config
    "load_config_file",
    "parse_int_list",
    # IO
    "sort_any",
    "csv_to_dicts",
    "tsv_to_dicts",
    "save_dicts_to_tsv",
]
