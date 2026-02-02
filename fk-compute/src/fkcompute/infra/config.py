"""
Configuration file parsing for FK computation.

This module provides functions for loading and parsing configuration
files in JSON and YAML formats.
"""

import json
import os
from typing import Any, Dict, List, Optional

# Optional YAML support
HAS_YAML = False
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def load_config_file(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from JSON or YAML file.

    Parameters
    ----------
    config_path
        Path to the configuration file.

    Returns
    -------
    dict
        Parsed configuration data.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    ImportError
        If YAML file is requested but PyYAML is not installed.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        if config_path.endswith((".yml", ".yaml")):
            if not HAS_YAML:
                raise ImportError(
                    "PyYAML is required for YAML config files. Install with: pip install PyYAML"
                )
            return yaml.safe_load(f)
        else:
            return json.load(f)


def parse_int_list(s: Optional[str]) -> Optional[List[int]]:
    """
    Parse a string into a list of ints.

    Accepts:
    - JSON-style: "[1, -2, 3]"
    - Comma-separated: "1,-2,3"
    - Space-separated: "1 -2 3"

    Parameters
    ----------
    s
        String to parse.

    Returns
    -------
    list[int] or None
        Parsed list of integers, or None if input is None/empty.
    """
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    # Try JSON first
    try:
        val = json.loads(s)
        if isinstance(val, list):
            return [int(x) for x in val]
    except Exception:
        pass
    # Fallback: split on commas/spaces
    parts = [p for p in s.replace(",", " ").split() if p]
    return [int(p) for p in parts]


def load_inversion_file(inversion_file: str) -> Dict[str, Any]:
    """
    Load inversion data from a JSON file and ensure integer keys.

    Parameters
    ----------
    inversion_file
        Path to the inversion data file.

    Returns
    -------
    dict
        Inversion data with integer keys.
    """
    with open(inversion_file, "r") as f:
        inversion = json.load(f)
    inversion["inversion_data"] = {
        int(k): v for k, v in inversion["inversion_data"].items()
    }
    return inversion


def load_ilp_file(ilp_file: str) -> str:
    """
    Load ILP data from a file.

    Parameters
    ----------
    ilp_file
        Path to the ILP data file.

    Returns
    -------
    str
        ILP data as a string.
    """
    with open(ilp_file, "r") as f:
        return f.read()
