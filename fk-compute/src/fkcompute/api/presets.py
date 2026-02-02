"""
Preset configurations for FK computation.

This module provides predefined configuration presets for common
use cases in FK computation.
"""

import os
from typing import Any, Dict


def _get_optimal_workers() -> int:
    """Get optimal number of worker processes for this system."""
    cpu_count = os.cpu_count() or 1
    # Reserve 1 core for system responsiveness on multi-core systems
    return max(1, cpu_count - 1) if cpu_count > 2 else cpu_count


def _get_optimal_threads() -> int:
    """Get optimal number of threads for C++ computation."""
    cpu_count = os.cpu_count() or 1
    # Reserve 1 core for system responsiveness on multi-core systems
    return max(1, cpu_count - 1) if cpu_count > 2 else cpu_count


# Preset configurations
PRESETS: Dict[str, Dict[str, Any]] = {
    "single thread": {
        "max_workers": 1,
        "chunk_size": 1 << 12,  # Smaller chunks for faster start
        "include_flip": False,  # Disable flip symmetry
        "max_shifts": None,  # Don't limit shifts
        "verbose": False,
        "save_data": False,
        "threads": 1,
    },
    "parallel": {
        "max_workers": _get_optimal_workers(),  # Auto-detect optimal workers
        "chunk_size": 1 << 14,  # Balanced chunk size
        "include_flip": False,
        "max_shifts": None,
        "verbose": True,
        "save_data": False,
        "threads": _get_optimal_threads(),  # Auto-detect optimal threads
    },
}
