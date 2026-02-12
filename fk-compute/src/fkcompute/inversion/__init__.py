"""
Sign assignment (inversion) package for FK computation.

Phase 1 of the FK computation pipeline: finds a valid sign assignment
for a braid at a given degree.

The main entry point is find_sign_assignment() from api.py.
"""

from .api import find_sign_assignment, InversionResult
from .search import (
    BraidType,
    braid_type,
    parallel_try_sign_assignments,
)

__all__ = [
    "find_sign_assignment",
    "InversionResult",
    "BraidType",
    "braid_type",
    "parallel_try_sign_assignments",
]
