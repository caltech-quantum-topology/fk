"""
Variable assignment logic for FK computation.

This subpackage contains:
- assignment: Variable assignment functions
- symbolic_constraints: Shared symbolic constraint processing
"""

from .assignment import (
    symbolic_variable_assignment,
    extend_variable_assignment,
    equivalence_assignment,
    find_expressions,
    minimal_free,
)

__all__ = [
    "symbolic_variable_assignment",
    "extend_variable_assignment",
    "equivalence_assignment",
    "find_expressions",
    "minimal_free",
]
