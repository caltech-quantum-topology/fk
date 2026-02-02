"""
Sign assignment logic for FK computation.

This subpackage contains:
- assignment: Variable assignment functions
- validation: Sign assignment validation and degree bounds
"""

from .assignment import (
    symbolic_variable_assignment,
    extend_variable_assignment,
    equivalence_assignment,
    find_expressions,
    minimal_free,
)
from .validation import (
    check_sign_assignment,
    czech_sign_assignment,
    violates_relation,
    violates_any_relation,
)

__all__ = [
    # Assignment
    "symbolic_variable_assignment",
    "extend_variable_assignment",
    "equivalence_assignment",
    "find_expressions",
    "minimal_free",
    # Validation
    "check_sign_assignment",
    "czech_sign_assignment",
    "violates_relation",
    "violates_any_relation",
]
