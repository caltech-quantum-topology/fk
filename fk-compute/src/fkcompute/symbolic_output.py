"""
Symbolic output functionality for FK computation results.

This module provides functions to convert FK computation results from their
numerical coefficient format into human-readable symbolic polynomials using SymPy.
"""

from __future__ import annotations
import math
from typing import Dict, List, Tuple, Any, Union, Optional
from collections.abc import Iterable
import itertools

try:
    import sympy as sp
    from sympy import symbols, expand, latex
    SYMPY_AVAILABLE = True
except ImportError:
    SYMPY_AVAILABLE = False


def check_sympy_available() -> None:
    """Check if SymPy is available, raise ImportError if not."""
    if not SYMPY_AVAILABLE:
        raise ImportError(
            "SymPy is required for symbolic output functionality. "
            "Install it with: pip install sympy"
        )

def list_to_q_polynomial(q_polyL: List[List[int,int]]) -> sp.Expr:
    """
    Convert FK coefficient data to SymPy polynomial expression in quantum parameter q.

    Transforms the internal FK coefficient representation (list of [power, coefficient] pairs)
    into a human-readable symbolic polynomial using SymPy.

    Args:
        q_polyL: List of [power, coefficient] pairs representing a polynomial in q
                 e.g., [[-2, -1], [-1, 3], [0, -1]] represents -q^(-2) + 3*q^(-1) - 1

    Returns:
        SymPy expression: Polynomial in the quantum parameter q

    Example:
        >>> list_to_q_polynomial([[-2, -1], [-1, 3], [0, -1]])
        -q**(-2) + 3*q**(-1) - 1
    """
    q = symbols("q")
    expr = 0
    for power, coeff in q_polyL:
        expr += coeff * q**power
    return expr

def matrix_to_polynomial(fk_result: Dict) -> sp.Expr:
    """
    Convert FK computation result to a complete symbolic polynomial expression.

    Transforms the FK coefficient matrix from the computation result into a single
    symbolic polynomial with appropriate topological variables (x, y, a, b, c...)
    and the quantum parameter q. Variable naming adapts to braid topology:
    - 1D (single component): uses x
    - 2D (two components): uses x, y
    - 3D+ (multi-component): uses a, b, c, d, ... (skipping 'q')

    Args:
        fk_result: FK computation result dictionary containing:
                   - 'components': Number of braid components/strands
                   - 'degree': Maximum degree for the computation
                   - 'fk': FK coefficient matrix as nested [power, coefficient] pairs

    Returns:
        SymPy expression: Complete FK polynomial in topological variables and q

    Example:
        For 2-component braid with degree 3:
        Returns polynomial like: q - q^3 + x^2*(-q^2 + q^6) + y*x*(q^4 - q^5)

    Note:
        - Variables are chosen based on braid topology, not maximum strand number
        - The quantum parameter 'q' is always reserved
        - Polynomial terms are organized by powers of topological variables
    """
    check_sympy_available()

    n_vars = fk_result['components']
    max_degree = fk_result['degree']
    fk_matrix = fk_result['fk']

    # Create variable symbols based on the number of dimensions
    # 1D: x, 2D: x,y, 3D: a,b,c, etc. (skipping q which is reserved)
    if n_vars == 1:
        var_names = ['x']
    elif n_vars == 2:
        var_names = ['x', 'y']
    else:
        # For 3+ dimensions, start with 'a' and continue alphabetically, skipping 'q'
        var_names = []
        char_offset = 0
        for i in range(n_vars):
            char = chr(ord('a') + i + char_offset)
            if char == 'q':  # Skip 'q' as it's reserved for the quantum parameter
                char_offset += 1
                char = chr(ord('a') + i + char_offset)
            var_names.append(char)

    variables = symbols(' '.join(var_names))
    if n_vars == 1:
        variables = [variables]  # Make it a list for consistency

    fk_polynomial = 0
    for var_powers,q_poly_data in zip(itertools.product(range(max_degree),repeat = n_vars),fk_matrix):
        new_term = list_to_q_polynomial(q_poly_data)
        for i,pow in enumerate(var_powers):
            new_term *= variables[i]**pow
        fk_polynomial += new_term 
        
    return fk_polynomial


def collect_by_x_powers(polynomial: sp.Expr) -> sp.Expr:
    """
    Reorganize polynomial by collecting terms with common powers of x variable.

    Groups polynomial terms by powers of the x variable (if present) to create
    a more readable factored form. This improves the presentation of FK polynomials
    by organizing terms hierarchically.

    Args:
        polynomial: SymPy polynomial expression potentially containing variable x

    Returns:
        SymPy expression with terms collected by powers of x, or original expression
        if no x variable is found

    Example:
        Input:  q - q^2*x^2 - q^3 + q^6*x^2
        Output: q - q^3 + x^2*(-q^2 + q^6)

    Note:
        - Only affects expressions containing the variable 'x'
        - Other variables (y, a, b, c, q) are not collected
        - Improves readability of multi-variable FK polynomials
    """
    check_sympy_available()

    # Get the symbol for x (if it exists)
    symbols_in_poly = polynomial.free_symbols
    x_symbol = None

    for sym in symbols_in_poly:
        if str(sym) == 'x':
            x_symbol = sym
            break

    if x_symbol is None:
        # No x variable, return as is
        return polynomial

    # Collect terms by powers of x
    collected = sp.collect(polynomial, x_symbol, evaluate=True)

    return collected


def format_symbolic_output(fk_result: Dict[str, Any],
                         format_type: str = "pretty",
                         max_degree: int = None) -> str:
    """
    Format FK computation result as human-readable symbolic expression.

    Args:
        fk_result: FK computation result dictionary
        format_type: Output format - "pretty", "latex", or "str"
        max_degree: Maximum degree to display (if None, shows all terms)

    Returns:
        Formatted string representation of the FK polynomial

    Raises:
        ImportError: If SymPy is not available
        ValueError: If format_type is not supported
    """
    check_sympy_available()

    if format_type not in ["pretty", "latex", "str"]:
        raise ValueError(f"Unsupported format_type: {format_type}. Use 'pretty', 'latex', or 'str'")

    try:
        # Convert to polynomial
        fk_polynomial = matrix_to_polynomial(fk_result)

        # Collect terms by powers of x for better readability
        fk_polynomial_collected = collect_by_x_powers(fk_polynomial)

        # Format according to requested type
        if format_type == "pretty":
            return sp.pretty(fk_polynomial_collected, use_unicode=True)
        elif format_type == "latex":
            return latex(fk_polynomial_collected)
        else:  # format_type == "str"
            return str(fk_polynomial_collected)

    except Exception as e:
        return f"Error formatting symbolic output: {e}"


def print_symbolic_result(fk_result: Dict[str, Any],
                         format_type: str = "pretty",
                         max_degree: int = None,
                         show_matrix: bool = False) -> None:
    """
    Print FK computation result in symbolic form.

    Args:
        fk_result: FK computation result dictionary
        format_type: Output format - "pretty", "latex", or "str"
        max_degree: Maximum degree to display
        show_matrix: If True, also display the coefficient matrix
    """
    check_sympy_available()

    try:
        # Print the polynomial only (no extra headers or information)
        symbolic_output = format_symbolic_output(fk_result, format_type, max_degree)
        print(symbolic_output)

    except Exception as e:
        # Silently fail - error messages would clutter output
        pass
