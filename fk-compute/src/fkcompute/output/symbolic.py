"""
Symbolic output functionality for FK computation results.

This module provides functions to convert FK computation results from their
numerical coefficient format into human-readable symbolic polynomials using SymPy.
"""

from __future__ import annotations

from typing import Dict, List, Any, Union

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


def list_to_q_polynomial(q_polyL: List[List[Union[int, str]]]) -> 'sp.Expr':
    """
    Convert FK coefficient data to SymPy polynomial expression in quantum parameter q.

    Parameters
    ----------
    q_polyL
        List of [power, coefficient] pairs representing a polynomial in q.

    Returns
    -------
    sympy.Expr
        Polynomial in the quantum parameter q.
    """
    check_sympy_available()
    q = symbols("q")
    expr = 0
    for power, coeff in q_polyL:
        if isinstance(coeff, str):
            if len(coeff.lstrip('-')) > 10000:
                raise ValueError(
                    f"Coefficient too large for symbolic computation: {len(coeff.lstrip('-'))} digits."
                )
            try:
                coeff = int(coeff)
            except ValueError as e:
                raise ValueError(f"Invalid coefficient string '{coeff}': {e}")
        expr += coeff * q**power
    return expr


def matrix_to_polynomial(fk_result: Dict) -> 'sp.Expr':
    """
    Convert FK computation result to a complete symbolic polynomial expression.

    Parameters
    ----------
    fk_result
        FK computation result dictionary.

    Returns
    -------
    sympy.Expr
        Complete FK polynomial in topological variables and q.
    """
    check_sympy_available()

    n_vars = fk_result['metadata']['num_x_variables']
    terms_data = fk_result["terms"]

    # Create variable symbols based on the number of dimensions
    if n_vars == 1:
        var_names = ['x']
    elif n_vars == 2:
        var_names = ['x', 'y']
    else:
        var_names = []
        char_offset = 0
        for i in range(n_vars):
            char = chr(ord('a') + i + char_offset)
            if char == 'q':
                char_offset += 1
                char = chr(ord('a') + i + char_offset)
            var_names.append(char)

    variables = symbols(' '.join(var_names))
    if n_vars == 1:
        variables = [variables]

    fk_polynomial = 0
    for term in terms_data:
        q_poly_data = [[q_term['q'], q_term['c']] for q_term in term['q_terms']]
        new_term = list_to_q_polynomial(q_poly_data)

        x_powers = term['x']
        for i, pow_val in enumerate(x_powers):
            if i < len(variables):
                new_term *= variables[i]**pow_val

        fk_polynomial += new_term

    return fk_polynomial


def collect_by_variables(polynomial: 'sp.Expr') -> 'sp.Expr':
    """
    Reorganize polynomial by collecting terms with common powers of all
    topological variables (everything except q).

    Parameters
    ----------
    polynomial
        SymPy polynomial expression.

    Returns
    -------
    sympy.Expr
        Polynomial with terms collected by all non-q variables.
    """
    check_sympy_available()

    symbols_in_poly = polynomial.free_symbols
    collect_vars = []
    for sym in sorted(symbols_in_poly, key=str):
        if str(sym) != 'q':
            collect_vars.append(sym)

    if not collect_vars:
        return polynomial

    return sp.collect(polynomial, collect_vars, evaluate=True)


# Backward-compatible alias
collect_by_x_powers = collect_by_variables


def format_symbolic_output(
    fk_result: Dict[str, Any],
    format_type: str = "pretty",
    max_degree: int = None
) -> str:
    """
    Format FK computation result as human-readable symbolic expression.

    Parameters
    ----------
    fk_result
        FK computation result dictionary.
    format_type
        Output format - "pretty", "latex", "inline", "mathematica", or "str".
    max_degree
        Maximum degree to display (if None, shows all terms).

    Returns
    -------
    str
        Formatted string representation of the FK polynomial.
    """
    check_sympy_available()

    valid_formats = ("pretty", "latex", "str", "inline", "mathematica")
    if format_type not in valid_formats:
        raise ValueError(f"Unsupported format_type: {format_type}. Choose from: {', '.join(valid_formats)}")

    try:
        fk_polynomial = matrix_to_polynomial(fk_result)
        fk_polynomial_collected = collect_by_variables(fk_polynomial)

        if format_type == "pretty":
            return sp.pretty(fk_polynomial_collected, use_unicode=True)
        elif format_type == "latex":
            return latex(fk_polynomial_collected)
        elif format_type == "mathematica":
            from sympy.printing.mathematica import mathematica_code
            return mathematica_code(fk_polynomial_collected)
        else:  # "inline" or "str"
            return str(fk_polynomial_collected)

    except Exception as e:
        return f"Error formatting symbolic output: {e}"


def print_symbolic_result(
    fk_result: Dict[str, Any],
    format_type: str = "pretty",
    max_degree: int = None,
    show_matrix: bool = False
) -> None:
    """
    Print FK computation result in symbolic form.

    Parameters
    ----------
    fk_result
        FK computation result dictionary.
    format_type
        Output format - "pretty", "latex", "inline", "mathematica", or "str".
    max_degree
        Maximum degree to display.
    show_matrix
        If True, also display the coefficient matrix.
    """
    check_sympy_available()

    try:
        symbolic_output = format_symbolic_output(fk_result, format_type, max_degree)
        print(symbolic_output)
    except Exception:
        pass
