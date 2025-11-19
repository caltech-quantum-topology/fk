#!/usr/bin/env python3
"""
Script to read JSON files from results folder and format them as nice sympy expressions.
The JSON files represent series in x with coefficients that are Laurent polynomials in q.
"""

import json
import sys
import argparse
from pathlib import Path
from sympy import symbols, sympify, latex, pretty
from sympy.core.add import Add
from sympy.core.mul import Mul


def parse_json_to_sympy(json_data, ascending=True):
    """
    Convert JSON data to a sympy expression.

    Args:
        json_data: Dictionary containing 'terms' with x powers and q Laurent polynomial coefficients
        ascending: If True, sort by ascending x powers; if False, sort by descending x powers

    Returns:
        sympy expression representing the series
    """
    x, q = symbols('x q')

    expression_terms = []

    # Sort terms by x power
    sorted_terms = sorted(json_data['terms'], key=lambda term: term['x'][0], reverse=not ascending)

    for term in sorted_terms:
        # Get the x power (assuming single variable for now)
        x_power = term['x'][0]

        # Build the q Laurent polynomial coefficient
        q_poly_terms = []
        for q_term in term['q_terms']:
            q_power = q_term['q']
            coefficient = q_term['c']

            if coefficient != 0:
                if q_power == 0:
                    q_poly_terms.append(coefficient)
                else:
                    q_poly_terms.append(coefficient * (q ** q_power))

        # Sum up all q terms for this x power
        if q_poly_terms:
            q_polynomial = sum(q_poly_terms)

            # Create the full term: (q_polynomial) * x^power
            if x_power == 0:
                full_term = q_polynomial
            elif x_power == 1:
                full_term = q_polynomial * x
            else:
                full_term = q_polynomial * (x ** x_power)

            expression_terms.append(full_term)

    # Sum all terms to get the final expression
    if expression_terms:
        return sum(expression_terms)
    else:
        return sympify(0)


def format_expression(expr, format_type='pretty'):
    """
    Format the sympy expression in the specified format.

    Args:
        expr: sympy expression
        format_type: 'pretty', 'latex', or 'str'

    Returns:
        formatted string
    """
    if format_type == 'pretty':
        return pretty(expr, use_unicode=True)
    elif format_type == 'latex':
        return latex(expr)
    elif format_type == 'str':
        return str(expr)
    else:
        return str(expr)


def main():
    parser = argparse.ArgumentParser(description='Format FK computation results as sympy expressions')
    parser.add_argument('json_file', help='Path to JSON file to format')
    parser.add_argument('--format', '-f', choices=['pretty', 'latex', 'str'],
                       default='pretty', help='Output format (default: pretty)')
    parser.add_argument('--metadata', '-m', action='store_true',
                       help='Also display metadata from the JSON file')
    parser.add_argument('--order', '-o', choices=['ascending', 'descending'],
                       default='ascending', help='Order of x powers in output (default: ascending)')

    args = parser.parse_args()

    # Read and parse JSON file
    try:
        with open(args.json_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File '{args.json_file}' not found")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in file '{args.json_file}'")
        sys.exit(1)

    # Convert to sympy expression
    try:
        ascending_order = args.order == 'ascending'
        expr = parse_json_to_sympy(data, ascending=ascending_order)
        formatted_expr = format_expression(expr, args.format)

        print(f"Series from {Path(args.json_file).name}:")
        print("=" * 50)
        print(formatted_expr)

        if args.metadata and 'metadata' in data:
            print("\nMetadata:")
            print("-" * 20)
            metadata = data['metadata']
            for key, value in metadata.items():
                print(f"{key}: {value}")

    except Exception as e:
        print(f"Error processing expression: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()