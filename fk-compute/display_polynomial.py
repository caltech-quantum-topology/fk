#!/usr/bin/env python3
"""
Script to load and display multivariable polynomials from the results folder
in a human-readable format with highlighted q-polynomials for each x power.
"""

import json
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple

def load_polynomial(filepath: str) -> Dict:
    """Load polynomial data from JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)

def format_q_polynomial(q_terms: List[Dict]) -> str:
    """Format q-polynomial terms into a readable string."""
    if not q_terms:
        return "0"

    # Sort by q power in ascending order
    sorted_terms = sorted(q_terms, key=lambda x: x['q'])

    terms = []
    for i, term in enumerate(sorted_terms):
        q_power = term['q']
        coeff = term['c']

        if coeff == 0:
            continue

        # Format coefficient
        if i == 0:  # First term
            if coeff == 1:
                coeff_str = ""
            elif coeff == -1:
                coeff_str = "-"
            else:
                coeff_str = str(coeff)
        else:  # Subsequent terms
            if coeff == 1:
                coeff_str = " + "
            elif coeff == -1:
                coeff_str = " - "
            elif coeff > 0:
                coeff_str = f" + {coeff}"
            else:
                coeff_str = f" - {abs(coeff)}"

        # Format q power
        if q_power == 0:
            q_str = ""
        elif q_power == 1:
            q_str = "q"
        else:
            q_str = f"q^{q_power}"

        # Combine coefficient and q power
        if q_power == 0:
            term_str = coeff_str if coeff_str not in ["", " + ", " - "] else (coeff_str + "1")
        else:
            term_str = coeff_str + q_str

        terms.append(term_str)

    result = "".join(terms)
    return result if result else "0"

def format_x_power(x_powers: List[int]) -> str:
    """Format x powers into a readable string."""
    if not x_powers:
        return ""

    terms = []
    for i, power in enumerate(x_powers):
        if power == 0:
            continue
        elif power == 1:
            terms.append(f"x_{i}")
        else:
            terms.append(f"x_{i}^{power}")

    return "".join(terms) if terms else "1"

def display_polynomial(data: Dict) -> None:
    """Display polynomial in human-readable format."""
    terms = data.get('terms', [])
    metadata = data.get('metadata', {})

    print("Polynomial:")
    print("=" * 60)

    if not terms:
        print("0")
        return

    # Sort terms by x powers in ascending order
    # For single variable, sort by x[0], for multivariable, sort lexicographically
    def sort_key(term):
        x_powers = term.get('x', [])
        # Pad with zeros to ensure consistent comparison
        max_vars = metadata.get('num_x_variables', 1)
        padded = x_powers + [0] * (max_vars - len(x_powers))
        return padded

    sorted_terms = sorted(terms, key=sort_key)

    for i, term in enumerate(sorted_terms):
        x_powers = term.get('x', [])
        q_terms = term.get('q_terms', [])

        # Format x part
        x_str = format_x_power(x_powers)
        if not x_str:
            x_str = "1"

        # Format q polynomial
        q_str = format_q_polynomial(q_terms)

        # Combine with highlighting
        if i == 0:
            prefix = ""
        else:
            prefix = "\n+ "

        print(f"{prefix}({q_str}) * {x_str}")

    # Display metadata
    print("\n" + "=" * 60)
    print("Metadata:")
    print(f"Number of x variables: {metadata.get('num_x_variables', 'Unknown')}")
    print(f"Maximum x degrees: {metadata.get('max_x_degrees', 'Unknown')}")
    print(f"Storage type: {metadata.get('storage_type', 'Unknown')}")
    print(f"Components: {metadata.get('components', 'Unknown')}")

def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) != 2:
        print("Usage: python display_polynomial.py <filename>")
        print("\nAvailable files in results/:")
        results_dir = Path("results")
        if results_dir.exists():
            for file in sorted(results_dir.glob("*.json")):
                print(f"  {file.name}")
        sys.exit(1)

    filename = sys.argv[1]

    # Check if it's just a filename (look in results/) or a full path
    if not os.path.sep in filename and not filename.startswith('.'):
        filepath = Path("results") / filename
    else:
        filepath = Path(filename)

    if not filepath.exists():
        print(f"Error: File '{filepath}' not found.")
        sys.exit(1)

    try:
        data = load_polynomial(filepath)
        print(f"Displaying polynomial from: {filepath}")
        print()
        display_polynomial(data)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in file '{filepath}': {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()