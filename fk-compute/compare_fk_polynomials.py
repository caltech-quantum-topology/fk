#!/usr/bin/env python3
"""
FK Polynomial Comparison Tool

Compares FK polynomials stored in two different data sources:
1. fibered_table.json - source of truth with fk_x_coeffs format
2. results/*.json - computed results in explicit term format

Polynomials are compared up to the minimum x-degree and can differ by a factor of -1.
"""

import json
import re
import os
from collections import defaultdict
from typing import Dict, Tuple, List, Optional


def parse_laurent_polynomial(expr: str) -> Dict[int, float]:
    """
    Parse a Laurent polynomial in q from string format into dictionary.

    Args:
        expr: String like "7 + 28*q - q^2" or "q^(-2) + 22/q"

    Returns:
        Dictionary mapping q-power to coefficient: {q_power: coefficient}
    """
    import re

    # Handle zero polynomial
    if expr.strip() == "0":
        return {}

    # Dictionary to store coefficients for each power of q
    coeffs = defaultdict(float)

    # Remove spaces and normalize the expression
    expr = expr.replace(' ', '')

    # Define regex patterns for different term types
    patterns = [
        # Coefficient * q^power patterns (including negative powers)
        (r'([+-]?\d*\.?\d*)\*q\^\((-?\d+)\)', lambda m: (float(m.group(1) or '1'), int(m.group(2)))),
        (r'([+-]?\d*\.?\d*)\*q\^(-?\d+)', lambda m: (float(m.group(1) or '1'), int(m.group(2)))),

        # q^power patterns
        (r'([+-]?)q\^\((-?\d+)\)', lambda m: (1.0 if m.group(1) != '-' else -1.0, int(m.group(2)))),
        (r'([+-]?)q\^(-?\d+)', lambda m: (1.0 if m.group(1) != '-' else -1.0, int(m.group(2)))),

        # Coefficient * q patterns
        (r'([+-]?\d*\.?\d*)\*q(?![0-9^])', lambda m: (float(m.group(1) or '1'), 1)),

        # Just q patterns
        (r'([+-]?)q(?![0-9^])', lambda m: (1.0 if m.group(1) != '-' else -1.0, 1)),

        # Fraction patterns: coeff/q^power
        (r'([+-]?\d*\.?\d*)/q\^(-?\d+)', lambda m: (float(m.group(1) or '1'), -int(m.group(2)))),

        # Fraction patterns: coeff/q
        (r'([+-]?\d*\.?\d*)/q(?![0-9^])', lambda m: (float(m.group(1) or '1'), -1)),

        # Constant terms (no q)
        (r'([+-]?\d+\.?\d*)(?![*/])', lambda m: (float(m.group(1)), 0)),
    ]

    # Make a copy of the expression to work with
    remaining = expr

    # Remove leading '+' if present
    if remaining.startswith('+'):
        remaining = remaining[1:]

    # Split into terms while preserving signs
    terms = re.findall(r'[+-]?[^+-]+', remaining)

    for term in terms:
        term = term.strip()
        if not term:
            continue

        # Try each pattern
        matched = False
        for pattern, extractor in patterns:
            match = re.fullmatch(pattern, term)
            if match:
                try:
                    coeff, power = extractor(match)
                    # Handle empty coefficient strings
                    if isinstance(coeff, str) and coeff in ['', '+', '-']:
                        coeff = 1.0 if coeff != '-' else -1.0
                    coeffs[power] += coeff
                    matched = True
                    break
                except (ValueError, IndexError):
                    continue

        if not matched:
            # Try to handle special cases manually
            if 'q^(' in term and ')' in term:
                # Handle q^(-n) format
                try:
                    if term.startswith('-'):
                        sign = -1
                        term = term[1:]
                    elif term.startswith('+'):
                        sign = 1
                        term = term[1:]
                    else:
                        sign = 1

                    if term.startswith('q^(') and term.endswith(')'):
                        power_str = term[3:-1]
                        power = int(power_str)
                        coeffs[power] += sign * 1
                        matched = True
                except (ValueError, IndexError):
                    pass

    # Remove zero coefficients and return
    return {k: v for k, v in coeffs.items() if abs(v) > 1e-10}


def parse_fk_from_fibered_table(entry: Dict) -> Dict[Tuple[int, int], float]:
    """
    Parse FK polynomial from fibered_table.json entry.

    Args:
        entry: Dictionary containing 'fk_x_coeffs' field

    Returns:
        Dictionary mapping (x_degree, q_power) to coefficient
    """
    fk_coeffs_str = entry['fk_x_coeffs']

    # Remove curly braces and split by commas
    coeffs_str = fk_coeffs_str.strip('{}')
    coeffs_list = [c.strip() for c in coeffs_str.split(',')]

    polynomial = {}

    for x_degree, coeff_expr in enumerate(coeffs_list):
        if coeff_expr.strip() == '0':
            continue

        # Parse the Laurent polynomial for this x-degree
        q_coeffs = parse_laurent_polynomial(coeff_expr)

        for q_power, coeff in q_coeffs.items():
            polynomial[(x_degree, q_power)] = coeff

    return polynomial


def parse_fk_from_result(result_data: Dict) -> Dict[Tuple[int, int], float]:
    """
    Parse FK polynomial from results/*.json format.

    Args:
        result_data: Dictionary containing 'terms' field

    Returns:
        Dictionary mapping (x_degree, q_power) to coefficient
    """
    polynomial = {}

    for term in result_data['terms']:
        x_degree = term['x'][0]  # x-degree

        for q_term in term['q_terms']:
            q_power = q_term['q']
            coeff = q_term['c']
            polynomial[(x_degree, q_power)] = coeff

    return polynomial


def get_min_x_degree(poly1: Dict[Tuple[int, int], float],
                    poly2: Dict[Tuple[int, int], float]) -> int:
    """Get the minimum x-degree present in both polynomials."""
    x_degrees1 = {x for x, q in poly1.keys()}
    x_degrees2 = {x for x, q in poly2.keys()}

    if not x_degrees1 or not x_degrees2:
        return 0

    return min(max(x_degrees1), max(x_degrees2))


def truncate_to_min_x_degree(polynomial: Dict[Tuple[int, int], float],
                           min_x_degree: int) -> Dict[Tuple[int, int], float]:
    """Truncate polynomial to only include terms up to min_x_degree."""
    return {(x, q): c for (x, q), c in polynomial.items() if x <= min_x_degree}


def compare_polynomials(poly1: Dict[Tuple[int, int], float],
                       poly2: Dict[Tuple[int, int], float]) -> Tuple[bool, str, Optional[int]]:
    """
    Compare two FK polynomials up to minimum x-degree.

    Returns:
        (matches, status_message, sign_factor)
        - matches: True if polynomials match (possibly with sign flip)
        - status_message: Description of the comparison result
        - sign_factor: 1 or -1 if they match with sign flip, None otherwise
    """
    # Get minimum x-degree
    min_x_deg = get_min_x_degree(poly1, poly2)

    # Truncate both polynomials
    trunc_poly1 = truncate_to_min_x_degree(poly1, min_x_deg)
    trunc_poly2 = truncate_to_min_x_degree(poly2, min_x_deg)

    # Get all keys from both polynomials
    all_keys = set(trunc_poly1.keys()) | set(trunc_poly2.keys())

    if not all_keys:
        return True, "Both polynomials are zero", 1

    # Check for exact match
    exact_match = True
    for key in all_keys:
        coeff1 = trunc_poly1.get(key, 0)
        coeff2 = trunc_poly2.get(key, 0)
        if abs(coeff1 - coeff2) > 1e-10:
            exact_match = False
            break

    if exact_match:
        return True, "Exact match", 1

    # Check for match with sign flip
    sign_match = True
    for key in all_keys:
        coeff1 = trunc_poly1.get(key, 0)
        coeff2 = trunc_poly2.get(key, 0)
        if abs(coeff1 + coeff2) > 1e-10:
            sign_match = False
            break

    if sign_match:
        return True, "Match with sign flip (-1 factor)", -1

    # Calculate differences for reporting
    max_diff = 0
    for key in all_keys:
        coeff1 = trunc_poly1.get(key, 0)
        coeff2 = trunc_poly2.get(key, 0)
        max_diff = max(max_diff, abs(coeff1 - coeff2))

    return False, f"No match (max difference: {max_diff:.6f})", None


def main():
    """Main comparison function."""
    print("FK Polynomial Comparison Tool")
    print("=" * 50)

    # Load fibered_table.json
    print("Loading fibered_table.json...")
    try:
        with open('fibered_table.json', 'r') as f:
            fibered_data = json.load(f)
    except FileNotFoundError:
        print("Error: fibered_table.json not found!")
        return

    # Get list of result files
    results_dir = 'results'
    if not os.path.exists(results_dir):
        print(f"Error: {results_dir} directory not found!")
        return

    result_files = [f for f in os.listdir(results_dir) if f.endswith('.json')]
    print(f"Found {len(result_files)} result files")

    # Create lookup for fibered data by name
    fibered_lookup = {entry['name']: entry for entry in fibered_data}

    # Track results
    total_comparisons = 0
    exact_matches = 0
    sign_matches = 0
    mismatches = 0
    missing_in_fibered = 0
    missing_in_results = 0

    print("\nComparison Results:")
    print("-" * 70)

    for result_file in sorted(result_files):
        knot_name = result_file.replace('.json', '')

        # Load result data
        try:
            with open(os.path.join(results_dir, result_file), 'r') as f:
                result_data = json.load(f)
        except Exception as e:
            print(f"❌ {knot_name:<12} Error loading result file: {e}")
            continue

        # Check if knot exists in fibered data
        if knot_name not in fibered_lookup:
            print(f"❓ {knot_name:<12} Not found in fibered_table.json")
            missing_in_fibered += 1
            continue

        total_comparisons += 1

        try:
            # Parse polynomials
            fibered_poly = parse_fk_from_fibered_table(fibered_lookup[knot_name])
            result_poly = parse_fk_from_result(result_data)

            # Compare polynomials
            matches, status, sign_factor = compare_polynomials(fibered_poly, result_poly)

            if matches:
                if sign_factor == 1:
                    print(f"✅ {knot_name:<12} {status}")
                    exact_matches += 1
                else:
                    print(f"🔄 {knot_name:<12} {status}")
                    sign_matches += 1
            else:
                print(f"❌ {knot_name:<12} {status}")
                mismatches += 1

                # Print details for specific knots to debug
                if knot_name in ['10_105', '10_125', '10_141']:
                    min_x_deg = get_min_x_degree(fibered_poly, result_poly)
                    trunc_fibered = truncate_to_min_x_degree(fibered_poly, min_x_deg)
                    trunc_result = truncate_to_min_x_degree(result_poly, min_x_deg)

                    print(f"   Debugging {knot_name} (up to x^{min_x_deg}):")

                    # Show fibered polynomial structure
                    print(f"   Fibered: {fibered_lookup[knot_name]['fk_x_coeffs']}")

                    all_keys = set(trunc_fibered.keys()) | set(trunc_result.keys())
                    print(f"   Comparison of terms:")
                    count = 0
                    for key in sorted(all_keys):
                        if count >= 6:  # Limit output
                            break
                        coeff_f = trunc_fibered.get(key, 0)
                        coeff_r = trunc_result.get(key, 0)
                        if abs(coeff_f - coeff_r) > 1e-10:
                            print(f"     x^{key[0]}*q^{key[1]}: fibered={coeff_f}, result={coeff_r}")
                            count += 1

        except Exception as e:
            print(f"❌ {knot_name:<12} Error during comparison: {e}")
            mismatches += 1

    # Check for knots in fibered_table but not in results
    result_names = {f.replace('.json', '') for f in result_files}
    for entry in fibered_data:
        if entry['name'] not in result_names:
            missing_in_results += 1

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total comparisons performed: {total_comparisons}")
    print(f"✅ Exact matches: {exact_matches}")
    print(f"🔄 Sign-flipped matches: {sign_matches}")
    print(f"❌ Mismatches: {mismatches}")
    print(f"❓ Missing in fibered_table: {missing_in_fibered}")
    print(f"❓ Missing in results: {missing_in_results}")

    if total_comparisons > 0:
        success_rate = (exact_matches + sign_matches) / total_comparisons * 100
        print(f"Overall success rate: {success_rate:.1f}%")


if __name__ == "__main__":
    main()