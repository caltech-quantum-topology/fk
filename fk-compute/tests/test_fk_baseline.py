"""
Test FK computation against baseline data.

This test computes FK homology for all knots up to 11 crossings up to degree 10
and compares the results against pre-computed baseline data.
"""

import json
import re
import multiprocessing as mp
from pathlib import Path

import pytest
import sympy
import yaml
from sympy.core.expr import Expr

import fkcompute


# Constants
MAX_CROSSINGS = 11
MAX_DEGREE = 10
TESTS_DIR = Path(__file__).parent
CONFIGS_DIR = TESTS_DIR / "configs" / "config_10"
BASELINE_DIR = TESTS_DIR / "data_baseline"


def data_file_to_sympy(fk_data: dict) -> Expr:
    """Convert FK computation result (dict) to sympy expression."""
    q = sympy.Symbol("q")

    meta = fk_data.get("metadata", {})
    num_x_vars = int(meta.get("num_x_variables", 1) or 1)

    if num_x_vars == 1:
        x_vars = (sympy.Symbol("x"),)
    else:
        x_vars = tuple(sympy.Symbol(f"x{i}") for i in range(num_x_vars))

    def _coeff_from_str(s: str) -> Expr:
        s = str(s)
        return sympy.Rational(s) if "/" in s else sympy.Integer(s)

    expr = sympy.Integer(0)
    for term in fk_data.get("terms", []):
        degrees = list(term.get("x", []))
        if len(degrees) < num_x_vars:
            degrees.extend([0] * (num_x_vars - len(degrees)))

        x_monom = sympy.Integer(1)
        for i, deg in enumerate(degrees[:num_x_vars]):
            deg_i = int(deg)
            if deg_i:
                x_monom *= x_vars[i] ** deg_i

        q_poly = sympy.Integer(0)
        for qt in term.get("q_terms", []):
            q_exp = int(qt.get("q"))
            c = _coeff_from_str(qt.get("c"))
            q_poly += c * (q ** q_exp)

        expr += x_monom * q_poly

    return sympy.expand(expr)


def truncate_to_x_degree(expr: Expr, max_x_degree: int) -> Expr:
    """Truncate expression to maximum x degree."""
    x = sympy.Symbol("x")
    if expr == 0:
        return sympy.Integer(0)

    poly = sympy.Poly(expr, x, domain="EX")

    out = sympy.Integer(0)
    for monom, coeff in poly.terms():
        x_deg = int(monom[0])
        if x_deg <= max_x_degree:
            out += coeff * (x ** x_deg)

    return sympy.expand(out)


def max_x_degree(expr: Expr) -> int:
    """Get maximum x degree of expression."""
    x = sympy.Symbol("x")
    if expr == 0:
        return -1
    return int(sympy.Poly(expr, x, domain="EX").degree())


def get_crossing_number(knot_name: str) -> int | None:
    """Extract crossing number from knot name (e.g., '10_100' -> 10, '11a_3' -> 11)."""
    try:
        prefix = knot_name.split('_')[0]
        match = re.match(r'^(\d+)', prefix)
        if match:
            return int(match.group(1))
        return None
    except (ValueError, IndexError):
        return None


def get_knot_configs() -> list[tuple[str, dict]]:
    """Load all knot configurations up to MAX_CROSSINGS."""
    configs = []

    for config_file in CONFIGS_DIR.glob("config_*.yaml"):
        knot_name = config_file.stem.replace("config_", "")
        crossing_num = get_crossing_number(knot_name)

        if crossing_num is not None and crossing_num <= MAX_CROSSINGS:
            baseline_file = BASELINE_DIR / f"{knot_name}.json"
            if baseline_file.exists():
                with open(config_file) as f:
                    config = yaml.safe_load(f)
                configs.append((knot_name, config))

    return sorted(configs, key=lambda x: (get_crossing_number(x[0]) or 0, x[0]))


def compute_fk(knot_name: str, config: dict) -> dict:
    """Compute FK polynomial for a knot."""
    braid = config.get("braid")
    inversion_data = config.get("inversion")
    inversion = {"inversion_data": inversion_data, "braid": braid}

    result = fkcompute.fk(
        braid,
        degree=MAX_DEGREE,
        name=knot_name,
        threads=mp.cpu_count(),
        save_data=False,
        inversion=inversion,
    )

    return result


def load_baseline(knot_name: str) -> dict:
    """Load baseline FK data for a knot."""
    baseline_file = BASELINE_DIR / f"{knot_name}.json"
    with open(baseline_file) as f:
        return json.load(f)


# Collect test parameters
KNOT_CONFIGS = get_knot_configs()


@pytest.mark.parametrize("knot_name,config", KNOT_CONFIGS, ids=[k[0] for k in KNOT_CONFIGS])
def test_fk_matches_baseline(knot_name: str, config: dict):
    """Test that computed FK matches baseline for a single knot."""
    # Compute FK
    computed_result = compute_fk(knot_name, config)
    computed_expr = data_file_to_sympy(computed_result)

    # Load baseline
    baseline_result = load_baseline(knot_name)
    baseline_expr = data_file_to_sympy(baseline_result)

    # Find common degree and truncate
    computed_degree = max_x_degree(computed_expr)
    baseline_degree = max_x_degree(baseline_expr)
    common_degree = min(computed_degree, baseline_degree)

    if common_degree < 0:
        # Both are zero
        assert computed_expr == 0 and baseline_expr == 0, \
            f"Knot {knot_name}: one expression is zero, the other is not"
        return

    truncated_computed = truncate_to_x_degree(computed_expr, common_degree)
    truncated_baseline = truncate_to_x_degree(baseline_expr, common_degree)

    # Compare
    diff = sympy.expand(truncated_computed - truncated_baseline)
    assert diff == 0, \
        f"Knot {knot_name}: FK mismatch at common degree {common_degree}\n" \
        f"Computed: {truncated_computed}\n" \
        f"Baseline: {truncated_baseline}\n" \
        f"Difference: {diff}"


if __name__ == "__main__":
    # Run tests directly for debugging
    print(f"Testing FK computation against baseline")
    print(f"Max crossings: {MAX_CROSSINGS}")
    print(f"Max degree: {MAX_DEGREE}")
    print(f"Number of knots: {len(KNOT_CONFIGS)}")
    print("=" * 60)

    passed = 0
    failed = 0

    for knot_name, config in KNOT_CONFIGS:
        try:
            test_fk_matches_baseline(knot_name, config)
            print(f"[PASS] {knot_name}")
            passed += 1
        except AssertionError as e:
            print(f"[FAIL] {knot_name}: {e}")
            failed += 1
        except Exception as e:
            print(f"[ERROR] {knot_name}: {e}")
            failed += 1

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
