import json
from pathlib import Path

import pytest
import sympy as sp

from bar_natan import T, bar_natan_Z
from bar_natan_links import bar_natan_link_Z, bar_natan_link_Z_from_braid
from test_bar_natan import braid_to_pd, normalize_pd


def _fk_coefficients(data: dict, *, derivative: bool = False) -> dict[tuple[sp.Rational, ...], sp.Expr]:
    metadata = data["metadata"]
    n_variables = int(metadata["num_x_variables"])
    overall_x = [
        sp.Rational(str(value))
        for value in metadata.get("overall_x_powers", [0] * n_variables)
    ]
    overall_q = sp.Rational(str(metadata.get("overall_q_power", 0)))

    coefficients = {}
    for term in data["terms"]:
        key = tuple(sp.Rational(power) + shift for power, shift in zip(term["x"], overall_x))
        coefficient = sp.Integer(0)
        for q_term in term["q_terms"]:
            q_power = sp.Rational(q_term["q"]) + overall_q
            q_coefficient = sp.Integer(q_term["c"])
            coefficient += q_power * q_coefficient if derivative else q_coefficient
        if coefficient:
            coefficients[key] = coefficients.get(key, sp.Integer(0)) + coefficient
    return coefficients


def _series_coefficients(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    orders: list[int],
) -> dict[tuple[sp.Rational, ...], sp.Expr]:
    series = expression
    for variable, order in reversed(list(zip(variables, orders))):
        series = sp.series(series, variable, 0, order + 3).removeO()
    series = sp.expand(series)

    coefficients = {}
    for term in sp.Add.make_args(series):
        powers = term.as_powers_dict()
        exponent = []
        coefficient = term
        for variable in variables:
            power = sp.Rational(powers.get(variable, 0))
            exponent.append(power)
            coefficient /= variable**power
        key = tuple(exponent)
        coefficients[key] = coefficients.get(key, sp.Integer(0)) + sp.simplify(coefficient)
    return {
        key: sp.simplify(value)
        for key, value in coefficients.items()
        if sp.simplify(value) != 0
    }


def _series_orders(data: dict) -> list[int]:
    metadata = data["metadata"]
    maxima = list(metadata["max_x_degrees"])
    for term in data["terms"]:
        for index, power in enumerate(term["x"]):
            maxima[index] = max(maxima[index], abs(int(power)))
    return maxima


def _expanded_series_orders(data: dict, margin: int = 6) -> list[int]:
    return [degree + margin for degree in _series_orders(data)]


def _matches_up_to_monomial(
    expected: dict[tuple[sp.Rational, ...], sp.Expr],
    actual: dict[tuple[sp.Rational, ...], sp.Expr],
) -> tuple[bool, tuple[sp.Rational, ...] | None, sp.Expr | None]:
    if not actual:
        return False, None, None
    if not expected:
        return False, None, None

    anchor_key, anchor_coefficient = next(iter(actual.items()))
    for expected_key, expected_coefficient in expected.items():
        if expected_coefficient == 0:
            continue
        shift = tuple(a - e for a, e in zip(anchor_key, expected_key))
        scale = sp.simplify(anchor_coefficient / expected_coefficient)
        if all(
            sp.simplify(scale * expected.get(tuple(k - s for k, s in zip(key, shift)), 0) - coefficient) == 0
            for key, coefficient in actual.items()
        ):
            return True, shift, scale
    return False, None, None


def test_link_z_specializes_to_knot_z_for_one_component_pd():
    pd_crossings = normalize_pd(braid_to_pd([1, 1, 1]))

    assert bar_natan_link_Z(pd_crossings, Tnames=[T]) == bar_natan_Z(pd_crossings, Tname=T)
    assert bar_natan_link_Z(pd_crossings) == bar_natan_Z(pd_crossings, Tname=T)


def test_link_z_from_braid_known_alexander_examples():
    x, y = sp.symbols("x y")

    hopf_p0, hopf_p1 = bar_natan_link_Z_from_braid([1, 1], Tnames=[x, y])
    assert sp.simplify(hopf_p0 - 1) == 0
    assert hopf_p1 != 0

    torus_p0, torus_p1 = bar_natan_link_Z_from_braid([1, 1, 1, 1], Tnames=[x, y])
    assert sp.simplify(torus_p0 * sp.sqrt(x) * sp.sqrt(y) - (1 + x * y)) == 0
    assert torus_p1 != 0

    mixed_p0, mixed_p1 = bar_natan_link_Z_from_braid([2, 2, 1, -2, 1], Tnames=[x, y])
    assert sp.simplify(mixed_p0 * x * sp.sqrt(y) - (x + y)) == 0
    assert mixed_p1 != 0


def test_bar_natan_p0_p1(json_file: Path):
    with json_file.open() as f:
        data = json.load(f)

    metadata = data["metadata"]
    components = int(metadata["components"])
    variables = tuple(sp.symbols(f"x0:{components}"))

    p0, p1 = bar_natan_link_Z_from_braid(metadata["braid"], Tnames=variables)
    assert p0 != 0, f"{json_file.name}: P0 vanished"

    orders = _expanded_series_orders(data)

    # P0 check: 1/P0 matches FK(q=1) coefficients up to a monomial normalisation
    expected_p0 = _series_coefficients(1 / p0, variables, orders)
    actual_p0 = _fk_coefficients(data)
    ok0, shift, scale = _matches_up_to_monomial(expected_p0, actual_p0)
    assert ok0, f"{json_file.name}: P0 mismatch; shift={shift}, scale={scale}"

    # P1 check: -P1/P0³ matches (d/dq FK)|_{q=1} with the same normalisation
    actual_p1 = _fk_coefficients(data, derivative=True)
    if not actual_p1:
        return
    expected_p1 = _series_coefficients(-p1 / p0**3, variables, orders)
    ok1, shift1, scale1 = _matches_up_to_monomial(expected_p1, actual_p1)
    assert ok1, f"{json_file.name}: P1 mismatch; shift={shift1}, scale={scale1}"


def test_bar_natan_links_folder(json_file: Path):
    with json_file.open() as f:
        data = json.load(f)

    metadata = data["metadata"]
    components = int(metadata["components"])
    variables = tuple(sp.symbols(f"x0:{components}"))

    p0, p1 = bar_natan_link_Z_from_braid(metadata["braid"], Tnames=variables)
    assert p0 != 0, f"{json_file.name}: P0 vanished"
    assert p1 is not None, f"{json_file.name}: P1 was not computed"

    expected = _series_coefficients(1 / p0, variables, _expanded_series_orders(data))
    actual = _fk_coefficients(data)
    ok, shift, scale = _matches_up_to_monomial(expected, actual)
    assert ok, f"{json_file.name}: q=1 term did not match 1/P0 up to monomial/sign; shift={shift}, scale={scale}"


if __name__ == "__main__":
    _LINKS_DIR = Path(__file__).parent / "links"
    _EXAMPLES = [
        "L2a1{0}_s0.json",
        "L4a1{0}_s0.json",
        "L6a2{0}_s0.json",
        "L6a4{0;1}_s0.json",
        "L8n3{0;1}_s0.json",
    ]

    passed = 0
    failed = 0
    for name in _EXAMPLES:
        json_file = _LINKS_DIR / name
        with json_file.open() as f:
            data = json.load(f)
        metadata = data["metadata"]
        components = int(metadata["components"])
        variables = tuple(sp.symbols(f"x0:{components}"))
        orders = _expanded_series_orders(data)

        p0, p1 = bar_natan_link_Z_from_braid(metadata["braid"], Tnames=variables)

        exp_p0 = _series_coefficients(1 / p0, variables, orders)
        act_p0 = _fk_coefficients(data)
        ok0, shift, scale = _matches_up_to_monomial(exp_p0, act_p0)

        act_p1 = _fk_coefficients(data, derivative=True)
        if act_p1:
            exp_p1 = _series_coefficients(-p1 / p0**3, variables, orders)
            ok1, shift1, scale1 = _matches_up_to_monomial(exp_p1, act_p1)
        else:
            ok1 = True

        status = "PASS" if (ok0 and ok1) else (
            f"FAIL (P0={'ok' if ok0 else 'FAIL'}, P1={'ok' if ok1 else 'FAIL'})"
        )
        print(f"{name}: {status}")
        if ok0 and ok1:
            passed += 1
        else:
            failed += 1

    print(f"\n{passed} passed, {failed} failed")
