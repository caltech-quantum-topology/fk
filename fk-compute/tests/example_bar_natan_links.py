"""
Example: compute Bar-Natan P1 for selected link braids and compare to FK JSON.

Run from this directory:

    python example_bar_natan_links.py

Useful options:

    python example_bar_natan_links.py --case L8a10 --print-p1
    python example_bar_natan_links.py --case L10n10 --factor
    python example_bar_natan_links.py --max-degree 6
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

from bar_natan_links import bar_natan_link_Z_from_braid


@dataclass(frozen=True)
class LinkExample:
    key: str
    display_name: str
    braid: list[int]
    preferred_json_names: tuple[str, ...]


EXAMPLES = [
    LinkExample(
        key="L8a10",
        display_name="L8a10{1}_s0",
        braid=[-1, 2, 3, 2, 2, 2, -1, 2, -3, 2],
        preferred_json_names=("L8a10{0}_s0.json", "L8a10{1}_s0.json"),
    ),
    LinkExample(
        key="L10n10",
        display_name="L10n10{0}_s0",
        braid=[1, 2, -3, 2, -1, -3, -2, -3, -2, -3, -2, -3],
        preferred_json_names=("L10n10{0}_s0.json",),
    ),
]


def load_fk_json(example: LinkExample, links_dir: Path) -> tuple[Path, dict]:
    for name in example.preferred_json_names:
        candidate = links_dir / name
        if candidate.exists():
            with candidate.open() as f:
                return candidate, json.load(f)

    for candidate in sorted(links_dir.glob("*.json")):
        if candidate.name.endswith("_inversion.json"):
            continue
        with candidate.open() as f:
            data = json.load(f)
        if data.get("metadata", {}).get("braid") == example.braid:
            return candidate, data

    raise FileNotFoundError(
        f"No FK JSON in {links_dir} matched {example.display_name} or braid {example.braid}"
    )


def fk_coefficients(data: dict, *, derivative: bool = False) -> dict[tuple[sp.Rational, ...], sp.Expr]:
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


def series_coefficients(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    max_degree: int,
) -> dict[tuple[sp.Rational, ...], sp.Expr]:
    series = expression
    for variable in reversed(variables):
        series = sp.series(series, variable, 0, max_degree + 1).removeO()
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


def truncate_coefficients(
    coefficients: dict[tuple[sp.Rational, ...], sp.Expr],
    max_degree: int,
) -> dict[tuple[sp.Rational, ...], sp.Expr]:
    return {
        key: value
        for key, value in coefficients.items()
        if all(-max_degree <= exponent <= max_degree for exponent in key)
    }


def matches_up_to_monomial(
    expected: dict[tuple[sp.Rational, ...], sp.Expr],
    actual: dict[tuple[sp.Rational, ...], sp.Expr],
) -> tuple[bool, tuple[sp.Rational, ...] | None, sp.Expr | None]:
    if not expected or not actual:
        return False, None, None

    anchor_key, anchor_coefficient = next(iter(actual.items()))
    for expected_key, expected_coefficient in expected.items():
        if expected_coefficient == 0:
            continue
        shift = tuple(a - e for a, e in zip(anchor_key, expected_key))
        scale = sp.simplify(anchor_coefficient / expected_coefficient)
        if all(
            sp.simplify(
                scale * expected.get(tuple(k - s for k, s in zip(key, shift)), 0)
                - coefficient
            )
            == 0
            for key, coefficient in actual.items()
        ):
            return True, shift, scale
    return False, None, None


def format_coefficients(coefficients: dict[tuple[sp.Rational, ...], sp.Expr]) -> str:
    if not coefficients:
        return "0"
    pieces = []
    for exponent in sorted(coefficients):
        pieces.append(f"{exponent}: {sp.simplify(coefficients[exponent])}")
    return "\n".join(pieces)


def shifted_scaled_coefficients(
    expected: dict[tuple[sp.Rational, ...], sp.Expr],
    actual_keys: set[tuple[sp.Rational, ...]],
    shift: tuple[sp.Rational, ...],
    scale: sp.Expr,
) -> dict[tuple[sp.Rational, ...], sp.Expr]:
    shifted = {}
    for key in actual_keys:
        expected_key = tuple(k - s for k, s in zip(key, shift))
        shifted[key] = sp.simplify(scale * expected.get(expected_key, 0))
    return shifted


def expression_summary(name: str, expression: sp.Expr, *, print_full: bool) -> None:
    if print_full:
        print(f"{name} = {sp.factor(expression)}")
        return

    print(f"{name}:")
    print(f"  operations: {sp.count_ops(expression)}")
    print(f"  string length: {len(str(expression))}")


def compare_series(
    label: str,
    expected_expression: sp.Expr,
    actual_coefficients: dict[tuple[sp.Rational, ...], sp.Expr],
    variables: tuple[sp.Symbol, ...],
    max_degree: int,
    expansion_margin: int,
    print_series: bool,
) -> None:
    expected_coefficients = series_coefficients(
        expected_expression,
        variables,
        max_degree + expansion_margin,
    )
    actual_coefficients = truncate_coefficients(actual_coefficients, max_degree)
    ok, shift, scale = matches_up_to_monomial(expected_coefficients, actual_coefficients)
    status = "MATCH" if ok else "NO MATCH"
    print(f"  {label}: {status}")
    if ok:
        print(f"    shift={shift}, scale={scale}")
    if print_series:
        print("    actual FK coefficients:")
        print(format_coefficients(actual_coefficients))
        if ok:
            aligned_expected = shifted_scaled_coefficients(
                expected_coefficients,
                set(actual_coefficients),
                shift,
                scale,
            )
            print("    expected coefficients after shift/scale:")
        else:
            aligned_expected = {
                key: expected_coefficients[key]
                for key in sorted(expected_coefficients)
                if all(-max_degree <= exponent <= max_degree for exponent in key)
            }
            print("    expected coefficients without alignment:")
        print(format_coefficients(aligned_expected))


def run_example(example: LinkExample, args: argparse.Namespace) -> None:
    json_path, fk_data = load_fk_json(example, args.links_dir)
    components = int(fk_data["metadata"]["components"])
    variables = tuple(sp.symbols(f"x0:{components}"))

    print(f"\n{example.display_name}", flush=True)
    print(f"  FK JSON: {json_path}", flush=True)
    if json_path.name not in example.preferred_json_names:
        print(f"  resolved by exact braid match", flush=True)
    print(f"  braid: {example.braid}", flush=True)
    print("  computing P0 and P1...", flush=True)

    p0, p1 = bar_natan_link_Z_from_braid(
        example.braid,
        Tnames=variables,
        factor_result=args.factor,
    )

    expression_summary("  P0", p0, print_full=args.print_p0)
    expression_summary("  P1", p1, print_full=args.print_p1)

    q1_coefficients = fk_coefficients(fk_data)
    derivative_coefficients = fk_coefficients(fk_data, derivative=True)

    print(f"  FK comparison through degree {args.max_degree}:")
    compare_series(
        "q=1 vs 1/P0",
        1 / p0,
        q1_coefficients,
        variables,
        args.max_degree,
        args.expansion_margin,
        args.print_series,
    )
    if args.compare_p1:
        compare_series(
            "d/dq|q=1 vs -P1/P0^3",
            -p1 / p0**3,
            derivative_coefficients,
            variables,
            args.max_degree,
            args.expansion_margin,
            args.print_series,
        )
        compare_series(
            "d/dq|q=1 vs +P1/P0^3",
            p1 / p0**3,
            derivative_coefficients,
            variables,
            args.max_degree,
            args.expansion_margin,
            args.print_series,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=["all", *(example.key for example in EXAMPLES)],
        default="all",
        help="Which example to run.",
    )
    parser.add_argument(
        "--links-dir",
        type=Path,
        default=Path("links"),
        help="Directory containing FK JSON files.",
    )
    parser.add_argument(
        "--max-degree",
        type=int,
        default=6,
        help="Degree window used for FK series comparisons.",
    )
    parser.add_argument(
        "--expansion-margin",
        type=int,
        default=6,
        help="Extra formal-series degree used to avoid false negatives after monomial shifts.",
    )
    parser.add_argument(
        "--factor",
        action="store_true",
        help="Ask bar_natan_links to factor P0 and P1 before returning them.",
    )
    parser.add_argument(
        "--print-p0",
        action="store_true",
        help="Print full factored P0 expression.",
    )
    parser.add_argument(
        "--print-p1",
        action="store_true",
        help="Print full factored P1 expression. This can be large.",
    )
    parser.add_argument(
        "--compare-p1",
        action="store_true",
        help="Also compare d/dq at q=1 with the two P1/P0^3 sign conventions.",
    )
    parser.add_argument(
        "--print-series",
        action="store_true",
        help="Print the coefficient series being compared.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = EXAMPLES if args.case == "all" else [example for example in EXAMPLES if example.key == args.case]
    for example in selected:
        run_example(example, args)


if __name__ == "__main__":
    main()
