#!/usr/bin/env python3
import json
import argparse
import os
import sympy as sp


def load_json_poly(json_path):
    """Load polynomial from the JSON file and return a dict: x_deg -> sympy expr in q."""
    q, x = sp.symbols('q x')
    with open(json_path, 'r') as f:
        data = json.load(f)

    poly = {}
    for term in data["terms"]:
        x_deg = term["x"][0]
        coeff_expr = 0
        for qt in term["q_terms"]:
            q_deg = qt["q"]
            c = qt["c"]
            coeff_expr += c * q**q_deg
        poly[x_deg] = sp.simplify(coeff_expr)

    return poly, q, x


def parse_string_poly(poly_str, q, x):
    """Parse polynomial string into sympy expression."""
    poly_str = poly_str.replace('^', '**')
    expr = sp.sympify(poly_str, locals={'x': x, 'q': q})
    return sp.expand(expr)


def extract_global_sign_q_power(ratio, q):
    """Determine whether ratio = ± q**n for some integer n."""
    ratio = sp.simplify(ratio)

    if ratio == 1:
        return True, 1, 0
    if ratio == -1:
        return True, -1, 0

    try:
        r_at_1 = sp.simplify(ratio.subs(q, 1))
    except Exception:
        return False, None, None

    if not r_at_1.is_number or r_at_1 not in (1, -1):
        return False, None, None

    sign = int(r_at_1)
    t = sp.simplify(ratio / sign)

    if t == 1:
        return True, sign, 0
    if t == q:
        return True, sign, 1
    if isinstance(t, sp.Pow) and t.base == q and t.exp.is_integer:
        return True, sign, int(t.exp)

    return False, None, None


def compare_polynomials(json_poly, expr, q, x):
    """
    Check whether JSON polynomial P satisfies:
      P = ± q**n * fk_terms
    up to the minimum x-degree of the two polynomials.
    """
    max_x_json = max(json_poly.keys()) if json_poly else -sp.oo
    poly_in_x = sp.Poly(expr, x)
    max_x_expr = poly_in_x.degree()

    max_compare = min(max_x_json, max_x_expr)

    ratio = None

    for k in range(0, max_compare + 1):
        c_json = sp.simplify(json_poly.get(k, 0))
        c_expr = sp.simplify(expr.coeff(x, k))

        if c_json == 0 and c_expr == 0:
            continue

        if c_json == 0 or c_expr == 0:
            return False, None, None

        r = sp.simplify(c_json / c_expr)

        if ratio is None:
            ratio = r
        else:
            if sp.simplify(r - ratio) != 0:
                return False, None, None

    if ratio is None:
        # Everything zero → trivial match
        return True, 1, 0

    ok, sign, n = extract_global_sign_q_power(ratio, q)
    if not ok:
        return False, None, None

    return True, sign, n


def main():
    parser = argparse.ArgumentParser(
        description="Compare all JSON polynomials in data/ to fk_terms in fibered_table.json."
    )
    parser.add_argument("--data-dir", default="results",
                        help="Directory containing <knot>.json files (default: results)")
    parser.add_argument("--table-path", default="fibered_table.json",
                        help="Path to fibered_table.json")
    args = parser.parse_args()

    # Load the entire table
    with open(args.table_path, 'r') as f:
        table = json.load(f)

    results = []  # for summary output

    for entry in table:
        name = entry["name"]
        fk_terms = entry["fk_terms"]

        json_path = os.path.join(args.data_dir, f"{name}.json")
        if not os.path.exists(json_path):
            results.append((name, "MISSING_JSON", None, None))
            continue

        try:
            json_poly, q, x = load_json_poly(json_path)
            expr = parse_string_poly(fk_terms, q, x)
            match, sign, n = compare_polynomials(json_poly, expr, q, x)

            if match:
                results.append((name, "MATCH", sign, n))
            else:
                results.append((name, "NO_MATCH", None, None))
        except Exception as e:
            results.append((name, f"ERROR: {e}", None, None))

    # Print summary
    print("\n=== SUMMARY OF ALL KNOTS ===")
    for name, status, sign, n in results:
        if status == "MATCH":
            print(f"{name}: MATCH  (P_json = ({sign}) * q^{n} * P_table)")
        elif status == "NO_MATCH":
            print(f"{name}: NO MATCH")
        elif status == "MISSING_JSON":
            print(f"{name}: JSON FILE NOT FOUND")
        else:
            print(f"{name}: {status}")

    print("\nDone.")


if __name__ == "__main__":
    main()
'''
#!/usr/bin/env python3
import json
import argparse
import os
import sympy as sp


def load_json_poly(json_path):
    """Load polynomial from the JSON file and return a dict: x_deg -> sympy expr in q."""
    q, x = sp.symbols('q x')
    with open(json_path, 'r') as f:
        data = json.load(f)

    poly = {}
    for term in data["terms"]:
        x_deg = term["x"][0]
        coeff_expr = 0
        for qt in term["q_terms"]:
            q_deg = qt["q"]
            c = qt["c"]
            coeff_expr += c * q**q_deg
        poly[x_deg] = sp.simplify(coeff_expr)

    return poly, q, x


def load_fk_terms_from_table(table_path, knot_name):
    """Load fk_terms (string polynomial in x,q) for a given knot name."""
    with open(table_path, 'r') as f:
        table = json.load(f)

    for entry in table:
        if entry.get("name") == knot_name:
            return entry["fk_terms"]

    raise ValueError(f"Knot name {knot_name!r} not found in {table_path}")


def parse_string_poly(poly_str, q, x):
    """Parse polynomial string into sympy expression."""
    poly_str = poly_str.replace('^', '**')
    expr = sp.sympify(poly_str, locals={'x': x, 'q': q})
    return sp.expand(expr)


def extract_global_sign_q_power(ratio, q):
    """
    Given a sympy expression `ratio`, try to write it as ± q**n
    with integer n. If possible, return (True, sign, n); else (False, None, None).
    """
    ratio = sp.simplify(ratio)

    # Trivial case
    if ratio == 1:
        return True, 1, 0
    if ratio == -1:
        return True, -1, 0

    # Evaluate at q = 1 to get numerical ±1 factor
    try:
        r_at_1 = sp.simplify(ratio.subs(q, 1))
    except Exception:
        return False, None, None

    if not r_at_1.is_number:
        return False, None, None

    if r_at_1 not in (1, -1):
        return False, None, None

    sign = int(r_at_1)
    t = sp.simplify(ratio / sign)  # should be q**n

    if t == 1:
        return True, sign, 0
    if t == q:
        return True, sign, 1

    if isinstance(t, sp.Pow) and t.base == q and t.exp.is_integer:
        return True, sign, int(t.exp)

    # Not a pure monomial in q
    return False, None, None


def compare_up_to_min_degree(json_poly, expr, q, x):
    """
    Compare JSON polynomial and expression up to min(max_x_json, max_x_expr).

    We accept success if there exists a single factor R = ± q**n such that
    for every k in [0, max_compare]:
        coeff_json(k) = R * coeff_expr(k)

    Returns:
        exact_match (bool),
        sign_match (bool),          # overall factor -1 only
        sign_q_match (bool),        # overall factor ± q**n
        sign (int or None),
        n (int or None),
        max_compare (int),
        diffs (list of (k, c_json, c_expr)),
        max_x_json (int),
        max_x_expr (int)
    """
    max_x_json = max(json_poly.keys()) if json_poly else -sp.oo
    poly_in_x = sp.Poly(expr, x)
    max_x_expr = poly_in_x.degree()

    max_compare = min(max_x_json, max_x_expr)

    diffs = []
    ratio = None  # global ratio c_json / c_expr
    all_equal = True

    for k in range(0, max_compare + 1):
        c_json = sp.simplify(json_poly.get(k, 0))
        c_expr = sp.simplify(expr.coeff(x, k))

        if sp.simplify(c_json - c_expr) != 0:
            all_equal = False
            diffs.append((k, c_json, c_expr))

        # Determine global ratio using nonzero coefficients
        if c_json == 0 and c_expr == 0:
            # no information from this term
            continue

        if c_expr == 0 or c_json == 0:
            # can't have a single global factor if one side is zero and the other is not
            return False, False, False, None, None, max_compare, diffs, max_x_json, max_x_expr

        r = sp.simplify(c_json / c_expr)

        if ratio is None:
            ratio = r
        else:
            if sp.simplify(r - ratio) != 0:
                # ratios differ for different k → no single global factor
                return False, False, False, None, None, max_compare, diffs, max_x_json, max_x_expr

    # If everything was zero up to max_compare, treat as equal with trivial ratio 1
    if ratio is None:
        ratio = sp.Integer(1)

    # Check if ratio is ± q**n
    ok, sign, n = extract_global_sign_q_power(ratio, q)
    if not ok:
        # ratio exists but is not ± q**n
        return all_equal, False, False, None, None, max_compare, diffs, max_x_json, max_x_expr

    # At this point, we know JSON = (± q**n) * fk_terms up to max_compare
    sign_q_match = True
    sign_only = (n == 0)  # special case: just ±1

    return all_equal, sign_only, sign_q_match, sign, n, max_compare, diffs, max_x_json, max_x_expr


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Given a knot name (e.g. 10_105), load its JSON polynomial from "
            "data/<name>.json and its fk_terms from fibered_table.json, and check "
            "if they match up to the maximum x-degree of the lower-degree polynomial. "
            "Match is accepted if P_json = ± q**n * P_table for some single integer n."
        )
    )
    parser.add_argument("knot_name", help="Knot name, e.g. '10_105'")
    parser.add_argument("--data-dir", default="results",
                        help="Directory containing <knot_name>.json (default: results)")
    parser.add_argument("--table-path", default="fibered_table.json",
                        help="Path to fibered_table.json (default: fibered_table.json)")
    args = parser.parse_args()

    json_path = os.path.join(args.data_dir, f"{args.knot_name}.json")
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"JSON polynomial file not found: {json_path}")

    # Load data
    json_poly, q, x = load_json_poly(json_path)
    fk_terms_str = load_fk_terms_from_table(args.table_path, args.knot_name)
    expr = parse_string_poly(fk_terms_str, q, x)

    # Compare
    (exact_match,
     sign_only_match,
     sign_q_match,
     sign,
     n,
     max_compare,
     diffs,
     max_x_json,
     max_x_expr) = compare_up_to_min_degree(json_poly, expr, q, x)

    print(f"Knot: {args.knot_name}")
    print(f"JSON polynomial degree in x: {max_x_json}")
    print(f"fk_terms degree in x:        {max_x_expr}")
    print(f"Comparing up to x-degree:    {max_compare}")
    print()

    if exact_match:
        print("✅ Exact match up to that degree.")
    elif sign_only_match:
        print(f"✅ Match up to overall sign: JSON = ({sign}) * fk_terms (n = 0).")
    elif sign_q_match:
        print(f"✅ Match up to overall factor ± q**n: JSON = ({sign}) * q**{n} * fk_terms.")
    else:
        print("❌ No match of the form JSON = ± q**n * fk_terms up to that degree.")
        if diffs:
            print("Differences in raw coefficients (without rescaling):")
            for k, c_json, c_expr in diffs:
                print(f"  x^{k}:")
                print(f"    JSON:    {c_json}")
                print(f"    fk_terms:{c_expr}")


if __name__ == "__main__":
    main()
'''
