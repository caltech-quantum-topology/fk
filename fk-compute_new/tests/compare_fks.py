import pandas as pd
import sympy
import json
from pathlib import Path
from sympy.core.expr import Expr
from sympy.parsing.sympy_parser import parse_expr

def truncate_to_x_degree(expr: Expr, max_x_degree: int) -> Expr:
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
    x = sympy.Symbol("x")
    if expr == 0:
        return -1
    return int(sympy.Poly(expr, x, domain="EX").degree())

def str_fk_to_sympy(fk_str: object | None) -> Expr:
    q, x = sympy.symbols("q x")

    if fk_str is None:
        return sympy.Integer(0)

    s = str(fk_str).strip()
    if not s or s == "0":
        return sympy.Integer(0)

    # topology_fyi uses ^ for exponentiation; SymPy expects Python's **.
    s = s.replace("^", "**")

    try:
        expr = parse_expr(s, local_dict={"q": q, "x": x}, evaluate=True)
    except Exception as e:
        raise ValueError(f"Failed to parse fk string: {fk_str!r}") from e

    return sympy.expand(expr)

def data_file_to_sympy(fk_data_file: Path) -> Expr:
    q = sympy.Symbol("q")

    with fk_data_file.open("r", encoding="utf-8") as f:
        data = json.load(f)

    meta = data.get("metadata", {})
    num_x_vars = int(meta.get("num_x_variables", 1) or 1)

    if num_x_vars == 1:
        x_vars = (sympy.Symbol("x"),)
    else:
        x_vars = tuple(sympy.Symbol(f"x{i}") for i in range(num_x_vars))

    def _coeff_from_str(s: str) -> Expr:
        s = str(s)
        return sympy.Rational(s) if "/" in s else sympy.Integer(s)

    expr = sympy.Integer(0)
    for term in data.get("terms", []):
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


if __name__ == "__main__":

    topology_fyi_data = pd.read_json("topology_fyi_data.json")
    data_folder = Path("data")

    data_contents = {p.name for p in data_folder.glob("*[0-9].json")}
    # This will iterate over each row as (index, Series)
    for index, knot in topology_fyi_data.iterrows():
        knot_file_name = f"{knot['knotinfo_id']}.json"
        if knot_file_name in data_contents:
            print(f"Comparing {knot['knotinfo_id']}...")
            knot_file = data_folder / knot_file_name
            data_fk = data_file_to_sympy(knot_file)
            fyi_fk = str_fk_to_sympy(knot["fk_polynomial"])

            common_degree = min(max_x_degree(data_fk), max_x_degree(fyi_fk))
            data_fk = truncate_to_x_degree(data_fk, common_degree)
            fyi_fk = truncate_to_x_degree(fyi_fk, common_degree)

            if sympy.expand(data_fk - fyi_fk) != 0:
                print(f"Mismatch: {knot['knotinfo_id']}")
                print("data folder")
                print(data_fk)
                print("topology")
                print(fyi_fk)
                exit()

            print(f"{knot['knotinfo_id']} passed")
    print("All knots passed!")
