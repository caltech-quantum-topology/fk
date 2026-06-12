import csv
import json
from ast import literal_eval
from pathlib import Path

import pytest
import sympy as sp
from bar_natan import bar_natan_Z


def load_fk_json(file_path: str) -> sp.Expr:
    """Load FK polynomial from JSON file and return as sympy expression.

    Args:
        file_path: Path to JSON file containing FK data.

    Returns:
        Sympy expression in variables x and q.
    """
    x, q = sp.symbols("x q")

    with open(file_path, "r") as f:
        data = json.load(f)

    result = sp.Integer(0)

    for term in data["terms"]:
        x_exp = term["x"][0]
        q_part = sp.Integer(0)

        for q_term in term["q_terms"]:
            q_exp = q_term["q"]
            coeff = sp.Integer(int(q_term["c"]))
            q_part += coeff * q**q_exp

        result += q_part * x**x_exp

    return sp.expand(result)


def braid_to_pd(word):
    """
    Convert an n-strand braid word into a PD code for its standard closure.

    word: list of ints, +i = sigma_i, -i = sigma_i^{-1}, 1 <= i < n
    returns: list of 4-tuples (a,b,c,d) meaning X[a,b,c,d]
    """
    n = max(map(abs, word)) + 1
    if n < 2:
        raise ValueError("n must be >= 2")
    for g in word:
        i = abs(g)
        if i < 1 or i >= n:
            raise ValueError(f"generator {g} out of range for n={n}")

    # Current arc label at each strand position (top endpoints labeled 1..n)
    cur = list(range(1, n + 1))
    next_arc = n + 1

    pd = []
    bottoms = [None] * n  # arc labels at bottom endpoints

    def new_arc():
        nonlocal next_arc
        a = next_arc
        next_arc += 1
        return a

    # Each crossing creates two new downward arcs (one per strand position).
    # PD convention used here:
    #  - For a positive generator sigma_i (left over right):
    #      X[left_in, right_in, right_out, left_out]
    #  - For a negative generator sigma_i^{-1} (right over left):
    #      X[right_in, left_in, left_out, right_out]
    #
    # This is a consistent local rule; if you need a different PD ordering,
    # swap these tuple constructions.
    for g in word:
        i = abs(g) - 1  # zero-based strand index for (i, i+1)
        left_in, right_in = cur[i], cur[i + 1]
        left_out, right_out = new_arc(), new_arc()

        if g > 0:
            pd.append((right_in, left_out, right_out, left_in))
        else:
            pd.append((left_in, right_out, left_out, right_in))

        # After the crossing, strands swap positions in the braid
        # (because the braid generators permute endpoints).
        cur[i], cur[i + 1] = right_out, left_out

    # Record the arc labels at the bottom of each strand position
    for k in range(n):
        bottoms[k] = cur[k]

    # Close the braid: connect bottom position k to top position k.
    # This means "identify" arc label bottoms[k] with top label (k+1).
    # We'll implement identification by renaming via a union-find-like map.
    parent = {a: a for a in range(1, next_arc)}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for k in range(n):
        union(bottoms[k], k + 1)

    # Apply identifications to PD labels and then relabel to 1..m consecutively.
    pd2 = []
    for (a, b, c, d) in pd:
        pd2.append((find(a), find(b), find(c), find(d)))

    # Compress labels to consecutive integers
    uniq = {}
    out = []
    nxt = 1
    for quad in pd2:
        q = []
        for x in quad:
            if x not in uniq:
                uniq[x] = nxt
                nxt += 1
            q.append(uniq[x])
        out.append(tuple(q))

    return out

def normalize_pd(crossings):
    """Relabel a PD code so arcs are numbered 1, 2, ..., 2n sequentially
    along the knot orientation.

    For each crossing X[a,b,c,d]: succ[a]=c (understrand) and succ[d]=b (overstrand).
    """
    succ = {}
    for a, b, c, d in crossings:
        succ[a] = c
        succ[d] = b

    # Trace the knot component starting from arc 1
    order = [1]
    cur = succ[1]
    while cur != 1:
        order.append(cur)
        cur = succ[cur]

    # Build relabeling: old label -> new sequential label
    relabel = {old: new for new, old in enumerate(order, start=1)}

    return [[relabel[a], relabel[b], relabel[c], relabel[d]]
            for a, b, c, d in crossings]


def load_braid(file_path):
    with open(file_path, "r") as f:
        data = json.load(f)
    return data["metadata"]["braid"]


def _fk_coeffs_q1(data, x_max):
    """Coefficients of fk(x, q=1), ordered from degree x_max down to 0."""
    by_degree = {}
    for term in data["terms"]:
        d = term["x"][0]
        by_degree[d] = sum(int(qt["c"]) for qt in term["q_terms"])
    return [by_degree.get(d, 0) for d in range(x_max, -1, -1)]


def _fk_deriv_coeffs_q1(data, x_max):
    """Coefficients of (∂fk/∂q)|_{q=1}, ordered from degree x_max down to 0."""
    by_degree = {}
    for term in data["terms"]:
        d = term["x"][0]
        by_degree[d] = sum(int(qt["q"]) * int(qt["c"]) for qt in term["q_terms"])
    return [by_degree.get(d, 0) for d in range(x_max, -1, -1)]


def _load_test_cases():
    tests_dir = Path(__file__).parent

    # Parse knotinfo.csv → {name: pd_crossings}
    knotinfo = {}
    with open(tests_dir / "knotinfo.csv", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["Name"]
            pd_cell = row["PD Notation"]
            knotinfo[name] = literal_eval(pd_cell.replace(";", ","))

    cases = []
    for json_file in sorted((tests_dir / "data").glob("*.json")):
        if "ilp" in json_file.name or "inversion" in json_file.name:
            continue
        with open(json_file) as f:
            fk_data = json.load(f)
        if fk_data.get("metadata", {}).get("components", 1) != 1:
            continue
        knot_name = json_file.stem
        if knot_name in knotinfo:
            pd_crossings = knotinfo[knot_name]
        else:
            braid = fk_data["metadata"]["braid"]
            pd_crossings = normalize_pd(braid_to_pd(braid))
        cases.append((knot_name, pd_crossings, fk_data))

    return cases


_TEST_CASES = _load_test_cases()

@pytest.mark.parametrize(
    "knot_name,pd_crossings,fk_data",
    _TEST_CASES,
    ids=[c[0] for c in _TEST_CASES],
)
def test_bar_natan_p0_p1(knot_name, pd_crossings, fk_data):
    T = sp.Symbol("T")
    P0, P1 = bar_natan_Z(pd_crossings, Tname=T)
    x_max = fk_data["metadata"]["max_x_degrees"][0]

    # Normalize P0: require the leading coefficient of the numerator polynomial to
    # be positive (Alexander polynomial is defined up to ±T^k; knotinfo PD may give
    # either sign depending on arc-labeling orientation).
    numer = sp.Poly(sp.numer(sp.together(P0)), T)
    if numer.LC() < 0:
        P0 = -P0

    # P0 check (Alexander polynomial)
    alex0_poly = sp.Poly(sp.series((1 - 1/T)/P0, T, 0, x_max + 1).removeO(), T)
    assert [int(alex0_poly.nth(d)) for d in range(x_max, -1, -1)] == _fk_coeffs_q1(fk_data, x_max), \
        f"{knot_name}: P0 mismatch"

    # P1 check (first-order correction; P1 is unaffected by the P0 sign flip since
    # it is computed from Delta^2, which is sign-invariant)
    alex1_poly = sp.Poly(sp.series(-(1 - 1/T)*P1/P0**3, T, 0, x_max + 1).removeO(), T)
    assert [int(alex1_poly.nth(d)) for d in range(x_max, -1, -1)] == _fk_deriv_coeffs_q1(fk_data, x_max), \
        f"{knot_name}: P1 mismatch"


if __name__ == "__main__":
    T = sp.Symbol("T")
    passed = 0
    failed = 0
    for knot_name, pd_crossings, fk_data in _load_test_cases():
        x_max = fk_data["metadata"]["max_x_degrees"][0]
        P0, P1 = bar_natan_Z(pd_crossings, Tname=T)
        numer = sp.Poly(sp.numer(sp.together(P0)), T)
        if numer.LC() < 0:
            P0 = -P0
        alex0 = sp.Poly(sp.series((1 - 1/T)/P0, T, 0, x_max + 1).removeO(), T)
        ok0 = [int(alex0.nth(d)) for d in range(x_max, -1, -1)] == _fk_coeffs_q1(fk_data, x_max)
        alex1 = sp.Poly(sp.series(-(1 - 1/T)*P1/P0**3, T, 0, x_max + 1).removeO(), T)
        ok1 = [int(alex1.nth(d)) for d in range(x_max, -1, -1)] == _fk_deriv_coeffs_q1(fk_data, x_max)
        status = "PASS" if (ok0 and ok1) else f"FAIL (P0={'ok' if ok0 else 'FAIL'}, P1={'ok' if ok1 else 'FAIL'})"
        print(f"{knot_name}: {status}")
        if ok0 and ok1:
            passed += 1
        else:
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
