"""
Bar-Natan's Z-function for knot invariants.

Pure Python + SymPy translation of the Mathematica implementation
from Davide-MMR-Check.nb.

References:
    https://drorbn.net/AcademicPensieve/Talks/Nara-2308/Cars.pdf
    https://arxiv.org/pdf/1708.04853
    https://arxiv.org/pdf/2509.18456
"""

from sympy import Symbol, Rational, eye, factor
import json
import sympy as sp

T = Symbol("T")


def _positive_q(a, b, c, d, num_edges):
    """
    Determine if crossing X[a,b,c,d] is positive.

    In KnotTheory PD convention, the over-strand connects arcs b
    (position 2) and d (position 4).  A positive crossing has the
    over-strand oriented d -> b, which means b = d + 1 (mod 2n).
    """
    return (b - d) % num_edges == 1


def _rot(crossings):
    """
    Compute crossing data and rotation numbers from a planar diagram.

    Translates the Mathematica ``Rot[pd_PD]`` function.

    Parameters
    ----------
    crossings : list of (a, b, c, d) tuples
        Crossings in KnotTheory PD notation ``X[a,b,c,d]``.

    Returns
    -------
    Cs : list of (sign, i, j) triples
    phi : list of int, length 2n (rotation numbers)
    """
    n = len(crossings)
    num_edges = 2 * n
    rots = [0] * (num_edges + 1)  # 1-indexed: rots[1] .. rots[2n]

    # Build xs list: Xp(i,j) or Xm(i,j) for each crossing
    xs = []
    for a, b, c, d in crossings:
        if _positive_q(a, b, c, d, num_edges):
            xs.append(("p", d, a))   # Xp[x[[4]], x[[1]]]
        else:
            xs.append(("m", b, a))   # Xm[x[[2]], x[[1]]]

    front = [1]

    for k in range(1, num_edges + 1):
        if -k not in front:
            # Compute replacement: apply rules to ALL crossings in xs.
            # Each crossing either matches a rule (producing elements)
            # or produces nothing.
            replacement = []
            for typ, ci, cj in xs:
                if (typ == "p" and ci == k) or (typ == "m" and cj == k):
                    # Xp[k, l] | Xm[l, k]  :>  {l+1, k+1, -l}
                    l = cj if typ == "p" else ci
                    replacement.extend([l + 1, k + 1, -l])
                elif (typ == "p" and cj == k) or (typ == "m" and ci == k):
                    # Xp[l, k] | Xm[k, l]  :>  (++rots[[l]]; {-l, k+1, l+1})
                    l = ci if typ == "p" else cj
                    rots[l] += 1
                    replacement.extend([-l, k + 1, l + 1])
                # else: crossing doesn't involve k -> {} (nothing)

            # Replace every occurrence of k in front with replacement
            new_front = []
            for item in front:
                if item == k:
                    new_front.extend(replacement)
                else:
                    new_front.append(item)
            front = new_front
        else:
            # -k IS in front.
            # Cases[front, k|-k] /. {k,-k} :> --rots[[k]]
            # Only decrement if both k and -k appear, with k before -k.
            cases = [x for x in front if x == k or x == -k]
            if cases == [k, -k]:
                rots[k] -= 1

    Cs = [
        (1, ci, cj) if typ == "p" else (-1, ci, cj)
        for typ, ci, cj in xs
    ]
    phi = rots[1:]  # length 2n
    return Cs, phi


def bar_natan_Z(crossings, Tname = T):
    """
    Compute Bar-Natan's Z = (P0, P1) for a knot given in PD notation.

    This is a faithful SymPy translation of the Mathematica ``Z[K_]``
    function from Davide-MMR-Check.nb.

    Given a knot projection, Z returns ``(P0, P1)`` such that the MMR is
    ``(1 - x^{-1}) * (1/P0 - (q-1)*P1/P0^3 + O((q-1)^2))``,
    and ``P0`` is the Alexander polynomial.

    Parameters
    ----------
    crossings : list of (a, b, c, d) tuples
        Crossings in KnotTheory PD notation ``X[a,b,c,d]``.
        For example the trefoil is::

            [(1, 4, 2, 5), (3, 6, 4, 1), (5, 2, 6, 3)]

    Returns
    -------
    (P0, P1) : tuple of SymPy expressions in ``T``
    """
    T = Tname
    Cs, phi = _rot(crossings)
    n = len(Cs)
    size = 2 * n + 1

    # Build (2n+1) x (2n+1) identity matrix
    A = eye(size)

    # A[[{i,j}, {i+1,j+1}]] += {{-T^s, T^s-1}, {0, -1}}
    for s, i, j in Cs:
        Ts = T ** s
        # Mathematica 1-indexed -> Python 0-indexed
        A[i - 1, i] += -Ts
        A[i - 1, j] += Ts - 1
        # A[j-1, i] += 0   (no-op)
        A[j - 1, j] += -1

    # Delta = T^((-sum(phi) - sum(signs)) / 2) * det(A)
    sum_phi = sum(phi)
    sum_signs = sum(s for s, _, _ in Cs)
    exp = Rational(-(sum_phi + sum_signs), 2)

    Delta = T ** exp * A.det()

    G = A.inv()

    # R1(s, i, j) = s * (g_{j,i}*(g_{j+1,j} + g_{j,j+1} - g_{i,j})
    #                   - g_{i,i}*(g_{j,j+1} - 1) - 1/2)
    # where g_{a,b} = G[a-1, b-1]  (converting 1-indexed to 0-indexed)
    def _r1(s, i, j):
        gji = G[j - 1, i - 1]
        gjp_j = G[j, j - 1]          # g_{j+1, j}
        gj_jp = G[j - 1, j]          # g_{j, j+1}
        gij = G[i - 1, j - 1]        # g_{i, j}
        gii = G[i - 1, i - 1]        # g_{i, i}
        return s * (
            gji * (gjp_j + gj_jp - gij) - gii * (gj_jp - 1) - Rational(1, 2)
        )

    rho1 = sum(_r1(s, i, j) for s, i, j in Cs)
    rho1 -= sum(
        phi[k - 1] * (G[k - 1, k - 1] - Rational(1, 2))
        for k in range(1, 2 * n + 1)
    )

    P0 = factor(Delta)
    P1 = factor(Delta ** 2 * rho1)

    return P0, P1



def json_to_sympy_expr(obj, *, q_symbol="q", x_prefix="x"):
    """
    Convert your JSON structure to a SymPy expression.

    JSON schema (as in your example):
      obj["terms"] = [
        {"x": [deg0, deg1, ...], "q_terms": [{"q": qp, "c": "..."}, ...]},
        ...
      ]
      obj["metadata"]["num_x_variables"] optionally specifies number of x vars.

    Returns:
      (expr, q, xs)
        expr: SymPy expression in q and xs
        q: SymPy Symbol
        xs: tuple of SymPy Symbols (x0, x1, ...)
    """
    q = sp.Symbol(q_symbol)

    # Determine number of x variables
    nvars = None
    max_x_pow = 0
    if isinstance(obj, dict):
        nvars = obj.get("metadata", {}).get("num_x_variables", None)
        max_x_pow= obj.get("metadata", {}).get("max_x_degrees", None)
        if max_x_pow is not None:
            max_x_pow = max_x_pow[0]

    # Fallback: infer from first term if metadata absent
    if nvars is None:
        first = obj["terms"][0]["x"]
        nvars = len(first)

    xs = sp.symbols(f"{x_prefix}0:{nvars}")  # x0, x1, ..., x(nvars-1)

    expr = sp.Integer(0)

    for term in obj.get("terms", []):
        x_degs = term.get("x", [])
        if len(x_degs) != nvars:
            raise ValueError(f"Expected {nvars} x-degrees, got {len(x_degs)} in term {term}")

        # monomial in x's: Π_i x_i**deg_i
        monom_x = sp.Integer(1)
        for xi, di in zip(xs, x_degs):
            monom_x *= xi ** int(di)

        # polynomial/series part in q: Σ (Integer(c) * q**q_power)
        q_part = sp.Integer(0)
        for qt in term.get("q_terms", []):
            qp = int(qt["q"])
            c = sp.Integer(qt["c"])   # arbitrary precision from string
            q_part += c * (q ** qp)

        expr += monom_x * q_part

    return sp.simplify(expr), q, xs, max_x_pow


if __name__ == "__main__":
    import pandas as pd
    from ast import literal_eval
    ki_pd = pd.read_csv("knotinfo.csv")

    def parse_braid(s):
        return literal_eval(s.replace(";", ","))

    ki_pd = pd.read_csv(
        "knotinfo.csv",
        converters={"PD Notation": parse_braid}
    )

    P0, P1 = bar_natan_Z(ki_pd["PD Notation"][0])
    with open("data/3_1.json") as f:
        fk, q, xs, x_max = json_to_sympy_expr(json.load(f))
    print(f"P0: {P0}\nP1: {P1}\nfk: {fk}")
    Alex0series = sp.series((1-1/T)/P0, T, 0, x_max+1).removeO()
    Alex0list = sp.Poly(Alex0series,T).all_coeffs()
    print(Alex0list)
    fk0 = fk.subs(q,1)
    fk0list = sp.Poly(fk0,xs[0]).all_coeffs()
    print(fk0list)
    print(fk0list==Alex0list)
     
    Alex1series = sp.series(-(1-1/T)*P1/P0**3,T,0,x_max+1).removeO()
    print(Alex1series)
    Alex1list = sp.Poly(Alex1series,T).all_coeffs()
    print(Alex1list)
    fk1 = fk.diff(q).subs(q,1)
    fk1list = sp.Poly(fk1,xs[0]).all_coeffs()
    print(fk1list)
    print(fk1list==Alex1list)
