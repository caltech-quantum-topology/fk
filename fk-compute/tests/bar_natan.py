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


def bar_natan_Z(crossings):
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

    for _, row in ki_pd.iterrows():
        print(f"Knot: {row['Name']}")
        print(f"P0, P1: {bar_natan_Z(row['PD Notation'])}")
