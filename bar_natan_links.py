"""
Bar-Natan's Z-function extended to oriented links.

This module keeps the matrix formula used by ``bar_natan.py`` and replaces the
single variable ``T`` by a component variable.  At a crossing, the local
``T**sign`` factor is taken from the over-strand component.  For a one-component
diagram this specializes to ``bar_natan.bar_natan_Z``.

The primary entry point is ``bar_natan_link_Z`` for PD data.  A braid-closure
helper, ``bar_natan_link_Z_from_braid``, is included because most link data in
this test tree is stored as braid words.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp
from sympy import Rational, eye, factor

from bar_natan import T as DEFAULT_T
from bar_natan import _rot as _knot_rot
from fkcompute.domain.braid.topology import BraidTopology
from fkcompute.inversion.permutations import label_crossings, new_rot


Crossing = tuple[int, int, int, int]
RotCrossing = tuple[int, int, int]


def _as_crossings(crossings: Sequence[Sequence[int]]) -> list[Crossing]:
    out = []
    for crossing in crossings:
        if len(crossing) != 4:
            raise ValueError(f"Expected PD crossing with four labels, got {crossing!r}")
        out.append(tuple(int(x) for x in crossing))
    return out


def _component_orders_from_pd(crossings: Sequence[Crossing]) -> list[list[int]]:
    """Return oriented arc cycles using the same PD successor convention as bar_natan.py."""
    successors: dict[int, int] = {}
    labels: set[int] = set()
    for a, b, c, d in crossings:
        successors[a] = c
        successors[d] = b
        labels.update((a, b, c, d))

    if labels != set(successors):
        missing = sorted(labels - set(successors))
        raise ValueError(f"PD code is missing successor data for labels {missing}")

    orders: list[list[int]] = []
    seen: set[int] = set()
    for start in sorted(labels):
        if start in seen:
            continue
        order = [start]
        seen.add(start)
        current = successors[start]
        while current != start:
            if current in seen:
                raise ValueError("PD successor graph is not a disjoint union of cycles")
            order.append(current)
            seen.add(current)
            current = successors[current]
        orders.append(order)
    return orders


def _arc_components_from_orders(component_orders: Sequence[Sequence[int]]) -> dict[int, int]:
    arc_components: dict[int, int] = {}
    for component, labels in enumerate(component_orders):
        for label in labels:
            label = int(label)
            if label in arc_components:
                raise ValueError(f"Arc label {label} appears in more than one component")
            arc_components[label] = component
    return arc_components


def _next_arc_by_order(component_orders: Sequence[Sequence[int]]) -> dict[int, int]:
    next_arc: dict[int, int] = {}
    for labels in component_orders:
        labels = [int(label) for label in labels]
        if not labels:
            continue
        for index, label in enumerate(labels):
            next_arc[label] = labels[(index + 1) % len(labels)]
    return next_arc


def _rot_link_from_pd(
    crossings: Sequence[Crossing],
    *,
    component_orders: Sequence[Sequence[int]] | None,
    crossing_signs: Sequence[int] | None,
) -> tuple[list[RotCrossing], list[int], dict[int, int]]:
    """
    Compute Bar-Natan crossing triples for a link PD.

    ``component_orders`` gives the oriented cyclic order of arc labels on each
    component.  If omitted it is inferred from the PD successor convention
    ``a -> c`` and ``d -> b``.  ``crossing_signs`` may be supplied when the PD
    labelling is not component-sequential enough to infer signs from arc order.
    """
    if component_orders is None:
        component_orders = _component_orders_from_pd(crossings)
    arc_components = _arc_components_from_orders(component_orders)
    next_arc = _next_arc_by_order(component_orders)

    if crossing_signs is not None and len(crossing_signs) != len(crossings):
        raise ValueError("crossing_signs must have one entry per crossing")

    cs: list[RotCrossing] = []
    for index, (a, b, _c, d) in enumerate(crossings):
        if crossing_signs is None:
            if next_arc.get(d) == b:
                sign = 1
            elif next_arc.get(b) == d:
                sign = -1
            else:
                raise ValueError(
                    "Cannot infer crossing sign from component_orders; "
                    "pass crossing_signs explicitly"
                )
        else:
            sign = 1 if int(crossing_signs[index]) > 0 else -1

        # Matches bar_natan._rot: Xp[d, a] and Xm[b, a].
        over_label = d if sign == 1 else b
        under_label = a
        cs.append((sign, over_label, under_label))

    max_label = max(arc_components) if arc_components else 0
    phi = [0] * max_label
    return cs, phi, arc_components


def _normalization_factor(
    cs: Sequence[RotCrossing],
    phi: Sequence[int],
    variables: Sequence[sp.Expr],
    arc_components: Mapping[int, int],
) -> sp.Expr:
    sum_phi = [0] * len(variables)
    for label, rot in enumerate(phi, start=1):
        component = arc_components.get(label)
        if component is not None:
            sum_phi[component] += int(rot)

    sum_signs = [0] * len(variables)
    for sign, over_label, _ in cs:
        sum_signs[arc_components[over_label]] += int(sign)

    out = sp.Integer(1)
    for variable, phi_sum, sign_sum in zip(variables, sum_phi, sum_signs):
        out *= variable ** Rational(-(phi_sum + sign_sum), 2)
    return out


def _compute_z(
    cs: Sequence[RotCrossing],
    phi: Sequence[int],
    variables: Sequence[sp.Expr],
    arc_components: Mapping[int, int],
) -> tuple[sp.Expr, sp.Expr]:
    if not cs:
        return sp.Integer(1), sp.Integer(0)

    max_label = max(max(i, j) for _, i, j in cs)
    size = max_label + 1
    a_matrix = eye(size)

    for sign, i, j in cs:
        variable = variables[arc_components[i]]
        local_t = variable ** sign
        a_matrix[i - 1, i] += -local_t
        a_matrix[i - 1, j] += local_t - 1
        a_matrix[j - 1, j] += -1

    delta = _normalization_factor(cs, phi, variables, arc_components) * a_matrix.det()
    green = a_matrix.inv()

    def _r1(sign: int, i: int, j: int) -> sp.Expr:
        gji = green[j - 1, i - 1]
        gjp_j = green[j, j - 1]
        gj_jp = green[j - 1, j]
        gij = green[i - 1, j - 1]
        gii = green[i - 1, i - 1]
        return sign * (
            gji * (gjp_j + gj_jp - gij)
            - gii * (gj_jp - 1)
            - Rational(1, 2)
        )

    rho1 = sum(_r1(sign, i, j) for sign, i, j in cs)
    rho1 -= sum(
        rot * (green[label - 1, label - 1] - Rational(1, 2))
        for label, rot in enumerate(phi, start=1)
        if label < size
    )

    return factor(delta), factor(delta**2 * rho1)


def bar_natan_link_Z(
    crossings: Sequence[Sequence[int]],
    *,
    Tnames: Sequence[sp.Expr] | None = None,
    component_orders: Sequence[Sequence[int]] | None = None,
    crossing_signs: Sequence[int] | None = None,
) -> tuple[sp.Expr, sp.Expr]:
    """
    Compute the link extension ``(P0, P1)`` for a PD diagram.

    Parameters
    ----------
    crossings
        PD crossings ``X[a,b,c,d]``.
    Tnames
        One SymPy variable per link component.  If omitted, variables are named
        ``T0, T1, ...``.  For one component, ``DEFAULT_T`` is used.
    component_orders
        Oriented cyclic arc labels for each component.  Supplying this is the
        safest way to use arbitrary link PD labels.
    crossing_signs
        Optional crossing signs.  Use this when signs cannot be inferred from
        component order.
    """
    crossings = _as_crossings(crossings)

    if component_orders is None and crossing_signs is None:
        inferred_orders = _component_orders_from_pd(crossings)
        if len(inferred_orders) == 1:
            variables = [DEFAULT_T] if Tnames is None else list(Tnames)
            if len(variables) != 1:
                raise ValueError(f"Expected 1 component variable, got {len(variables)}")
            cs, phi = _knot_rot(crossings)
            arc_components = {label: 0 for label in range(1, len(phi) + 1)}
            return _compute_z(cs, phi, variables, arc_components)
        component_orders = inferred_orders

    cs, phi, arc_components = _rot_link_from_pd(
        crossings,
        component_orders=component_orders,
        crossing_signs=crossing_signs,
    )
    n_components = max(arc_components.values(), default=-1) + 1
    if Tnames is None:
        variables = [DEFAULT_T] if n_components == 1 else list(sp.symbols(f"T0:{n_components}"))
    else:
        variables = list(Tnames)
        if len(variables) != n_components:
            raise ValueError(f"Expected {n_components} component variables, got {len(variables)}")

    return _compute_z(cs, phi, variables, arc_components)


def _braid_arc_components(braid: Sequence[int]) -> dict[int, int]:
    topology = BraidTopology(list(braid))
    labels = label_crossings(list(braid))
    arc_components: dict[int, int] = {}

    for index, row in enumerate(labels):
        top_component = topology.top_crossing_components[index]
        bottom_component = topology.bottom_crossing_components[index]
        if top_component is None or bottom_component is None:
            raise ValueError("Unexpected missing braid component assignment")

        for label in (row[0], row[3]):
            arc_components[int(label)] = int(top_component)
        for label in (row[1], row[2]):
            arc_components[int(label)] = int(bottom_component)

    return arc_components


def bar_natan_link_Z_from_braid(
    braid: Sequence[int],
    *,
    Tnames: Sequence[sp.Expr] | None = None,
    phi: Sequence[int] | None = None,
) -> tuple[sp.Expr, sp.Expr]:
    """
    Compute the link extension for the standard closure of a braid word.

    ``phi`` defaults to zero rotation on every braid arc.  Pass explicit
    rotation data if the braid closure has been drawn with additional cuaps.
    """
    braid = [int(generator) for generator in braid]
    if not braid:
        return sp.Integer(1), sp.Integer(0)

    cs = [tuple(int(x) for x in crossing) for crossing in new_rot(braid)]
    arc_components = _braid_arc_components(braid)
    n_components = max(arc_components.values(), default=-1) + 1

    if Tnames is None:
        variables = [DEFAULT_T] if n_components == 1 else list(sp.symbols(f"T0:{n_components}"))
    else:
        variables = list(Tnames)
        if len(variables) != n_components:
            raise ValueError(f"Expected {n_components} component variables, got {len(variables)}")

    max_label = max(arc_components) if arc_components else 0
    rotations = [0] * max_label if phi is None else [int(x) for x in phi]
    if len(rotations) < max_label:
        raise ValueError(f"Expected phi to have at least {max_label} entries")

    return _compute_z(cs, rotations, variables, arc_components)


__all__ = [
    "bar_natan_link_Z",
    "bar_natan_link_Z_from_braid",
]
