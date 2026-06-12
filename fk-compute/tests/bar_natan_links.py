"""
Bar-Natan's Z-function extended to oriented links.

The knot implementation in ``bar_natan.py`` works with one long component and
one variable.  For links we also need to know which component each arc segment
belongs to, so the crossing matrix can use the correct component variable.
"""

from __future__ import annotations

from collections import Counter
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
    """Return oriented arc cycles using the same PD successor convention as ``bar_natan.py``."""
    successors = {}
    labels = set()
    for a, b, c, d in crossings:
        successors[a] = c
        successors[d] = b
        labels.update((a, b, c, d))

    if labels != set(successors):
        missing = sorted(labels - set(successors))
        raise ValueError(f"PD code is missing successor data for labels {missing}")

    orders = []
    seen = set()
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
    arc_components = {}
    for component, labels in enumerate(component_orders):
        for label in labels:
            label = int(label)
            if label in arc_components:
                raise ValueError(f"Arc label {label} appears in more than one component")
            arc_components[label] = component
    return arc_components


def _next_arc_by_order(component_orders: Sequence[Sequence[int]]) -> dict[int, int]:
    next_arc = {}
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
    component_orders: Sequence[Sequence[int]] | None = None,
    crossing_signs: Sequence[int] | None = None,
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

    cs = []
    for index, (a, b, _c, d) in enumerate(crossings):
        if crossing_signs is None:
            if next_arc.get(d) == b:
                sign = 1
            elif next_arc.get(b) == d:
                sign = -1
            else:
                raise ValueError(
                    "Cannot infer crossing sign from component_orders; pass crossing_signs explicitly"
                )
        else:
            sign = 1 if int(crossing_signs[index]) > 0 else -1

        over_label = d if sign == 1 else b
        under_label = a
        cs.append((sign, over_label, under_label))

    max_label = max(arc_components) if arc_components else 0
    phi = [0] * max_label
    return cs, phi, arc_components


def _rotation_normalization_factor(
    cs: Sequence[RotCrossing],
    phi: Sequence[int],
    variables: Sequence[sp.Expr],
    arc_components: Mapping[int, int],
) -> sp.Expr:
    sum_phi = [0] * len(variables)
    for label, rot in enumerate(phi, start=1):
        component = arc_components.get(label)
        if component is None:
            continue
        sum_phi[component] += int(rot)

    sum_signs = [0] * len(variables)
    for sign, over_label, _ in cs:
        sum_signs[arc_components[over_label]] += int(sign)

    out = sp.Integer(1)
    for variable, phi_sum, sign_sum in zip(variables, sum_phi, sum_signs):
        out *= variable ** Rational(-(phi_sum + sign_sum), 2)
    return out


def _braid_alexander_normalization(
    braid: Sequence[int],
    variables: Sequence[sp.Expr],
) -> tuple[sp.Expr, int]:
    topology = BraidTopology(list(braid))
    open_component = int(topology.closed_strand_components[0])
    strand_counts = [0] * len(variables)
    for component in topology.closed_strand_components:
        strand_counts[int(component)] += 1

    linking_sums = [sp.Integer(0)] * len(variables)
    for sign, over_component, under_component in zip(
        topology.crossing_signs,
        topology.top_crossing_components,
        topology.bottom_crossing_components,
    ):
        over_component = int(over_component)
        under_component = int(under_component)
        if over_component == under_component:
            continue
        linking_sums[over_component] += Rational(sign, 2)
        linking_sums[under_component] += Rational(sign, 2)

    factor_out = sp.Integer(1)
    for variable, strand_count, linking_sum in zip(variables, strand_counts, linking_sums):
        factor_out *= variable ** Rational(-(-strand_count + linking_sum), 2)
    return factor_out, open_component


def _compute_z(
    cs: Sequence[RotCrossing],
    phi: Sequence[int],
    variables: Sequence[sp.Expr],
    arc_components: Mapping[int, int],
    *,
    delta_factor: sp.Expr | None = None,
    denominator_component: int | None = None,
    label_map: Mapping[int, int] | None = None,
    factor_result: bool = True,
) -> tuple[sp.Expr, sp.Expr]:
    if not cs:
        return sp.Integer(1), sp.Integer(0)

    if label_map is None:
        label_map = {}

    def _target(label: int) -> int:
        return label_map.get(label, label)

    max_label = max(
        max(_target(i), _target(j), _target(i + 1), _target(j + 1))
        for _, i, j in cs
    )
    size = max_label + 1
    a_matrix = eye(size)

    for sign, i, j in cs:
        i_next = _target(i + 1)
        j_next = _target(j + 1)
        over_variable = variables[arc_components[i]]
        under_variable = variables[arc_components[j]]
        local_t = under_variable**sign
        a_matrix[i - 1, i_next - 1] += -local_t
        if sign == 1:
            a_matrix[i - 1, j_next - 1] += over_variable - 1
        else:
            a_matrix[i - 1, j_next - 1] += (1 - over_variable) / under_variable
        a_matrix[j - 1, j_next - 1] += -1

    if delta_factor is None:
        delta_factor = _rotation_normalization_factor(cs, phi, variables, arc_components)
    denominator = sp.Integer(1)
    if denominator_component is not None:
        denominator = 1 - variables[denominator_component]

    delta = delta_factor * a_matrix.det() / denominator
    inverse_matrix = a_matrix.inv()

    def _green(row_label: int, col_label: int) -> sp.Expr:
        row_label = _target(row_label)
        col_label = _target(col_label)
        row_component = arc_components.get(row_label)
        col_component = arc_components.get(col_label)
        entry = inverse_matrix[row_label - 1, col_label - 1]
        if row_component is None or col_component is None:
            return entry
        return entry * (1 - variables[col_component]) / (1 - variables[row_component])

    def _r1(sign: int, i: int, j: int) -> sp.Expr:
        gji = _green(j, i)
        gjp_j = _green(j + 1, j)
        gj_jp = _green(j, j + 1)
        gij = _green(i, j)
        gii = _green(i, i)
        return sign * (
            gji * (gjp_j + gj_jp - gij) - gii * (gj_jp - 1) - Rational(1, 2)
        )

    rho1 = sum(_r1(sign, i, j) for sign, i, j in cs)
    rho1 -= sum(
        rot * (_green(label, label) - Rational(1, 2))
        for label, rot in enumerate(phi, start=1)
        if label < size
    )

    p0 = delta
    p1 = delta**2 * rho1
    if factor_result:
        p0 = factor(p0)
        p1 = factor(p1)
    return p0, p1


def bar_natan_link_Z(
    crossings: Sequence[Sequence[int]],
    *,
    Tnames: Sequence[sp.Expr] | None = None,
    component_orders: Sequence[Sequence[int]] | None = None,
    crossing_signs: Sequence[int] | None = None,
    factor_result: bool = True,
) -> tuple[sp.Expr, sp.Expr]:
    """
    Compute the link extension ``(P0, P1)`` for a PD diagram.

    ``Tnames`` contains one SymPy variable per component.  For one-component
    PD diagrams this delegates to the original knot rotation code, preserving
    exact compatibility with ``bar_natan_Z``.
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
            return _compute_z(cs, phi, variables, arc_components, factor_result=factor_result)
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

    return _compute_z(cs, phi, variables, arc_components, factor_result=factor_result)


def _braid_arc_components(braid: Sequence[int]) -> dict[int, int]:
    topology = BraidTopology(list(braid))
    labels = label_crossings(list(braid))
    arc_components = {}

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


def _partial_closure_label_map(
    braid: Sequence[int],
    arc_components: Mapping[int, int],
) -> dict[int, int]:
    topology = BraidTopology(list(braid))
    open_component = int(topology.closed_strand_components[0])
    labels = label_crossings(list(braid))
    counts = Counter(label for crossing in labels for label in crossing)
    terminals = sorted(label for label, count in counts.items() if count == 1)

    label_map = {}
    for index in range(0, len(terminals), 2):
        start = int(terminals[index])
        end = int(terminals[index + 1])
        if arc_components[start] == open_component:
            continue
        label_map[end] = start
    return label_map


def bar_natan_link_Z_from_braid(
    braid: Sequence[int],
    *,
    Tnames: Sequence[sp.Expr] | None = None,
    phi: Sequence[int] | None = None,
    factor_result: bool = True,
) -> tuple[sp.Expr, sp.Expr]:
    """
    Compute the link extension for the standard closure of a braid word.

    ``phi`` defaults to zero rotation on every braid arc.  Pass explicit
    rotation data if the braid closure is drawn with additional cuaps.
    """
    braid = [int(generator) for generator in braid]
    if not braid:
        return sp.Integer(1), sp.Integer(0)

    cs = [tuple(int(x) for x in crossing) for crossing in new_rot(braid)]
    arc_components = _braid_arc_components(braid)
    label_map = _partial_closure_label_map(braid, arc_components)
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

    delta_factor, denominator_component = _braid_alexander_normalization(braid, variables)
    return _compute_z(
        cs,
        rotations,
        variables,
        arc_components,
        delta_factor=delta_factor,
        denominator_component=denominator_component,
        label_map=label_map,
        factor_result=factor_result,
    )


__all__ = ["bar_natan_link_Z", "bar_natan_link_Z_from_braid"]

if __name__=="__main__":
    braid = [1,1]
    p0,p1 = bar_natan_link_Z_from_braid(braid)
    print(p0, p1)
