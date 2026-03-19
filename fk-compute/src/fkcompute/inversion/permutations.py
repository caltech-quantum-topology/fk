"""
Permutation generation and sign-diagram conversion for braid arc labels.

Ports the following Mathematica functions from ImproveSignDiagramSearch-Davide.nb:
  topstrandQ, firstAppearance, lastAppearanceQ, nextgeneratorpos,
  labelCrossings, NewRot, extendAndFilter, permsRot, permsRotClosed, permToSigns.

Arc labels are 1-indexed throughout (matching Mathematica semantics).
Crossing positions j are also 1-indexed.
"""

from __future__ import annotations

from bisect import bisect_right
from enum import IntEnum
from functools import reduce
from typing import Callable, Dict, Iterator, List, Optional, Tuple
from collections import Counter

from ..domain.braid.topology import BraidTopology


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


class _Corner(IntEnum):
    """Corner labels (Mathematica-compatible): 1..4."""

    BOTTOM_LEFT = 1
    BOTTOM_RIGHT = 2
    TOP_LEFT = 3
    TOP_RIGHT = 4


class _OccurrenceIndex:
    """Index generator occurrences for O(1)/O(log n) queries."""

    def __init__(self, braid: List[int]):
        positions: Dict[int, List[int]] = {}
        for j, g in enumerate(braid, 1):
            positions.setdefault(abs(g), []).append(j)
        self.positions = positions
        self.first = {gen: pos[0] for gen, pos in positions.items()}
        self.last = {gen: pos[-1] for gen, pos in positions.items()}

    def first_appearance(self, gen: int) -> int:
        g = abs(gen)
        if g not in self.first:
            raise ValueError(f"Generator {gen} not found in braid")
        return self.first[g]

    def is_last_appearance(self, gen: int, j: int) -> bool:
        return self.last[abs(gen)] == j

    def next_pos(self, gen: int, j: int) -> Optional[int]:
        pos = self.positions.get(abs(gen))
        if not pos:
            return None
        idx = bisect_right(pos, j)
        return None if idx >= len(pos) else pos[idx]


def _paired_corner(corner: _Corner) -> _Corner:
    """Move through a crossing: BL->TR, BR->TL."""
    if corner == _Corner.BOTTOM_LEFT:
        return _Corner.TOP_RIGHT
    if corner == _Corner.BOTTOM_RIGHT:
        return _Corner.TOP_LEFT
    raise ValueError(f"Invalid entry corner: {corner}")


def _is_topstrand(
    *,
    occ: _OccurrenceIndex,
    j: int,
    gen_j: int,
    corner: _Corner,
    max_gen: int,
) -> bool:
    """Return True if the (j, corner) position is on the top strand."""
    if occ.last[gen_j] != j:
        return False
    if corner == _Corner.TOP_LEFT and gen_j > 1:
        if occ.last.get(gen_j - 1, 0) > j:
            return False
    if corner == _Corner.TOP_RIGHT and gen_j < max_gen:
        if occ.last.get(gen_j + 1, 0) > j:
            return False
    return True


def _next_crossing(
    *,
    occ: _OccurrenceIndex,
    j: int,
    corner: _Corner,
    gen_j: int,
    max_gen: int,
) -> int:
    """Return the next 1-indexed crossing position j_next."""
    if _is_topstrand(occ=occ, j=j, gen_j=gen_j, corner=corner, max_gen=max_gen):
        # Top strand: wrap around to the first appearance of adjacent generators.
        if corner == _Corner.TOP_LEFT:
            if gen_j > 1:
                return min(
                    occ.first_appearance(gen_j),
                    occ.first_appearance(gen_j - 1),
                )
            return occ.first_appearance(gen_j)

        # corner == TOP_RIGHT
        if gen_j < max_gen:
            return min(
                occ.first_appearance(gen_j),
                occ.first_appearance(gen_j + 1),
            )
        return occ.first_appearance(gen_j)

    # Not top strand: find next crossing of same or adjacent generator.
    candidates: List[int] = []
    if not occ.is_last_appearance(gen_j, j):
        pos = occ.next_pos(gen_j, j)
        if pos is not None:
            candidates.append(pos)
    if corner == _Corner.TOP_LEFT and gen_j > 1:
        pos = occ.next_pos(gen_j - 1, j)
        if pos is not None:
            candidates.append(pos)
    if corner == _Corner.TOP_RIGHT and gen_j < max_gen:
        pos = occ.next_pos(gen_j + 1, j)
        if pos is not None:
            candidates.append(pos)

    if not candidates:
        raise ValueError("No next crossing found (unexpected braid structure)")
    return min(candidates)


def _entry_corner(exit_corner: _Corner, same_generator: bool) -> _Corner:
    """Return entry corner (BL/BR) at the next crossing."""
    if exit_corner == _Corner.TOP_LEFT:
        return _Corner.BOTTOM_LEFT if same_generator else _Corner.BOTTOM_RIGHT
    if exit_corner == _Corner.TOP_RIGHT:
        return _Corner.BOTTOM_RIGHT if same_generator else _Corner.BOTTOM_LEFT
    raise ValueError(f"Invalid exit corner: {exit_corner}")


# ---------------------------------------------------------------------------
# Core labelling
# ---------------------------------------------------------------------------


def label_crossings(braid: List[int]) -> List[List[int]]:
    """
    Label arc segments at every corner of every crossing.

    Returns an n×4 matrix (list of lists) where entry [i][c] is the arc-segment
    label at corner c+1 of crossing i+1.

    Labels are 1-indexed. For a knot (single component) they run from 1 to
    2n+1 where n = len(braid). For links with multiple components the range
    extends to 2n+c, where c is the number of component cycles induced by the
    entry-corner successor map.

    Direct port of Mathematica ``labelCrossings``.
    """
    n = len(braid)
    if n == 0:
        return []

    max_gen = max(abs(g) for g in braid)
    occ = _OccurrenceIndex(braid)

    # 0-indexed rows, 0-indexed corners (corner c ↔ index c-1)
    labels_corners: List[List[int]] = [[0, 0, 0, 0] for _ in range(n)]

    def _find_next_start() -> Optional[Tuple[int, _Corner]]:
        # Prefer entry corners (bottom corners) only.
        for row_idx in range(n):
            if labels_corners[row_idx][0] == 0:
                return row_idx + 1, _Corner.BOTTOM_LEFT
            if labels_corners[row_idx][1] == 0:
                return row_idx + 1, _Corner.BOTTOM_RIGHT
        return None

    def _walk_component(start_j: int, start_corner: _Corner, k_start: int) -> int:
        """Walk one closed cycle, starting from an already-labelled entry corner."""
        j = start_j
        corner = start_corner
        k = k_start

        # Each loop iteration assigns one new arc-segment label (k += 1) and
        # propagates it to two corners: the exit (top) corner at j and the entry
        # (bottom) corner at j_next.
        steps = 0
        while True:
            steps += 1
            if steps > 2 * n + 5:
                raise ValueError("label_crossings exceeded step limit (cycle walk)")

            # Each step creates a new arc label.
            k += 1

            # Move through the crossing: BL->TR, BR->TL (same crossing, top side)
            exit_corner = _paired_corner(corner)
            labels_corners[j - 1][exit_corner.value - 1] = k

            # Determine next crossing j_next and entry corner
            gen_j = abs(braid[j - 1])
            j_next = _next_crossing(
                occ=occ,
                j=j,
                corner=exit_corner,
                gen_j=gen_j,
                max_gen=max_gen,
            )

            same_generator = abs(braid[j - 1]) == abs(braid[j_next - 1])
            entry_corner = _entry_corner(exit_corner, same_generator)

            # If the cycle closes, stop without overwriting the starting label.
            if j_next == start_j and entry_corner == start_corner:
                break

            if labels_corners[j_next - 1][entry_corner.value - 1] != 0:
                raise ValueError(
                    "Cycle walk hit an already-labelled entry corner (unexpected)"
                )

            j = j_next
            corner = entry_corner
            labels_corners[j - 1][corner.value - 1] = k

        return k

    # Walk each connected component cycle. The original Mathematica function
    # starts from the first appearance of generator 1 at BL.
    k = 1
    j_first = occ.first_appearance(1)  # 1-indexed
    corner_first = _Corner.BOTTOM_LEFT
    labels_corners[j_first - 1][corner_first.value - 1] = k
    k = _walk_component(j_first, corner_first, k)

    # For links, the arc graph has multiple disjoint cycles. Continue from the
    # next unlabelled entry corner until we've used all available arc labels.
    while True:
        nxt = _find_next_start()
        if nxt is None:
            break
        start_j, start_corner = nxt
        k += 1
        labels_corners[start_j - 1][start_corner.value - 1] = k
        k = _walk_component(start_j, start_corner, k)

    # After walking all cycles, all corners should typically be labelled.
    # Keep a defensive fill for any remaining zeros.
    last_label = k
    for row in labels_corners:
        for m in range(4):
            if row[m] == 0:
                row[m] = last_label

    return labels_corners


def new_rot(braid: List[int]) -> List[List[int]]:
    """
    Return a list of ``[sign, overstrand_label, understrand_label]`` per crossing.

    For a positive crossing the overstrand labels are corners 1,2; for a
    negative crossing they are corners 2,1 (swapped).

    Direct port of Mathematica ``NewRot``.
    """
    lc = label_crossings(braid)
    result: List[List[int]] = []
    for n_idx in range(len(braid)):
        sign = 1 if braid[n_idx] > 0 else -1
        if sign == 1:
            result.append([sign, lc[n_idx][0], lc[n_idx][1]])  # corners 1, 2
        else:
            result.append([sign, lc[n_idx][1], lc[n_idx][0]])  # corners 2, 1
    return result


# ---------------------------------------------------------------------------
# Permutation generation helpers
# ---------------------------------------------------------------------------


def _extend_and_filter(
    current: List[List[Tuple[int, int]]],
    next_options: List[Tuple[int, int]],
) -> List[List[Tuple[int, int]]]:
    """
    Cross-product extend and filter.

    Appends each option in *next_options* to each partial permutation in
    *current*, keeping only those where all destination values are distinct.

    Port of Mathematica ``extendAndFilter``.
    """
    extended = [partial + [rule] for partial in current for rule in next_options]
    return [x for x in extended if len({dst for _, dst in x}) == len(x)]


def _build_perm_options(braid: List[int], closed: bool) -> List[List[Tuple[int, int]]]:
    """Build the sorted list of option groups used by perms_rot / perms_rot_closed."""
    rot_vec = new_rot(braid)
    n = len(braid)
    perm_options: List[List[Tuple[int, int]]] = []

    for entry in rot_vec:
        _, label_a, label_b = entry
        perm_options.append(
            [(label_a, label_a), (label_a, label_a + 1), (label_a, label_b + 1)]
        )
        perm_options.append([(label_b, label_b), (label_b, label_b + 1)])

    last = 2 * n + 1
    if closed:
        labels = label_crossings(braid)
        counts = Counter(x for crossing in labels for x in crossing)
        terminals = sorted(x for x, c in counts.items() if c == 1)
        for i in range(0, len(terminals), 2):
            perm_options.append(
                [(terminals[i + 1], terminals[i]), (terminals[i + 1], terminals[i + 1])]
            )
    else:
        perm_options.append([(last, last)])

    perm_options.sort(key=lambda group: group[0][0])
    return perm_options


def _fold_to_perms(
    perm_options: List[List[Tuple[int, int]]], total: int
) -> List[List[int]]:
    """Fold option groups into valid permutation lists."""
    partial_perms = reduce(_extend_and_filter, perm_options, [[]])
    result: List[List[int]] = [[dst for _, dst in pp] for pp in partial_perms]
    return result


def _iter_fold_to_perms(perm_options: List[List[Tuple[int, int]]]) -> Iterator[List[int]]:
    """Yield valid permutation lists lazily.

    This is a generator equivalent of :func:`_fold_to_perms` that avoids
    constructing the full cross-product in memory.
    """

    def _recurse(i: int, used_dsts: set[int], out: List[int]) -> Iterator[List[int]]:
        if i >= len(perm_options):
            # out is already in the correct order (sorted by src label).
            yield list(out)
            return

        for _, dst in perm_options[i]:
            if dst in used_dsts:
                continue
            used_dsts.add(dst)
            out.append(dst)
            yield from _recurse(i + 1, used_dsts, out)
            out.pop()
            used_dsts.remove(dst)

    yield from _recurse(0, set(), [])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def perms_rot(braid: List[int]) -> List[List[int]]:
    """
    Generate all valid arc-label permutations for the *open* braid.

    Each permutation is a 1-indexed list of length 2n+1 (n = len(braid))
    where entry i is the image of arc label i under the permutation.

    The final arc 2n+1 maps only to itself (open-strand convention).

    Port of Mathematica ``permsRot``.
    """
    n = len(braid)
    opts = _build_perm_options(braid, closed=False)
    return list(_iter_fold_to_perms(opts))


def perms_rot_closed(braid: List[int]) -> List[List[int]]:
    """
    Generate all valid arc-label permutations for the *closed* braid.

    Same as :func:`perms_rot` except the final arc 2n+1 may also map to 1
    (closure allows the open strand to loop back).

    Port of Mathematica ``permsRotClosed``.
    """
    n = len(braid)
    opts = _build_perm_options(braid, closed=True)
    return list(_iter_fold_to_perms(opts))


def iter_perms_rot(braid: List[int]) -> Iterator[List[int]]:
    """Iterate valid arc-label permutations for the *open* braid."""
    opts = _build_perm_options(braid, closed=False)
    return _iter_fold_to_perms(opts)


def iter_perms_rot_closed(braid: List[int]) -> Iterator[List[int]]:
    """Iterate valid arc-label permutations for the *closed* braid."""
    opts = _build_perm_options(braid, closed=True)
    return _iter_fold_to_perms(opts)


def perm_to_signs_metadata(braid: List[int]) -> tuple[List[int], Dict[int, int], int]:
    """Precompute metadata needed for fast perm -> signs conversion."""
    rot = new_rot(braid)
    topo = BraidTopology(braid)
    arc_comp: Dict[int, int] = {}
    arc_labels: List[int] = []
    for j, (crossing_sign, label_a, label_b) in enumerate(rot):
        top_c = topo.top_crossing_components[j]
        bottom_c = topo.bottom_crossing_components[j]
        if top_c is None or bottom_c is None:
            raise ValueError("Unexpected missing component assignment")

        if crossing_sign == 1:
            arc_comp[label_a] = top_c
            arc_comp[label_b] = bottom_c
        else:
            arc_comp[label_a] = bottom_c
            arc_comp[label_b] = top_c
        arc_labels.append(label_a)
        arc_labels.append(label_b)

    arc_labels = sorted(set(arc_labels))
    if len(arc_labels) != 2 * len(braid):
        raise ValueError("Unexpected number of arc labels for sign assignment")
    return arc_labels, arc_comp, topo.n_components


def compile_perm_to_signs(braid: List[int]) -> Callable[[List[int]], Dict[int, List[int]]]:
    """Return a fast converter perm -> strand_signs for this braid."""
    arc_labels, arc_comp, n_components = perm_to_signs_metadata(braid)

    def _convert(perm: List[int]) -> Dict[int, List[int]]:
        if not arc_labels:
            return {c: [] for c in range(n_components)}
        if max(arc_labels) > len(perm):
            raise ValueError(
                "Permutation too short for braid arc labels: "
                f"need >= {max(arc_labels)}, got {len(perm)}"
            )
        result: Dict[int, List[int]] = {c: [] for c in range(n_components)}
        for label in arc_labels:
            result[arc_comp[label]].append(1 if perm[label - 1] == label else -1)
        return result

    return _convert


def perm_to_signs(perm: List[int], braid: List[int]) -> Dict[int, List[int]]:
    """
    Convert a 1-indexed permutation to a dict mapping component index → list of ±1 signs.

    Arc labels used by :func:`new_rot` are assigned to link components using
    ``BraidTopology``. Signs within each component are collected in arc-label
    order (ascending).

    Port of Mathematica ``permToSigns``, extended to multi-component links.
    """
    rot = new_rot(braid)
    topo = BraidTopology(braid)
    arc_comp: Dict[int, int] = {}
    arc_labels: List[int] = []
    for j, (crossing_sign, label_a, label_b) in enumerate(rot):
        top_c = topo.top_crossing_components[j]
        bottom_c = topo.bottom_crossing_components[j]
        if top_c is None or bottom_c is None:
            raise ValueError("Unexpected missing component assignment")

        # new_rot returns labels as (over, under). For positive crossings the
        # overstrand is the top strand; for negative crossings it's the bottom.
        if crossing_sign == 1:
            arc_comp[label_a] = top_c
            arc_comp[label_b] = bottom_c
        else:
            arc_comp[label_a] = bottom_c
            arc_comp[label_b] = top_c
        arc_labels.append(label_a)
        arc_labels.append(label_b)

    # The sign assignment is defined on the 2n arc labels appearing in new_rot.
    # For links, these labels need not be exactly 1..2n.
    arc_labels = sorted(set(arc_labels))
    if len(arc_labels) != 2 * len(braid):
        raise ValueError("Unexpected number of arc labels for sign assignment")
    if not arc_labels:
        return {c: [] for c in range(topo.n_components)}
    if max(arc_labels) > len(perm):
        raise ValueError(
            "Permutation too short for braid arc labels: "
            f"need >= {max(arc_labels)}, got {len(perm)}"
        )

    result: Dict[int, List[int]] = {c: [] for c in range(topo.n_components)}
    for label in arc_labels:
        result[arc_comp[label]].append(1 if perm[label - 1] == label else -1)
    return result
