"""
Braid variant generation for sign assignment search.

This module provides functions for generating braid variants
(rotations and flips) used during the parallel sign assignment search.
"""

from __future__ import annotations

from copy import copy
from typing import (
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    TypedDict,
)

from ..domain.braid.states import BraidStates

# Type definitions
PartialSignsType = Dict[int, List[Optional[int]]]
Coordinate = Tuple[int, int]


def flip(braid: List[int]) -> List[int]:
    """
    Mirror (flip) the braid's generator indices.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.

    Returns
    -------
    list[int]
        The flipped braid word.
    """
    n: int = BraidStates(braid).n_strands
    return [(2 * (x > 0) - 1) * (n - abs(x)) for x in braid]


def rotate(braid: List[int]) -> List[int]:
    """
    Rotate the braid word left by one generator.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.

    Returns
    -------
    list[int]
        The rotated braid word with the first generator moved to the end.
    """
    return braid[1:] + [braid[0]]


def flip_coordinate(coord: Coordinate, n_strands: int) -> Coordinate:
    """Mirror a lattice coordinate (i, j) horizontally for an n-strand braid."""
    i, j = coord
    return (n_strands - 1 - i, j)


def component_locations(braid: List[int]) -> List[List[Coordinate]]:
    """
    Traverse the braid closure and return a list of components,
    each as an ordered list of lattice coordinates (i,j).
    """
    n_cross = len(braid)
    n_strands = 1 + max(abs(g) for g in braid) if braid else 0

    top_inputs = {(abs(braid[j]) - 1, j) for j in range(n_cross)}
    bottom_inputs = {(abs(braid[j]), j) for j in range(n_cross)}

    visited = set()
    comps: List[List[Coordinate]] = []

    for s in range(n_strands):
        start = (s, 0)
        if start in visited:
            continue

        loc = list(start)
        path: List[Coordinate] = []
        seen = set()

        while True:
            t = tuple(loc)
            path.append(t)
            visited.add(t)
            seen.add(t)

            if t in bottom_inputs:
                loc = [loc[0] - 1, loc[1] + 1]
            elif t in top_inputs:
                loc = [loc[0] + 1, loc[1] + 1]
            elif loc[1] == n_cross:
                loc = [loc[0], 0]
            else:
                loc = [loc[0], loc[1] + 1]
            if tuple(loc) in seen:
                break

        comps.append(path)

    return comps


def build_coord_signs(braid: List[int], strand_signs: Dict[int, List[int]]) -> Dict[Coordinate, int]:
    """Expand strand_signs into a coordinate->sign dictionary."""
    comps = component_locations(braid)
    coord_signs: Dict[Coordinate, int] = {}
    for comp_idx, comp_coords in enumerate(comps):
        signs = iter(strand_signs[comp_idx] + [strand_signs[comp_idx][0]])
        last_x = -1
        last_sign = 0
        for coord in comp_coords:
            if coord[0] != last_x:
                last_sign = next(signs)
                last_x = coord[0]
            coord_signs[coord] = last_sign
    return coord_signs


def collapse_coord_signs(coord_signs: Dict[Coordinate, int], braid: List[int]) -> Dict[int, List[int]]:
    """Collapse coordinate->signs into strand_signs (component-indexed cyclic lists)."""
    comps = component_locations(braid)
    strand_signs: Dict[int, List[int]] = {}
    for comp_idx, comp_coords in enumerate(comps):
        strand_signs[comp_idx] = []
        last_x = -1
        for c in comp_coords:
            if c[0] != last_x:
                strand_signs[comp_idx].append(coord_signs[c])
                last_x = c[0]
        strand_signs[comp_idx].pop()
    return strand_signs


def flip_partial_signs(braid: List[int], strand_signs: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """Perform a horizontal flip of an ansatz sign assignment."""
    n_strands = 1 + max(abs(g) for g in braid)
    coord_signs = build_coord_signs(braid, strand_signs)
    flipped_coord_signs = {flip_coordinate(c, n_strands): s for c, s in coord_signs.items()}
    braid_flipped = flip(braid)
    return collapse_coord_signs(flipped_coord_signs, braid_flipped)


def rotate_partial_signs(braid: List[int], strand_signs: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """Perform a rotation of an ansatz sign assignment."""
    coord_signs = build_coord_signs(braid, strand_signs)
    rotated_coord_signs = {(c[0], (c[1] - 1) % (len(braid) + 1)): s for c, s in coord_signs.items()}
    braid_rotated = braid[1:] + [braid[0]]
    return collapse_coord_signs(rotated_coord_signs, braid_rotated)


class VariantMeta(TypedDict):
    flipped: bool
    shift: int


def generate_braid_variants(
    braid_state: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
    include_flip: bool = True,
    max_shifts: Optional[int] = None,
) -> Iterator[Tuple[BraidStates, Optional[PartialSignsType], VariantMeta]]:
    """
    Enumerate braid variants by cyclic rotation (and optional flipping).

    Yields
    ------
    (BraidStates, partial_signs, meta)
        Tuples of (variant_braid_state, partial_signs, meta) where
        meta = {"flipped": bool, "shift": int}.
    """
    limit: int = len(braid_state.braid) if max_shifts is None else min(len(braid_state.braid), max_shifts)
    current: BraidStates = copy(braid_state)
    flipped: bool = False
    if partial_signs:
        for shift in range(limit):
            yield current, partial_signs, {"flipped": flipped, "shift": shift}
            if include_flip:
                partial_signs = flip_partial_signs(current.braid, partial_signs)
                current = BraidStates(flip(current.braid))
                flipped = not flipped
                yield current, partial_signs, {"flipped": flipped, "shift": shift}
            partial_signs = rotate_partial_signs(current.braid, partial_signs)
            current = BraidStates(rotate(current.braid))
    else:
        for shift in range(limit):
            yield current, None, {"flipped": flipped, "shift": shift}
            if include_flip:
                flipped = not flipped
                current = BraidStates(flip(current.braid))
                yield current, None, {"flipped": flipped, "shift": shift}
            current = BraidStates(rotate(current.braid))
