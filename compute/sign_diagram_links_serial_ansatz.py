from __future__ import annotations

from copy import copy
from enum import Enum
from typing import (
    Iterable,
    Iterator,
    List,
    Dict,
    Optional,
    Tuple,
    TypedDict,
    Union,
)

from braidstates_links import BraidStates
from relations_links import full_reduce, check_sign_assignment
from braids import is_homogeneous_braid
import json


# --- classification constants -------------------------------------------------
class BraidType(Enum):
    """Classification of braids used by this module."""

    HOMOGENEOUS = 0
    FIBERED = 1


# Backwards-compatibility with existing code that may import these names
BRAID_TYPE_HOMEGENOUS = BraidType.HOMOGENEOUS.value
BRAID_TYPE_FIBERED = BraidType.FIBERED.value


def braid_type(braid: List[int]) -> int:
    """
    Classify a braid as homogeneous or fibered.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.

    Returns
    -------
    int
        `BraidType.HOMOGENEOUS.value` if the braid is homogeneous,
        else `BraidType.FIBERED.value`.
    """
    return (
        BraidType.HOMOGENEOUS.value
        if is_homogeneous_braid(braid)
        else BraidType.FIBERED.value
    )


# --- simple transforms --------------------------------------------------------
def flip(braid: List[int]) -> List[int]:
    """
    Mirror (flip) the braid’s generator indices.

    The flip uses the number of strands determined by :class:`BraidStates`
    and maps each generator index ``x`` to ``(2*(x>0)-1) * (n - |x|)``,
    effectively reflecting across the middle strand while preserving sign.

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


# --- partial_signs transforms --------------------------------------------------

PartialSignsType = Dict[int, List[Optional[int]]]
Coordinate = Tuple[int, int]

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

    # Precompute top/bottom input locations
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
            elif loc[1] == n_cross:        # wrap around closure
                loc = [loc[0], 0]
            else:                          # go straight down
                loc = [loc[0], loc[1] + 1]
            if tuple(loc) in seen:
                break

        comps.append(path)

    return comps

def build_coord_signs(braid: List[int], strand_signs: Dict[int, List[int]]) -> Dict[Coordinate, int]:
    """
    Expand strand_signs into a coordinate→sign dictionary
    by following the component geometry of the braid closure.
    """
    comps = component_locations(braid)
    coord_signs: Dict[Coordinate, int] = {}
    for comp_idx, comp_coords in enumerate(comps):
        signs = iter(strand_signs[comp_idx] + [strand_signs[comp_idx][0]])
        L = len(comp_coords)
        last_x = -1
        last_sign = 0
        for coord in comp_coords:
            if coord[0] != last_x:
                last_sign = next(signs)
                last_x = coord[0]
            coord_signs[coord] = last_sign
    return coord_signs

def collapse_coord_signs(coord_signs: Dict[Coordinate, int],
                         braid: List[int]) -> Dict[int, List[int]]:
    """
    Collapse coordinate→signs into strand_signs (component-indexed cyclic lists).
    """
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

def flip_partial_signs(braid: List[int],
                            strand_signs: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """
    Perform a horizontal flip of an ansatz sign assignment.
    """
    n_strands = 1 + max(abs(g) for g in braid)
    # Step 1. Expand to coordinate-level dict
    coord_signs = build_coord_signs(braid, strand_signs)
    # Step 2. Flip coordinates
    flipped_coord_signs = {flip_coordinate(c, n_strands): s for c, s in coord_signs.items()}
    # Step 3. Flip braid word
    braid_flipped = flip(braid)
    # Step 4. Collapse back into strand_signs for flipped braid
    return collapse_coord_signs(flipped_coord_signs, braid_flipped)


def rotate_partial_signs(braid: List[int],
                            strand_signs: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """
    Perform a rotation of an ansatz sign assignment.
    """
    n_strands = 1 + max(abs(g) for g in braid)

    # Step 1. Expand to coordinate-level dict
    coord_signs = build_coord_signs(braid, strand_signs)
    # Step 2. Rotate coordinates
    rotated_coord_signs = {(c[0], (c[1]-1) % (len(braid)+1)): s for c, s in coord_signs.items()}
    # Step 3. Rotate braid word
    braid_rotated = braid[1:] + [braid[0]]

    # Step 4. Collapse back into strand_signs for rotated braid
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

    For ``shift = 0..limit-1`` (where ``limit`` is ``len(braid)`` or
    ``min(len(braid), max_shifts)``), yields the rotated braid and,
    if requested, its flipped counterpart.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    include_flip
        If ``True``, also yield the flipped braid for each rotation.
    max_shifts
        If provided, cap the number of rotations explored.

    Yields
    ------
    (list[int], dict)
        Pairs of ``(variant_braid, meta)`` where ``meta = {"flipped": bool,
        "shift": int}``.
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


# --- sign assignment utilities -----------------------------------------------
def sign_assignment_from_index(
    index: int,
    braid_states: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
) -> Dict[int, List[int]]:
    """
    Decode an integer into per-component ±1 sign lists.

    The binary expansion of ``index`` over ``braid_states.n_s_total`` bits
    is interpreted bitwise via ``sign = 2*bit - 1``, i.e., ``0 -> -1`` and
    ``1 -> +1``. The flat sign list is then split across components according
    to ``braid_states.n_s``.

    Parameters
    ----------
    index
        Integer in ``[0, 2**n_s_total - 1]`` selecting a sign pattern.
    braid_states
        Initialized :class:`BraidStates` for the target braid.

    Returns
    -------
    list[list[int]]
        A list of components; each component is a list of signs (±1).
         Shape: ``[n_components][n_s[component]]``.
    """

    fixed_n_s_total = (
        sum(
            abs(x) for values in partial_signs.values() for x in values if x is not None
        )
        if partial_signs
        else 0
    )
    bits: List[int] = list(
        map(int, list(bin(index)[2:].zfill(braid_states.n_s_total - fixed_n_s_total)))
    )
    signs: List[int] = [2 * b - 1 for b in bits]  # 0 -> -1, 1 -> +1
    out: Dict[int, List[int]] = {}
    if partial_signs:
        iter_signs = iter(signs)
        for component in range(braid_states.n_components):
            out[component] = [
                s if s else next(iter_signs) for s in partial_signs[component]
            ]
    else:
        cursor = 0
        for component in range(braid_states.n_components):
            k = braid_states.n_s[component]
            out[component] = signs[cursor : cursor + k]
            cursor += k
    return out


def check_assignment_for_braid(
    index: int,
    degree: int,
    braid_state: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
    """
    Test a single sign assignment (by index) for a given braid.

    The function mutates ``braid_state.strand_signs`` to match the sign
    assignment decoded from ``index``, runs the necessary validation and
    relation checks, and returns the successful braid/assignment pair if
    consistent.

    Parameters
    ----------
    index
        Integer in ``[0, 2**n_s_total - 1]`` selecting a sign pattern.
    degree
        Degree parameter forwarded to :func:`relations_links.check_sign_assignment`.
    braid_state
        :class:`BraidStates` instance for ``braid`` (will be mutated).

    Returns
    -------
    tuple[list[int], list[list[int]]] or None
        ``(braid, strand_signs)`` on success, otherwise ``None``.
    """
    # fill strand signs from index
    braid_state.strand_signs = sign_assignment_from_index(
        index, braid_state, partial_signs=partial_signs
    )

    # compute & validate
    braid_state.compute_matrices()
    if not braid_state.validate():
        return None

    braid_state.generate_position_assignments()
    all_relations = braid_state.get_state_relations()
    relations = full_reduce(all_relations)
    out = check_sign_assignment(degree, relations, braid_state)

    if out is not None:
        return braid_state.braid, braid_state.strand_signs
    return None


# --- serial search ------------------------------------------------------------
class FailureResult(TypedDict):
    inversion_data: str  # 'failure'


class SuccessResult(TypedDict):
    link_type: str  # 'homogeneous' or 'fibered'
    inversion_data: Dict[int, List[int]]
    braid: List[int]
    degree: int


ResultDict = Union[SuccessResult, FailureResult]


def try_sign_assignments_serial(
    degree: int,
    braid: List[int],
    braid_state: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
    include_flip: bool = True,
    max_shifts: Optional[int] = None,
    verbose: bool = False,
) -> Union[Tuple[List[int], Dict[int, List[int]]], bool]:
    """
    Search serially for a consistent sign assignment across braid variants.

    The search iterates braid variants produced by
    :func:`generate_braid_variants` in the order:
    original (and flipped), then 1-step rotation (and flipped), etc.,
    up to ``len(braid)`` rotations or ``max_shifts`` if provided. For each
    variant, it enumerates all ``2**n_s_total`` sign patterns until one
    validates.

    Parameters
    ----------
    degree
        Degree parameter forwarded to :func:`relations_links.check_sign_assignment`.
    braid
        Base braid word as a list of signed generator indices.
    braid_state
        :class:`BraidStates` built from ``braid`` (will be reused and mutated).
    partial_signs
        Optional partial sign assignment
    include_flip
        If ``True``, also test flipped variants of each rotation.
    max_shifts
        If provided, cap the number of rotations explored.
    verbose
        If ``True``, print progress information.

    Returns
    -------
    (list[int], list[list[int]]) or bool
        ``(braid_variant, strand_signs)`` on success. ``False`` if no
        assignment is found across all variants.
    """
    variants_total = 2 * len(braid) if include_flip else len(braid)
    variant_idx = -1
    fixed_n_s_total = sum(
        abs(x) for values in partial_signs.values() for x in values if x is not None
    ) if partial_signs else 0
    total = 2 ** (braid_state.n_s_total - fixed_n_s_total)
    if verbose:
        print(f"  Search space: {total:,} assignments")

    for variant_braid_state, partial_sign_variant, meta in generate_braid_variants(
        braid_state,
        partial_signs=partial_signs,
        include_flip=include_flip,
        max_shifts=max_shifts,
    ):
        variant_idx += 1
        if verbose:
            head = "flipped" if meta["flipped"] else "original"
            print(f"Trying rotation shift={meta['shift']} ({head})...")
            print(variant_braid_state.braid)
            print(partial_sign_variant)

        # serial enumeration
        for idx in range(total):
            if verbose:
                pct = round(
                    100.0
                    * (variant_idx / variants_total + idx / (variants_total * total)),
                    2,
                )
                print(f"    {pct}%", end="\r")
            hit = check_assignment_for_braid(
                idx,
                degree,
                variant_braid_state,
                partial_signs=partial_sign_variant,
            )
            if hit is not None:
                if verbose:
                    print("\n  ✓ Found valid assignment.")
                return hit

        if verbose:
            print("  No valid assignment for this variant.\n")

    # No success across all variants
    return False


def get_sign_assignment_serial(
    braid: List[int],
    degree: int = 10,
    verbose: bool = False,
    partial_signs: Optional[PartialSignsType] = None,
    **kwargs: object,
) -> ResultDict:
    """
    Compute a sign assignment for a braid without parallelism.

    This mirrors the original ``get_sign_assignment`` schema:
    homogeneous braids immediately return the initial strand signs from
    :class:`BraidStates`; otherwise, a serial search is performed across
    rotations (and optional flips) until a consistent assignment is found.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    degree
        Degree parameter forwarded to :func:`relations_links.check_sign_assignment`.
    verbose
        If ``True``, print progress during search.
    ansatz
        Ansatz for
    **kwargs
        Accepted but ignored (for drop-in replacement flexibility).

    Returns
    -------
    dict
        One of:
            * ``{'link_type': 'homogeneous', 'inversion_data': signs, 'braid': braid, 'degree': degree}``
            * ``{'link_type': 'fibered', 'inversion_data': signs, 'braid': braid_variant, 'degree': degree}``
            * ``{'inversion_data': 'failure'}`` if no assignment is found.
    """
    t = braid_type(braid)
    bs = BraidStates(braid)

    if t == BraidType.HOMOGENEOUS.value:
        # homogeneous: just return the initial strand_signs
        sign_assignment: Dict[int, List[int]] = bs.strand_signs
        return {
            "link_type": "homogeneous",
            "inversion_data": sign_assignment,
            "braid": bs.braid,
            "degree": degree,
        }

    # fibered / non-homogeneous: run serial search
    sol = try_sign_assignments_serial(
        degree, braid, bs, partial_signs=partial_signs, verbose=verbose
    )
    if sol is False:
        return {"inversion_data": "failure"}

    new_braid, sign_assignment = sol
    return {
        "link_type": "fibered",
        "inversion_data": sign_assignment,
        "braid": new_braid,
        "degree": degree,
    }


if __name__ == "__main__":
    import time

    braid_example: List[int] = [-1, 2, 3, -2, -1, -2, -1, -2, -3, 2]
    partial_signs: PartialSignsType = {0: [-1, -1, 1, None, None, None, None, None, -1, -1, 1, -1, -1, -1], 1: [-1, -1, None, None, None, -1]}
    st = time.time()
    assignment = get_sign_assignment_serial(
        braid=braid_example,
        partial_signs= partial_signs,
        degree=10,
        verbose=True,
        max_workers=None,
        chunk_size=1 << 14,
        include_flip=True,
        max_shifts=None,
    )
    end = time.time()
    print(assignment, end="\n\n")
    print(end - st)
