from __future__ import annotations

import multiprocessing as mp
import os
from enum import Enum
from typing import (
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    TypedDict,
    Union,
)

from copy import copy
from fkcompute.braidstates_links import BraidStates
from fkcompute.relations_links import full_reduce, check_sign_assignment
from fkcompute.braids import is_homogeneous_braid
import json


# --- classification -----------------------------------------------------------
class BraidType(Enum):
    """Classification of braids used by this module."""

    HOMOGENEOUS = 0
    FIBERED = 1


# Backwards-compatibility with previous constants
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
        `BraidType.HOMOGENEOUS.value` if homogeneous, else `BraidType.FIBERED.value`.
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
    and maps each generator index ``x`` to ``(2*(x>0)-1) * (n - |x|)``.

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

# --- utilities ----------------------------------------------------------------
def total_assignments_for_braid(braid_states: BraidStates, partial_signs: Optional[PartialSignsType] = None) -> int:
    """
    Return the total number of sign assignments for this braid.

    Parameters
    ----------
    braid_states
        Initialized :class:`BraidStates` for the target braid.

    Returns
    -------
    int
        ``2 ** braid_states.n_s_total``.
    """
    fixed_n_s_total = (
        sum(
            abs(x) for values in partial_signs.values() for x in values if x is not None
        )
        if partial_signs
        else 0
    )
    return 1 << (braid_states.n_s_total-fixed_n_s_total)


def split_ranges(n: int, chunk_size: int) -> Iterator[Tuple[int, int]]:
    """
    Yield half-open integer ranges ``[start, end)`` covering ``0..n-1``.

    Parameters
    ----------
    n
        Total count to cover.
    chunk_size
        Maximum size of each range.

    Yields
    ------
    (int, int)
        Consecutive ranges ``[start, end)``.
    """
    start = 0
    while start < n:
        end = min(start + chunk_size, n)
        yield start, end
        start = end

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

# --- multiprocessing search ---------------------------------------------------
class FailureResult(TypedDict):
    inversion_data: str  # 'failure'


class SuccessResult(TypedDict):
    link_type: str  # 'homogeneous' or 'fibered'
    inversion_data: Dict[int, List[int]]
    braid: List[int]
    degree: int


ResultDict = Union[SuccessResult, FailureResult]


def worker_on_range(
    args: Tuple[int, BraidStates, int, int, bool, Optional[PartialSignsType]],
) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
    """
    Worker for a *range* of indices; returns the first success in ``[start, end)``, else ``None``.

    Parameters
    ----------
    args
        Tuple ``(degree, braid, start, end, progress, braid_state)``.
        If ``progress`` is ``True``, prints periodic progress for this batch.

    Returns
    -------
    tuple[list[int], list[list[int]]] or None
        First successful ``(braid, strand_signs)`` found in the range, else ``None``.
    """
    degree, braid_state, start, end, progress, partial_signs = args
    n = end - start
    for i, idx in enumerate(range(start, end)):
        if progress and n > 0 and (i % max(1, n // 10) == 0):
            pct = round(100.0 * i / n, 1)
            print(f"[worker {os.getpid()}] {pct}% of batch", end="\r")
        hit = check_assignment_for_braid(idx, degree, braid_state, partial_signs)
        if hit is not None:
            return hit
    return None


def parallel_try_sign_assignments(
    degree: int,
    braid_state: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
    max_workers: Optional[int] = None,
    chunk_size: int = 1 << 14,
    include_flip: bool = True,
    max_shifts: Optional[int] = None,
    verbose: bool = False,
) -> Union[Tuple[List[int], Dict[int, List[int]]], bool]:
    """
    Parallel search for a consistent sign assignment across braid variants.

    Strategy
    --------
    1. Iterate braid variants produced by :func:`generate_braid_variants`
       (original, flipped, and rotations).
    2. For each variant, split the index space into ranges and process
       them in a pool of worker processes.
    3. Return the first ``(braid_variant, strand_signs)`` that validates;
       otherwise return ``False`` after exhausting all variants.

    Parameters
    ----------
    degree
        Degree parameter forwarded to :func:`relations_links.check_sign_assignment`.
    braid
        Base braid word as a list of signed generator indices.
    braid_state
        :class:`BraidStates` built from ``braid`` (will be reused and mutated).
        **Note:** This object is passed to workers; ensure it is safe for
        your platform’s start method and usage pattern.
    max_workers
        Size of the process pool. Defaults to ``os.cpu_count()`` (or 1 if unknown).
    chunk_size
        Number of indices per task (default ``16384``).
    include_flip
        If ``True``, also test flipped variants of each rotation.
    max_shifts
        If provided, cap the number of rotations explored.
    verbose
        If ``True``, print progress information.

    Returns
    -------
    (list[int], list[list[int]]) or bool
        ``(braid_variant, strand_signs)`` on success, else ``False``.
    """
    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 1))

    total = total_assignments_for_braid(braid_state, partial_signs)
    if verbose:
        print(f"Search space: {total:,} assignments")

    for variant_braid_state, variant_partial_signs, meta in generate_braid_variants(
        braid_state, partial_signs, include_flip=include_flip, max_shifts=max_shifts
    ):
        if verbose:
            print(f"Trying cyclic shift by {meta['shift']} on {"flipped" if meta["flipped"] else "original"}")
            print(variant_braid_state.braid)
            print(variant_partial_signs)


        tasks: List[Tuple[int, BraidStates, int, int, bool, Optional[PartialSignsType]]] = [
            (degree, variant_braid_state, start, end, verbose, variant_partial_signs)
            for (start, end) in split_ranges(total, chunk_size)
        ]

        with mp.Pool(processes=max_workers) as pool:
            async_results = [pool.apply_async(worker_on_range, (t,)) for t in tasks]

            # Gather until first success
            for r in async_results:
                out = r.get()
                if out is not None:
                    try:
                        pool.terminate()
                    finally:
                        pool.join()
                    return out

            pool.close()
            pool.join()

    return False


def get_sign_assignment_parallel(
    braid: List[int],
    partial_signs: Optional[PartialSignsType] = None,
    degree: int = 10,
    verbose: bool = False,
    **kwargs: object,
) -> ResultDict:
    """
    Parallel version mirroring the original return schema.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    degree
        Degree parameter forwarded to :func:`relations_links.check_sign_assignment`.
    verbose
        If ``True``, print progress during search.
    **kwargs
        Forwarded to :func:`parallel_try_sign_assignments`
        (e.g., ``max_workers``, ``chunk_size``, ``include_flip``, ``max_shifts``).

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
        sign_assignment: Dict[int, List[int]] = bs.strand_signs
        return {
            "link_type": "homogeneous",
            "inversion_data": sign_assignment,
            "braid": bs.braid,
            "degree": degree,
        }

    # fibered / non-homogeneous: run parallel search
    sol = parallel_try_sign_assignments(degree, bs, partial_signs, verbose=verbose, **kwargs)
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
    name = "L10n37_0m"
    print(name)
    braid_example: List[int] = [-a for a in [2, 2, 2, -1, -2, -3, -2, -3, -2, 1, 2, 2, 2, 3]]
    partial_signs: PartialSignsType = {0:[None for _ in range(6)] , 1: [1,1,None,None,-1,1,-1,1,-1,1,1,None,1,1,-1,1,-1,None,None,1]}
    st = time.time()
    assignment = get_sign_assignment_parallel(
        braid=braid_example,
        partial_signs= None,#partial_signs,
        degree=10,
        verbose=True,
        max_workers=None,  # or set an explicit int
        chunk_size=1 << 14,  # tune if you like
        include_flip=False,
        max_shifts=None,  # or a small int to limit shifts
    )
    end = time.time()
    with open(name, "w") as f:
        json.dump(assignment,f)
    print(assignment, end="\n\n")
    print(end - st)
