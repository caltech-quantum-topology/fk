"""
Parallel sign assignment search for FK computation.

This module provides functions for searching through sign assignment
configurations in parallel to find valid assignments.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from enum import Enum
from functools import lru_cache
from math import comb
from typing import (
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    Union,
)

from ..domain.braid.states import BraidStates
from ..domain.braid.word import is_homogeneous_braid
from ..domain.constraints.reduction import full_reduce
from .validation import check_sign_assignment
from .variants import (
    PartialSignsType,
    generate_braid_variants,
)


class BraidType(Enum):
    """Classification of braids used by this module."""
    HOMOGENEOUS = 0
    FIBERED = 1


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
        BraidType.HOMOGENEOUS.value if homogeneous, else BraidType.FIBERED.value.
    """
    return (
        BraidType.HOMOGENEOUS.value
        if is_homogeneous_braid(braid)
        else BraidType.FIBERED.value
    )


def total_assignments_for_braid(braid_states: BraidStates, partial_signs: Optional[PartialSignsType] = None) -> int:
    """Return the total number of sign assignments for this braid."""
    fixed_n_s_total = (
        sum(
            abs(x) for values in partial_signs.values() for x in values if x is not None
        )
        if partial_signs
        else 0
    )
    return 1 << (braid_states.n_s_total  - fixed_n_s_total)


def split_ranges(n: int, chunk_size: int) -> Iterator[Tuple[int, int]]:
    """Yield half-open integer ranges [start, end) covering 0..n-1."""
    start = 0
    while start < n:
        end = min(start + chunk_size, n)
        yield start, end
        start = end


@lru_cache(maxsize=None)
def _n_choose_k(n: int, k: int) -> int:
    """Cached binomial coefficient used by rank decoding."""
    if k < 0 or k > n:
        return 0
    return comb(n, k)


def _free_signs_from_rank(rank: int, n_free: int) -> List[int]:
    """Decode rank into a +/-1 list with the desired search ordering.

    Ordering (for n_free = 4):
        ++++, +++-, ++-+, +-++, -+++, ++--, ...

    This is: sort by number of -1 entries (fewest first), and within a fixed
    number of -1 entries, sort by the underlying bitstring (where 1 means -1)
    in increasing lexicographic / numeric order.
    """
    if n_free < 0:
        raise ValueError("n_free must be >= 0")
    if n_free == 0:
        if rank != 0:
            raise ValueError("rank out of range for n_free=0")
        return []

    total = 1 << n_free
    if rank < 0 or rank >= total:
        raise ValueError("rank out of range")

    # First decide k = number of '-'
    k = 0
    r = rank
    while k <= n_free:
        block = _n_choose_k(n_free, k)
        if r < block:
            break
        r -= block
        k += 1
    if k > n_free:
        raise ValueError("rank out of range")

    # Unrank within the k-block in lex order over bitstrings (MSB first).
    bits: List[int] = []
    ones_left = k
    r_within = r
    for i in range(n_free):
        if ones_left == 0:
            bits.extend([0] * (n_free - i))
            break
        if n_free - i == ones_left:
            bits.extend([1] * ones_left)
            break

        # If we put a 0 here, count how many completions exist.
        count_if_zero = _n_choose_k(n_free - i - 1, ones_left)
        if r_within < count_if_zero:
            bits.append(0)
        else:
            bits.append(1)
            r_within -= count_if_zero
            ones_left -= 1

    # Map bit 0 -> +1, bit 1 -> -1.
    return [1 if b == 0 else -1 for b in bits]


def sign_assignment_from_index(
    index: int,
    braid_states: BraidStates,
    partial_signs: Optional[PartialSignsType] = None,
) -> Dict[int, List[int]]:
    """
    Decode an integer into per-component +/-1 sign lists.

    The binary expansion of index is interpreted bitwise via sign = 2*bit - 1.
    """
    fixed_n_s_total = (
        sum(
            abs(x) for values in partial_signs.values() for x in values if x is not None
        )
        if partial_signs
        else 0
    )
    n_free = braid_states.n_s_total - fixed_n_s_total
    signs: List[int] = _free_signs_from_rank(index, n_free)
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
            out[component] = signs[cursor: cursor + k]
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

    Returns
    -------
    tuple[list[int], dict] or None
        (braid, strand_signs) on success, otherwise None.
    """
    braid_state.strand_signs = sign_assignment_from_index(
        index, braid_state, partial_signs=partial_signs
    )

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


def worker_on_range(
    args: Tuple[int, BraidStates, int, int, bool, Optional[PartialSignsType]],
) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
    """
    Worker for a range of indices; returns the first success in [start, end).
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

    Returns
    -------
    (list[int], dict) or bool
        (braid_variant, strand_signs) on success, else False.
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
            print(f"Trying cyclic shift by {meta['shift']} on {'flipped' if meta['flipped'] else 'original'}")
            print(variant_braid_state.braid)
            print(variant_partial_signs)

        tasks: List[Tuple[int, BraidStates, int, int, bool, Optional[PartialSignsType]]] = [
            (degree, variant_braid_state, start, end, verbose, variant_partial_signs)
            for (start, end) in split_ranges(total, chunk_size)
        ]

        with mp.Pool(processes=max_workers) as pool:
            async_results = [pool.apply_async(worker_on_range, (t,)) for t in tasks]

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
