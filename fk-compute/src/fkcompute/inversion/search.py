"""
Sign assignment search for FK computation.

This module searches through sign-assignment candidates (induced by
arc-label permutations) to find a valid assignment for a braid at a
given degree. The search runs in-process for a single worker and uses
a multiprocessing pool when more workers are requested.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from itertools import islice
from enum import Enum
from functools import lru_cache
from math import comb
from typing import (
    Dict,
    Iterator,
    List,
    Literal,
    Optional,
    Tuple,
    Union,
)

from ..domain.braid.states import BraidStates
from ..domain.braid.word import is_homogeneous_braid
from ..domain.constraints.reduction import full_reduce
from .permutations import iter_perms_rot_closed, perm_to_signs_metadata
from .validation import check_sign_assignment
from .variants import (
    PartialSignsType,
    generate_braid_variants,
    matches_partial_signs,
)


def _batched(it: Iterator[List[int]], batch_size: int) -> Iterator[List[List[int]]]:
    """Yield lists of items from *it* of length up to *batch_size*."""
    if batch_size <= 0:
        raise ValueError("batch_size must be > 0")
    while True:
        batch = list(islice(it, batch_size))
        if not batch:
            return
        yield batch


class _PermChecker:
    """Check permutation candidates against one braid variant.

    Holds the per-variant state (BraidStates, arc-label metadata) so that
    candidate checks reuse the expensive topology computation. Instances are
    used directly for serial search and as process-global state for workers.
    """

    def __init__(
        self,
        braid: List[int],
        degree: int,
        partial_signs: Optional[PartialSignsType],
        arc_labels: List[int],
        arc_comp: Dict[int, int],
        n_components: int,
        weight: Optional[int] = None,
    ):
        self.state = BraidStates(braid)
        self.degree = int(degree)
        self.weight = weight
        self.partial_signs = partial_signs
        self.arc_labels = arc_labels
        self.arc_comp = arc_comp
        self.n_components = int(n_components)
        self.max_arc_label = max(arc_labels) if arc_labels else 0

    def check_signs(
        self, signs: Dict[int, List[int]]
    ) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
        if self.partial_signs is not None and not matches_partial_signs(
            signs, self.partial_signs
        ):
            return None

        state = self.state
        state.strand_signs = signs
        state.compute_matrices()
        if not state.validate():
            return None

        state.generate_position_assignments()
        relations = full_reduce(state.get_state_relations())
        out = check_sign_assignment(self.degree, relations, state, weight=self.weight)
        if out is None:
            return None
        return state.braid, state.strand_signs

    def check_perm_batch(
        self, perms: List[List[int]]
    ) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
        for perm in perms:
            if self.max_arc_label > len(perm):
                continue
            signs: Dict[int, List[int]] = {c: [] for c in range(self.n_components)}
            for label in self.arc_labels:
                signs[self.arc_comp[label]].append(1 if perm[label - 1] == label else -1)
            hit = self.check_signs(signs)
            if hit is not None:
                return hit
        return None


_WORKER_CHECKER: Optional[_PermChecker] = None


def _init_perm_worker(*args) -> None:
    global _WORKER_CHECKER
    _WORKER_CHECKER = _PermChecker(*args)


def worker_on_perm_batch(perms: List[List[int]]) -> Optional[Tuple[List[int], Dict[int, List[int]]]]:
    if _WORKER_CHECKER is None:
        raise RuntimeError("Worker state not initialized")
    return _WORKER_CHECKER.check_perm_batch(perms)


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
    weight: Optional[int] = None,
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
    out = check_sign_assignment(degree, relations, braid_state, weight=weight)

    if out is not None:
        return braid_state.braid, braid_state.strand_signs
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
    weight: Optional[int] = None,
) -> Union[Tuple[List[int], Dict[int, List[int]]], Literal[False]]:
    """
    Search for a consistent sign assignment across braid variants.

    This search enumerates only sign assignments induced by permutation
    candidates from :func:`fkcompute.inversion.permutations.iter_perms_rot_closed`.
    The *chunk_size* parameter controls how many candidate permutations are
    bundled into each worker task.

    With max_workers == 1 the search runs in-process (no multiprocessing
    overhead); otherwise candidates are distributed over a process pool.

    Returns
    -------
    (list[int], dict) or bool
        (braid_variant, strand_signs) on success, else False.
    """
    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 1))

    # Candidates are enumerated lazily; cap batches to keep IPC payloads small.
    batch_size = max(1, min(int(chunk_size), 1024))

    for variant_braid_state, variant_partial_signs, meta in generate_braid_variants(
        braid_state, partial_signs, include_flip=include_flip, max_shifts=max_shifts
    ):
        if verbose:
            print(f"Trying cyclic shift by {meta['shift']} on {'flipped' if meta['flipped'] else 'original'}")
            print(variant_braid_state.braid)
            print(variant_partial_signs)

        # Ensure this braid admits permutation candidates before searching.
        try:
            arc_labels, arc_comp, n_components = perm_to_signs_metadata(variant_braid_state.braid)
        except ValueError:
            continue

        checker_args = (
            variant_braid_state.braid,
            degree,
            variant_partial_signs,
            arc_labels,
            arc_comp,
            n_components,
            weight,
        )
        perm_iter = iter_perms_rot_closed(variant_braid_state.braid)
        perm_batches = _batched(perm_iter, batch_size)

        if max_workers == 1:
            checker = _PermChecker(*checker_args)
            for batch in perm_batches:
                out = checker.check_perm_batch(batch)
                if out is not None:
                    return out
            continue

        with mp.Pool(
            processes=max_workers,
            initializer=_init_perm_worker,
            initargs=checker_args,
        ) as pool:
            for out in pool.imap_unordered(worker_on_perm_batch, perm_batches, chunksize=1):
                if out is not None:
                    pool.terminate()
                    pool.join()
                    return out

    return False
