"""
Entry point for Phase 1 sign assignment (inversion) computation.

Provides find_sign_assignment() as the single entry point for finding
a valid sign assignment for a braid at a given degree.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, List, Optional

from ..domain.braid.states import BraidStates
from ..domain.constraints.reduction import full_reduce
from .permutations import compile_perm_to_signs, iter_perms_rot_closed
from .search import (
    BraidType,
    braid_type,
    parallel_try_sign_assignments,
)
from .validation import check_sign_assignment
from .variants import PartialSignsType


@dataclass
class InversionResult:
    """Result of the sign assignment search.

    Attributes
    ----------
    success
        Whether a valid sign assignment was found.
    braid
        The (possibly transformed) braid word.
    sign_assignment
        The sign assignment dict, or None if not found.
    braid_type
        "homogeneous" or "fibered".
    degree
        The degree used for the search.
    """
    success: bool
    braid: List[int]
    sign_assignment: Optional[Dict[int, List[int]]]
    braid_type: str
    degree: int

    def to_result_dict(self) -> dict:
        """Convert to the legacy result dictionary format."""
        if not self.success:
            return {"inversion_data": "failure"}
        return {
            "link_type": self.braid_type,
            "inversion_data": self.sign_assignment,
            "braid": self.braid,
            "degree": self.degree,
        }


def find_sign_assignment(
    braid: List[int],
    degree: int = 10,
    partial_signs: Optional[PartialSignsType] = None,
    max_workers: int = 1,
    chunk_size: int = 1 << 14,
    include_flip: bool = True,
    max_shifts: Optional[int] = None,
    verbose: bool = False,
) -> InversionResult:
    """Find a valid sign assignment for the given braid and degree.

    For homogeneous braids, the sign assignment is determined automatically.
    For fibered braids, searches through braid variants and sign combinations.

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    degree
        Degree bound for the computation.
    partial_signs
        Optional partial sign constraints.
    max_workers
        Number of parallel workers for the search.
    chunk_size
        Number of permutation candidates per worker task (internally capped for IPC safety).
    include_flip
        Whether to include flipped braid variants.
    max_shifts
        Maximum number of cyclic shifts to try.
    verbose
        If True, print progress during search.

    Returns
    -------
    InversionResult
        Result with success status, braid, sign assignment, and metadata.
    """
    t = braid_type(braid)
    bs = BraidStates(braid)

    if t == BraidType.HOMOGENEOUS.value:
        return InversionResult(
            success=True,
            braid=bs.braid,
            sign_assignment=bs.strand_signs,
            braid_type="homogeneous",
            degree=degree,
        )

    sol = parallel_try_sign_assignments(
        degree, bs, partial_signs,
        max_workers=max_workers,
        chunk_size=chunk_size,
        include_flip=include_flip,
        max_shifts=max_shifts,
        verbose=verbose,
    )
    if sol is False:
        return InversionResult(
            success=False,
            braid=braid,
            sign_assignment=None,
            braid_type="fibered",
            degree=degree,
        )

    new_braid, sign_assignment = sol
    return InversionResult(
        success=True,
        braid=new_braid,
        sign_assignment=sign_assignment,
        braid_type="fibered",
        degree=degree,
    )


def find_sign_assignment_full(
    braid: List[int],
    degree: int = 10,
    partial_signs: Optional[PartialSignsType] = None,
    *,
    verbose: bool = False,
) -> List[InversionResult]:
    """Find all *unique* valid sign assignments induced by multicycle candidates.

    This enumerates only the sign diagrams coming from permutation candidates
    produced by :func:`fkcompute.inversion.permutations.iter_perms_rot_closed`.

    Differences vs :func:`find_sign_assignment`:
    - no braid variants are tried (no flips, no cyclic shifts)
    - does not stop at the first hit; returns all unique hits

    Validity predicate:
    - the sign diagram must pass :meth:`fkcompute.domain.braid.signed.SignedBraid.validate`
    - and :func:`fkcompute.inversion.validation.check_sign_assignment` for the given degree

    Parameters
    ----------
    braid
        Braid word as a list of signed generator indices.
    degree
        Degree bound for the feasibility check.
    partial_signs
        Optional per-component sign constraints (None entries mean free).
    verbose
        If True, print coarse progress information.

    Returns
    -------
    list[InversionResult]
        One result per unique valid sign assignment. Returns an empty list if
        no valid assignments exist.
    """

    def _matches_partial(
        signs: Dict[int, List[int]],
        partial: PartialSignsType,
    ) -> bool:
        for comp, templ in partial.items():
            cand = signs.get(comp)
            if cand is None or len(cand) != len(templ):
                return False
            for s_cand, s_req in zip(cand, templ):
                if s_req is not None and int(s_cand) != int(s_req):
                    return False
        return True

    bs = BraidStates(braid)
    t = braid_type(braid)
    braid_kind = "homogeneous" if t == BraidType.HOMOGENEOUS.value else "fibered"

    perm_to_signs = compile_perm_to_signs(braid)

    # Deduplicate by the actual sign diagram; multiple permutations can induce
    # the same strand_signs.
    seen: set[tuple[tuple[int, ...], ...]] = set()
    out: List[InversionResult] = []

    if verbose:
        print("Enumerating multicycle candidates (closed permutations)")

    for perm in iter_perms_rot_closed(braid):
        try:
            signs = perm_to_signs(perm)
        except ValueError:
            # Defensive: skip malformed candidates.
            continue

        if partial_signs is not None and not _matches_partial(signs, partial_signs):
            continue

        key = tuple(
            tuple(int(s) for s in signs.get(c, ())) for c in range(bs.n_components)
        )
        if key in seen:
            continue
        seen.add(key)

        bs.strand_signs = signs
        bs.compute_matrices()
        if not bs.validate():
            continue

        bs.generate_position_assignments()
        relations = full_reduce(bs.get_state_relations())
        if check_sign_assignment(degree, relations, bs) is None:
            continue

        out.append(
            InversionResult(
                success=True,
                braid=list(braid),
                sign_assignment=deepcopy(signs),
                braid_type=braid_kind,
                degree=degree,
            )
        )

    return out
