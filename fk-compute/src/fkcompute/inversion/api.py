"""
Entry point for Phase 1 sign assignment (inversion) computation.

Provides find_sign_assignment() as the single entry point for finding
a valid sign assignment for a braid at a given degree.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from ..domain.braid.states import BraidStates
from .search import (
    BraidType,
    braid_type,
    parallel_try_sign_assignments,
)
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
    degree: int,
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
        Size of index chunks for parallel workers.
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
