from __future__ import annotations


from fkcompute.domain.braid.states import BraidStates
from fkcompute.inversion.api import find_sign_assignment_full
from fkcompute.inversion.permutations import iter_perms_rot_closed, perm_to_signs


def _canon(signs: dict[int, list[int]], n_components: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(s) for s in signs.get(c, ())) for c in range(n_components))


def test_find_sign_assignment_full_dedupes_multicycle_candidates(monkeypatch) -> None:
    # This braid has multicycle permutations that map to duplicate sign diagrams.
    braid = [1, -1, -2]
    bs = BraidStates(braid)

    calls: list[int] = []

    def _ok_check_sign_assignment(degree: int, relations: list, braid_states: object, weight=None):
        # Record calls to ensure we dedupe before the expensive feasibility check.
        calls.append(int(degree))
        return {"ok": True}

    monkeypatch.setattr(
        "fkcompute.inversion.api.check_sign_assignment",
        _ok_check_sign_assignment,
    )

    got = find_sign_assignment_full(braid, degree=7)

    # Compute expected unique *validated* sign diagrams from multicycle candidates.
    expected: set[tuple[tuple[int, ...], ...]] = set()
    for perm in iter_perms_rot_closed(braid):
        signs = perm_to_signs(perm, braid)
        bs.strand_signs = signs
        bs.compute_matrices()
        assert bs.validate(), "sanity: this braid's multicycle candidates should validate"
        expected.add(_canon(signs, bs.n_components))

    got_keys = {_canon(r.sign_assignment or {}, bs.n_components) for r in got}
    assert got_keys == expected

    # Dedupe should happen before calling check_sign_assignment: one call per unique sign diagram.
    assert len(calls) == len(expected)
