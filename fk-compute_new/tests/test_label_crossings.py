import pytest


from fkcompute.inversion.permutations import label_crossings
from fkcompute.inversion.permutations import perm_to_signs


def test_label_crossings_trefoil_matrix():
    assert label_crossings([1, 1, 1]) == [
        [1, 4, 5, 2],
        [5, 2, 3, 6],
        [3, 6, 7, 4],
    ]


def test_label_crossings_figure_eight_matrix():
    assert label_crossings([1, -2, 1, -2]) == [
        [1, 4, 5, 2],
        [2, 7, 8, 3],
        [5, 8, 9, 6],
        [6, 3, 4, 7],
    ]


def test_label_crossings_requires_generator_one():
    with pytest.raises(ValueError, match=r"Generator 1 not found in braid"):
        label_crossings([2, 2, 2])


def test_label_crossings_link_like_braid_uses_full_range():
    braid = [1, 1]
    lc = label_crossings(braid)
    # Hopf-link closure has 2 components, so labels run to 2n+2.
    assert max(max(row) for row in lc) == 2 * len(braid) + 2
    assert min(min(row) for row in lc) >= 1


def test_perm_to_signs_does_not_assume_contiguous_arc_labels_for_links():
    braid = [-1, 2, -1, 2, -1]
    perm = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1]
    out = perm_to_signs(perm, braid)
    assert set(out.keys()) == {0, 1}
    assert sum(len(v) for v in out.values()) == 2 * len(braid)
