"""
Tests matching the examples in ImproveSignDiagramSearch-Davide.nb.

The notebook computes candidate sign diagrams for knots from the FK database
using ``weightPermsClosed``, which calls ``permsRotClosed`` then ``permToSigns``.
Each test case checks the exact set of sign diagrams produced from a given braid.

Reference outputs are taken directly from the notebook's evaluated Output cells.
"""

import pytest

from fkcompute.inversion.permutations import perm_to_signs, perms_rot_closed


def sign_diagrams(braid):
    """Return all candidate sign diagrams for a braid (closed-strand convention)."""
    return [perm_to_signs(p, braid) for p in perms_rot_closed(braid)]


# ---------------------------------------------------------------------------
# 3_1  (trefoil knot)
# braid: {1, 1, 1}
# Notebook output: 6 sign diagrams
# ---------------------------------------------------------------------------

EXPECTED_3_1 = {
    (1, 1, 1, 1, 1, 1),
    (1, 1, -1, -1, -1, 1),
    (-1, -1, -1, -1, -1, -1),
    (-1, -1, -1, 1, 1, 1),
    (-1, 1, 1, 1, -1, -1),
    (-1, 1, -1, 1, -1, 1),
}


def test_3_1_count():
    assert len(sign_diagrams([1, 1, 1])) == 6


def test_3_1_diagrams():
    actual = {tuple(d[0]) for d in sign_diagrams([1, 1, 1])}
    assert actual == EXPECTED_3_1


# ---------------------------------------------------------------------------
# 4_1  (figure-eight knot)
# braid: {1, -2, 1, -2}
# Notebook output: 9 sign diagrams
# ---------------------------------------------------------------------------

EXPECTED_4_1 = {
    (1, 1, 1, 1, 1, 1, 1, 1),
    (1, 1, -1, -1, -1, -1, -1, 1),
    (1, 1, -1, 1, 1, 1, -1, 1),
    (-1, -1, -1, -1, -1, -1, -1, -1),
    (-1, -1, -1, -1, -1, 1, 1, 1),
    (-1, -1, -1, 1, 1, 1, -1, -1),
    (-1, 1, 1, 1, -1, -1, -1, -1),
    (-1, 1, 1, 1, -1, 1, 1, 1),
    (-1, 1, -1, 1, -1, 1, -1, 1),
}


def test_4_1_count():
    assert len(sign_diagrams([1, -2, 1, -2])) == 9


def test_4_1_diagrams():
    actual = {tuple(d[0]) for d in sign_diagrams([1, -2, 1, -2])}
    assert actual == EXPECTED_4_1


# ---------------------------------------------------------------------------
# 8_20
# braid: {1, -2, -1, -1, 2, 2, -1, -2}
# Notebook output: 26 sign diagrams  (16 signs each, n=8)
# ---------------------------------------------------------------------------

EXPECTED_8_20 = {
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    (1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1),
    (1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1),
    (1, 1, 1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1),
    (1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1),
    (1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, 1),
    (1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1),
    (1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1),
    (1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1),
    (1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1),
    (1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1),
    (1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1),
    (1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1),
    (1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1),
    (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (-1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1),
    (-1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1),
    (-1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1),
    (-1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1),
    (-1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1),
    (-1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (-1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1),
    (-1, 1, 1, -1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1),
    (-1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1),
    (-1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1, -1),
    (-1, 1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1),
}


def test_8_20_count():
    assert len(sign_diagrams([1, -2, -1, -1, 2, 2, -1, -2])) == 26


def test_8_20_diagrams():
    actual = {tuple(d[0]) for d in sign_diagrams([1, -2, -1, -1, 2, 2, -1, -2])}
    assert actual == EXPECTED_8_20


# ---------------------------------------------------------------------------
# 10_132
# braid: {1, 1, -3, -3, -2, -1, 2, -1, -2, -3, 2}
# Notebook output: 54 sign diagrams  (22 signs each, n=11)
# ---------------------------------------------------------------------------


def test_10_132_count():
    braid = [1, 1, -3, -3, -2, -1, 2, -1, -2, -3, 2]
    diagrams = sign_diagrams(braid)
    assert len(diagrams) == 54


def test_10_132_all_ones_present():
    """The all-+1 and all--1 diagrams must always appear."""
    braid = [1, 1, -3, -3, -2, -1, 2, -1, -2, -3, 2]
    diagrams = sign_diagrams(braid)
    n = len(braid)
    assert {0: [1] * (2 * n)} in diagrams
    assert {0: [-1] * (2 * n)} in diagrams
