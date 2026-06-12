import json

import sympy as sp

from bar_natan import T, bar_natan_Z
from bar_natan_links import bar_natan_link_Z, bar_natan_link_Z_from_braid
from test_bar_natan import braid_to_pd, normalize_pd


def test_link_z_specializes_to_knot_z_for_one_component_pd():
    pd = normalize_pd(braid_to_pd([1, 1, 1]))

    assert bar_natan_link_Z(pd, Tnames=[T]) == bar_natan_Z(pd, Tname=T)
    assert bar_natan_link_Z(pd) == bar_natan_Z(pd, Tname=T)


def test_link_z_from_braid_smoke_for_hopf_link():
    x, y = sp.symbols("x y")

    p0, p1 = bar_natan_link_Z_from_braid([1, 1, 1, 1], Tnames=[x, y])

    assert sp.simplify(p0 - 1 / (x * y)) == 0
    assert sp.simplify(p1 + 2 / (x**2 * y**2)) == 0


def test_bar_natan_links_folder(json_file):
    with json_file.open() as f:
        data = json.load(f)

    metadata = data["metadata"]
    components = int(metadata["components"])
    variables = sp.symbols(f"T0:{components}")

    p0, p1 = bar_natan_link_Z_from_braid(metadata["braid"], Tnames=variables)

    assert p0 != 0, f"{json_file.name}: P0 vanished"
    assert p1 is not None, f"{json_file.name}: P1 was not computed"

