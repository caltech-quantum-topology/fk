from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import yaml

from fkcompute.inversion.permutations import perm_to_signs, perms_rot_closed


def _canonicalize_signs(signs: Mapping[int, Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    if not signs:
        return ()
    n_components = max(int(k) for k in signs.keys()) + 1
    return tuple(tuple(int(s) for s in signs.get(c, ())) for c in range(n_components))


def _load_entries(path: Path) -> dict[str, dict]:
    with path.open() as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise TypeError(f"Expected mapping in {path}, got {type(data).__name__}")
    return data


def test_perm_signs_knots_contains_yaml_inversion() -> None:
    entries = _load_entries(Path(__file__).parent / "all_knots.yaml")

    for name, cfg in entries.items():
        braid = list(cfg["braid"])
        expected = _canonicalize_signs(cfg["inversion"])

        found = False
        for perm in perms_rot_closed(braid):
            if _canonicalize_signs(perm_to_signs(perm, braid)) == expected:
                found = True
                break

        assert found, f"{name}: inversion not produced by perms_rot_closed/perm_to_signs"
