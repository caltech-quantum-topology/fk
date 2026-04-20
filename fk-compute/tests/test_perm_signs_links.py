from __future__ import annotations

from collections import defaultdict
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


def test_perm_signs_links_contains_yaml_inversion() -> None:
    entries = _load_entries(Path(__file__).parent / "all_links.yaml")

    # Many link entries share the same braid but have different inversions.
    expected_by_braid: dict[tuple[int, ...], dict[tuple[tuple[int, ...], ...], list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for name, cfg in entries.items():
        braid_t = tuple(int(x) for x in cfg["braid"])
        expected = _canonicalize_signs(cfg["inversion"])
        expected_by_braid[braid_t][expected].append(name)

    for braid_t, expected_map in expected_by_braid.items():
        braid = list(braid_t)
        expected_set = set(expected_map.keys())
        found: set[tuple[tuple[int, ...], ...]] = set()

        for perm in perms_rot_closed(braid):
            canon = _canonicalize_signs(perm_to_signs(perm, braid))
            if canon in expected_set:
                found.add(canon)
                if found == expected_set:
                    break

        missing = expected_set - found
        if missing:
            names = sorted({n for exp in missing for n in expected_map[exp]})
            sample = names[:10]
            extra = "" if len(names) <= 10 else f" (+{len(names) - 10} more)"
            raise AssertionError(
                "Inversion(s) not produced by perms_rot_closed/perm_to_signs for braid "
                f"{list(braid_t)}; failing entries: {sample}{extra}"
            )
