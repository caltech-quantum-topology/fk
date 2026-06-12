from __future__ import annotations

import os
import sys
from pathlib import Path


# Ensure tests exercise the in-repo src/ layout, not a site-packages install.
_SRC = (Path(__file__).resolve().parents[1] / "src").as_posix()
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


def pytest_addoption(parser):
    parser.addoption(
        "--bar-natan-links-dir",
        action="store",
        default=None,
        help="Run Bar-Natan link checks over FK JSON files in this directory.",
    )
    parser.addoption(
        "--bar-natan-links-limit",
        action="store",
        type=int,
        default=None,
        help="Limit the number of link JSON files checked.",
    )


_LINK_JSON_TEST_FUNCTIONS = {"test_bar_natan_links_folder", "test_bar_natan_p0_p1"}


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ not in _LINK_JSON_TEST_FUNCTIONS:
        return
    if "json_file" not in metafunc.fixturenames:
        return

    directory = metafunc.config.getoption("--bar-natan-links-dir")
    if directory is None:
        directory = os.environ.get("BAR_NATAN_LINKS_DIR")
    if not directory:
        metafunc.parametrize("json_file", [])
        return

    files = [
        path
        for path in sorted(Path(directory).glob("*.json"))
        if not path.name.endswith("_inversion.json")
    ]
    limit = metafunc.config.getoption("--bar-natan-links-limit")
    if limit is None and os.environ.get("BAR_NATAN_LINKS_LIMIT"):
        limit = int(os.environ["BAR_NATAN_LINKS_LIMIT"])
    if limit is not None:
        files = files[:limit]

    metafunc.parametrize("json_file", files, ids=[path.stem for path in files])
