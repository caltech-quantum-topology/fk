"""
C++ binary execution utilities for FK computation.

This module provides functions for locating and executing the
compiled FK computation binary.
"""

import os
import pathlib
import shutil
import subprocess
from importlib import resources
from typing import Optional


def binary_path(name: str) -> str:
    """
    Find the path to a required binary, checking PATH first then package.

    Parameters
    ----------
    name
        Name of the binary (without .exe extension).

    Returns
    -------
    str
        Path to the binary.

    Raises
    ------
    RuntimeError
        If the binary cannot be found.
    """
    exe = f"{name}.exe" if os.name == "nt" else name
    # 1) Prefer PATH if user has a system install
    found = shutil.which(exe)
    if found:
        return found
    # 2) Fall back to packaged binary inside fkcompute/_bin/
    try:
        return str(resources.files("fkcompute").joinpath("_bin", exe))
    except Exception as e:
        raise RuntimeError(
            f"Could not find required helper binary '{exe}' in PATH or package."
        ) from e


def safe_unlink(path: str) -> None:
    """
    Safely remove a file, ignoring errors.

    Parameters
    ----------
    path
        Path to the file to remove.
    """
    try:
        pathlib.Path(path).unlink(missing_ok=True)
    except Exception:
        pass


def run_fk_binary(
    ilp_file: str,
    output_file: str,
    threads: int = 1,
    verbose: bool = False,
) -> None:
    """
    Run the FK computation binary.

    Parameters
    ----------
    ilp_file
        Path to the ILP input file (without extension).
    output_file
        Path for the output file.
    threads
        Number of threads to use.
    verbose
        Whether to show binary output.

    Raises
    ------
    subprocess.CalledProcessError
        If the binary execution fails.
    """
    bin_path = binary_path("fk_main")
    cmd = [bin_path, ilp_file, output_file, "--threads", str(threads)]

    if verbose:
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
