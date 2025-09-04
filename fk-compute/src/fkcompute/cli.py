# src/fk-compute/cli.py
from __future__ import annotations
import sys
import argparse, json, logging, os, subprocess, sys, time, shutil, tempfile, pathlib
from typing import Optional, Dict, Any, List, Union
from importlib import resources

# Updated absolute imports (assuming these files live in the package):
from fkcompute.sign_diagram_links_parallel_ansatz import get_sign_assignment_parallel
from fkcompute.fklinks import FK as fk_ilp

logger = logging.getLogger("fk_logger")

def configure_logging(verbose: bool) -> None:
    """Set logging level based on verbose flag."""
    handler = logging.StreamHandler()
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(fmt)

    # Remove old handlers so repeated runs don’t duplicate output
    logger.handlers.clear()
    logger.addHandler(handler)

    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

def _binary_path(name: str) -> str:
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

def _safe_unlink(path: str):
    try:
        pathlib.Path(path).unlink(missing_ok=True)
    except Exception:
        pass

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
def _assess_inversion(inversion: Dict[str, Any]) -> None:
    """Validate the structure of inversion data."""
    if "inversion_data" not in inversion:
        logger.error("Inversion data missing in inversion dictionary")
    if inversion.get("inversion_data") == "failure":
        logger.error("Inversion data invalid (failure)")


def _load_inversion_file(inversion_file: str) -> Dict[str, Any]:
    """Load inversion data from a JSON file and ensure integer keys."""
    with open(inversion_file, "r") as f:
        inversion = json.load(f)
    inversion["inversion_data"] = {int(k): v for k, v in inversion["inversion_data"].items()}
    return inversion


def _load_ilp_file(ilp_file: str) -> str:
    """Load ILP data from a file."""
    with open(ilp_file, "r") as f:
        return f.read()


def _parse_int_list(s: Optional[str]) -> Optional[List[int]]:
    """Parse a string into a list of ints.

    Accepts:
      - JSON-style: "[1, -2, 3]"
      - Comma-separated: "1,-2,3"
      - Space-separated: "1 -2 3"
    """
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    # Try JSON first
    try:
        val = json.loads(s)
        if isinstance(val, list):
            return [int(x) for x in val]
    except Exception:
        pass
    # Fallback: split on commas/spaces
    parts = [p for p in s.replace(",", " ").split() if p]
    return [int(p) for p in parts]




def _parse_bool(default_true: bool):
    """Return a function that adds --foo / --no-foo toggle flags for a bool."""

    def add(parser: argparse.ArgumentParser, name: str, help_text: str):
        dest = name.replace("-", "_")
        group = parser.add_mutually_exclusive_group()
        if default_true:
            group.add_argument(f"--{name}", dest=dest, action="store_true", help=f"{help_text} (default)")
            group.add_argument(f"--no-{name}", dest=dest, action="store_false", help=f"Disable {help_text}")
            parser.set_defaults(**{dest: True})
        else:
            group.add_argument(f"--{name}", dest=dest, action="store_true", help=f"Enable {help_text}")
            group.add_argument(f"--no-{name}", dest=dest, action="store_false", help=f"{help_text} (default)")
            parser.set_defaults(**{dest: False})
        return group

    return add

# -------------------------------------------------------------------------
# Main computation function
# -------------------------------------------------------------------------
def fk(
    braid: List[int],
    degree: int,
    ilp: Optional[str] = None,
    ilp_file: Optional[str] = None,
    inversion: Optional[Dict[str, Any]] = None,
    inversion_file: Optional[str] = None,
    partial_signs: Optional[List[int]] = None,
    max_workers: int = 1,
    chunk_size: int = 1 << 14,
    include_flip: bool = True,
    max_shifts: Optional[int] = None,
    verbose: bool = False,
    save_data: bool = False,
    save_dir: str = "data",
    link_name: Optional[str] = None,
) -> Dict[str, Union[int, float]]:
    """
    Compute the FK invariant of a braid via inversion data and ILP reduction.
    """
    # --- Handle naming and save directories ---
    if save_data:
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        if not link_name:
            tstr = time.strftime("%H%M%S", time.localtime(time.time()))
            link_name = f"Unknown_{tstr}"
            logger.debug(f"No name provided for link/knot, saving as: {link_name}")
        else:
            logger.debug(f"Saving as: {link_name}")
        link_name = os.path.join(save_dir, link_name)
    else:
        link_name = "temp"

    # --- Step 1: Inversion data ---
    if ilp is None and ilp_file is None:
        logger.debug("No ilp data provided, I need to calculate it")

        if inversion is not None and inversion_file is not None:
            logger.error("Both inversion data and inversion file are passed through. Pass one or the other")

        if inversion_file is not None and inversion is None:
            logger.debug("Loading inversion from file")
            inversion = _load_inversion_file(inversion_file)

        if inversion is None and inversion_file is None:
            logger.debug("Calculating inversion data")
            inversion = get_sign_assignment_parallel(
                braid,
                partial_signs=partial_signs,
                degree=degree,
                verbose=verbose,
                max_workers=max_workers,
                chunk_size=chunk_size,
                include_flip=include_flip,
                max_shifts=max_shifts,
            )
            logger.debug("Inversion data calculated")

        _assess_inversion(inversion)

    if save_data:
        with open(link_name + "_inversion.json", "w") as f:
            json.dump(inversion, f)

    # --- Step 2: ILP data ---
    if ilp is not None and ilp_file is not None:
        logger.error("Both ilp data and ilp file are passed through. Pass one or the other")

    if ilp is None and ilp_file is not None:
        logger.debug("Loading ILP from file")
        ilp = _load_ilp_file(ilp_file)

    if ilp is None and ilp_file is None:
        logger.debug("Calculating ILP data")
        ilp = fk_ilp(
            braid,
            degree=degree,
            inversion_data=inversion["inversion_data"],
            outfile=link_name + "_ilp.csv",
        )
        logger.debug("ILP data calculated")

    # --- Step 3: FK invariant computation ---
    bin_path = _binary_path("fk_segments_links")
    if verbose:
        subprocess.run([bin_path, f"{link_name}_ilp", f"{link_name}"], check=True)
    else:
        subprocess.run([bin_path, f"{link_name}_ilp", f"{link_name}"], check=True, stdout=subprocess.DEVNULL)

    with open(link_name + ".json", "r") as f:
        result = json.load(f)

    result = result["coefficient_q_powers"]

    # Clean up intermediate files unless explicitly saving
    if not save_data:
        _safe_unlink(f"{link_name}_ilp.csv")
        _safe_unlink(f"{link_name}.json")

    return result


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute FK invariant from a braid using inversion and ILP stages."
    )

    p.add_argument(
        "--braid",
        type=str,
        required=True,
        help='Braid word as a list. Examples: "[ -2, -2, -3, 1 ]" or "-2,-2,-3,1" or "-2 -2 -3 1".',
    )
    p.add_argument("--degree", type=int, required=True, help="Degree of the invariant computation.")

    p.add_argument("--ilp-file", type=str, default=None, help="Path to a precomputed ILP file.")
    # Supplying raw ILP text via CLI is uncommon; keep for parity but optional:
    p.add_argument("--ilp", type=str, default=None, help="ILP data as a raw string (advanced).")

    p.add_argument("--inversion-file", type=str, default=None, help="Path to inversion JSON file.")
    p.add_argument(
        "--partial-signs",
        type=str,
        default=None,
        help='Optional partial sign assignments. Same formats as --braid.',
    )

    p.add_argument("--max-workers", type=int, default=1, help="Parallel workers for inversion calculation.")
    p.add_argument("--chunk-size", type=int, default=(1 << 14), help="Chunk size for parallel tasks.")
    p.add_argument("--max-shifts", type=int, default=None, help="Maximum shifts considered in inversion.")

    # Toggle flags
    add_verbose = _parse_bool(default_true=False)
    add_verbose(p, "verbose", "Verbose logging")

    add_flip = _parse_bool(default_true=True)
    add_flip(p, "include-flip", "Include flip symmetry in inversion")

    add_save = _parse_bool(default_true=False)
    add_save(p, "save-data", "Save intermediate files (inversion/ILP/JSON)")

    add_print = _parse_bool(default_true=False)
    add_print(p, "print-result", "Print the FK result to stdout")

    p.add_argument("--save-dir", type=str, default="data", help="Directory to store data when --save-data is set.")
    p.add_argument("--link-name", type=str, default=None, help="Base name for output files (without extension).")

    return p


def main(argv: Optional[List[str]]=None) -> None:
    if argv is None:
        argv = sys.argv
    parser = build_parser()
    args = parser.parse_args(argv[1:])

    braid = _parse_int_list(args.braid)
    if not braid:
        parser.error("Could not parse --braid into a non-empty list of integers.")

    partial_signs = _parse_int_list(args.partial_signs)

    # sync logger level with --verbose toggle
    
    configure_logging(args.verbose)
    result = fk(
        braid=braid,
        degree=args.degree,
        ilp=args.ilp,
        ilp_file=args.ilp_file,
        inversion=None,                      # prefer file path via --inversion-file, else compute
        inversion_file=args.inversion_file,
        partial_signs=partial_signs,
        max_workers=args.max_workers,
        chunk_size=args.chunk_size,
        include_flip=args.include_flip,
        max_shifts=args.max_shifts,
        verbose=args.verbose,
        save_data=args.save_data,
        save_dir=args.save_dir,
        link_name=args.link_name,
    )

    # Pretty-print result JSON to stdout
    if args.print_result:
        print(json.dumps(result, indent=2, sort_keys=True))
