"""
Core FK computation function.

This module provides the main fk() function for computing FK invariants.
"""

import json
import logging
import os
import tempfile
import time
from typing import Any, Dict, List, Optional, Union

from .batch import fk_from_config
from ..domain.braid.states import BraidStates
from ..inversion.api import find_sign_assignment
from ..infra.binary import run_fk_binary
from ..infra.config import load_inversion_file, load_ilp_file
from ..output.symbolic import format_symbolic_output, SYMPY_AVAILABLE

# Logger setup
logger = logging.getLogger("fk_logger")


class SignAssignmentError(RuntimeError):
    """Raised when no valid sign assignment exists for a braid at a degree."""


def configure_logging(verbose: bool) -> None:
    """Set logging level based on verbose flag."""
    handler = logging.StreamHandler()
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(fmt)

    logger.handlers.clear()
    logger.addHandler(handler)

    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)


def _assess_inversion(inversion: Dict[str, Any], braid: List[int], degree: int) -> None:
    """Validate the structure of inversion data, raising on failure."""
    if inversion.get("inversion_data") in (None, "failure"):
        raise SignAssignmentError(
            f"No valid sign assignment found for braid {braid} at degree {degree}. "
            "Try a higher degree, a different braid representative, or pass "
            "known inversion data via the 'inversion' parameter."
        )


# Defaults for every parameter accepted by _fk_compute beyond braid/degree.
_DEFAULT_PARAMS: Dict[str, Any] = {
    'verbose': False,
    'max_workers': 1,
    'chunk_size': 1 << 14,
    'include_flip': False,
    'max_shifts': None,
    'save_data': False,
    'save_dir': "data",
    'link_name': None,
    'symbolic': False,
    'threads': 1,
    'ilp': None,
    'ilp_file': None,
    'inversion': None,
    'inversion_file': None,
    'partial_signs': None,
    'weight': None,
    '_progress_callback': None,
}


def fk(
    braid_or_config: Union[List[int], str],
    degree: Optional[int] = None,
    config: Optional[str] = None,
    symbolic: bool = False,
    threads: Optional[int] = None,
    name: Optional[str] = None,
    weight: Optional[int] = None,
    **kwargs,
) -> Union[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """
    Unified FK invariant computation with intelligent interface detection.

    This function automatically detects the type of call and handles:
    1. Config file mode: fk("config.yaml") or fk(config="config.yaml")
    2. Simple mode: fk([1,-2,3], 2)

    Parameters
    ----------
    braid_or_config
        Either a braid list [1,-2,3] OR config file path "config.yaml"
    degree
        Computation degree (required unless using config file)
    config
        Config file path (alternative to passing as first argument)
    symbolic
        Generate symbolic polynomial representation using SymPy
    threads
        Number of threads for C++ FK computation
    name
        Name for saved files when save_data is enabled
    weight
        Optional weight parameter for stratified calculation
    **kwargs
        Additional parameters (verbose, max_workers, save_data, ...).

    Returns
    -------
    dict
        Dictionary containing comprehensive computation results.

    Raises
    ------
    SignAssignmentError
        If no valid sign assignment exists for the braid at this degree.

    Examples
    --------
    >>> fk([1,-2,3], 2)                      # Simple mode
    >>> fk([1,-2,3], 2, symbolic=True)       # With symbolic output
    >>> fk("config.yaml")                    # From configuration file
    """
    # Config file mode detection
    if isinstance(braid_or_config, str) or config is not None:
        config_path = config or braid_or_config
        if not isinstance(config_path, str):
            raise TypeError("Config path must be a string")
        return fk_from_config(config_path)

    # Simple mode - just braid and degree with defaults
    braid = braid_or_config
    if degree is None:
        raise ValueError("degree is required when providing a braid list")

    unknown = set(kwargs) - set(_DEFAULT_PARAMS)
    if unknown:
        raise TypeError(
            f"Unknown fk() parameter(s): {sorted(unknown)}. "
            f"Valid parameters: {sorted(_DEFAULT_PARAMS)}"
        )

    configure_logging(kwargs.get('verbose', False))

    params = dict(_DEFAULT_PARAMS)
    params['link_name'] = name
    params['symbolic'] = symbolic
    if threads is not None:
        params['threads'] = threads
    params['weight'] = weight
    params.update(kwargs)

    return _fk_compute(braid, degree, **params)


def _fk_compute(
    braid: List[int],
    degree: int,
    verbose: bool = True,
    max_workers: int = 1,
    chunk_size: int = 1 << 14,
    include_flip: bool = False,
    max_shifts: Optional[int] = None,
    save_data: bool = False,
    save_dir: str = "data",
    link_name: Optional[str] = None,
    symbolic: bool = False,
    threads: int = 1,
    ilp: Optional[str] = None,
    ilp_file: Optional[str] = None,
    inversion: Optional[Dict[str, Any]] = None,
    inversion_file: Optional[str] = None,
    partial_signs: Optional[List[int]] = None,
    weight: Optional[int] = None,
    _progress_callback: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Internal FK computation function.

    This function executes the complete FK computation pipeline:
    1. Sign assignment calculation (inversion step)
    2. ILP formulation
    3. FK invariant computation using compiled helper binary
    4. Optional symbolic polynomial conversion

    Intermediate files live next to the saved output when save_data is
    enabled; otherwise they are written to a private temporary directory
    that is removed afterwards (safe for concurrent runs).
    """
    with tempfile.TemporaryDirectory(prefix="fkcompute_") as tmpdir:
        if save_data:
            os.makedirs(save_dir, exist_ok=True)
            if not link_name:
                tstr = time.strftime("%H%M%S", time.localtime(time.time()))
                link_name = f"Unknown_{tstr}"
                logger.debug(f"No name provided for link/knot, saving as: {link_name}")
            else:
                logger.debug(f"Saving as: {link_name}")
            out_base = os.path.join(save_dir, link_name)
        else:
            out_base = os.path.join(tmpdir, "fk")

        # Step 1: Inversion data
        if ilp is None and ilp_file is None:
            logger.debug("No ilp data provided, I need to calculate it")

            inversion_task = None
            if _progress_callback:
                inversion_task = _progress_callback.start_inversion(braid, degree)

            if inversion is not None and inversion_file is not None:
                raise ValueError(
                    "Both inversion data and an inversion file were passed. "
                    "Pass one or the other."
                )

            if inversion_file is not None:
                logger.debug("Loading inversion from file")
                inversion = load_inversion_file(inversion_file)

            if inversion is None:
                logger.debug("Calculating inversion data")
                result = find_sign_assignment(
                    braid,
                    degree=degree,
                    partial_signs=partial_signs,
                    max_workers=max_workers,
                    chunk_size=chunk_size,
                    include_flip=include_flip,
                    max_shifts=max_shifts,
                    verbose=verbose,
                    weight=weight,
                )
                inversion = result.to_result_dict()
                logger.debug("Inversion data calculated")

            _assess_inversion(inversion, braid, degree)

            if _progress_callback and inversion_task is not None:
                if 'metadata' in inversion and 'components' in inversion['metadata']:
                    _progress_callback.complete_inversion(inversion['metadata']['components'])
                else:
                    _progress_callback.complete_inversion(len(braid))

            if save_data:
                with open(out_base + "_inversion.json", "w") as f:
                    json.dump(inversion, f)

        # Step 2: ILP data
        if ilp is not None and ilp_file is not None:
            raise ValueError(
                "Both ilp data and an ilp file were passed. Pass one or the other."
            )

        ilp_task = None
        if _progress_callback:
            ilp_task = _progress_callback.start_ilp_generation()

        if ilp is None and ilp_file is not None:
            logger.debug("Loading ILP from file")
            ilp = load_ilp_file(ilp_file)

        if ilp_file is None:
            ilp_file_base = f"{out_base}_ilp"
        else:
            ilp_file_base = os.path.splitext(ilp_file)[0]

        if ilp is None:
            logger.debug("Calculating ILP data")
            if inversion is None:
                raise ValueError("Inversion data is required for ILP calculation")
            ilp = _generate_ilp(
                braid,
                degree=degree,
                inversion_data=inversion["inversion_data"],
                outfile=ilp_file_base + ".csv",
                weight=weight,
            )
            if ilp is None:
                raise SignAssignmentError(
                    f"The sign assignment is not valid for braid {braid} at degree "
                    f"{degree}; the ILP system is unbounded or infeasible."
                )
            logger.debug("ILP data calculated")

            if _progress_callback and ilp_task is not None:
                estimated_constraints = len(braid) * degree * 10
                _progress_callback.complete_ilp_generation(estimated_constraints)
        elif ilp_file is None:
            # ILP passed as in-memory data: materialize it for the binary.
            with open(ilp_file_base + ".csv", "w") as f:
                f.write(ilp)

        # Step 3: FK invariant computation
        estimated_points = len(braid) * degree * 50

        fk_task = None
        if _progress_callback:
            fk_task = _progress_callback.start_fk_computation(threads, estimated_points)

        run_fk_binary(ilp_file_base, out_base, threads=threads, verbose=verbose)

        with open(out_base + ".json", "r") as f:
            fk_result = json.load(f)

        if _progress_callback and fk_task is not None:
            _progress_callback.update_fk_progress(estimated_points, estimated_points)
            _progress_callback.complete_fk_computation(len(fk_result.get("terms", [])))

        # Get the braid that was actually used
        computed_braid = inversion["braid"] if inversion and "braid" in inversion else braid

        # Extract number of components from the braid
        braid_states = BraidStates(computed_braid)
        components = braid_states.n_components

        # Prepare result dictionary
        if "metadata" not in fk_result:
            fk_result["metadata"] = {}
        fk_result["metadata"]["braid"] = computed_braid
        fk_result["metadata"]["inversion"] = inversion["inversion_data"] if inversion else None
        fk_result["metadata"]["components"] = components

        # Add symbolic representation if requested
        if symbolic and SYMPY_AVAILABLE:
            try:
                fk_result["metadata"]["symbolic"] = format_symbolic_output(fk_result, "pretty")
            except Exception as e:
                logger.warning(f"Could not generate symbolic representation: {e}")
        elif symbolic and not SYMPY_AVAILABLE:
            logger.warning(
                "SymPy not available for symbolic output. Install with: pip install sympy"
            )

        # Save updated metadata to file if save_data is enabled
        if save_data:
            with open(out_base + ".json", "w") as f:
                json.dump(fk_result, f, indent="\t")

        return fk_result


def _generate_ilp(braid: List[int], degree: int, inversion_data: Dict, outfile: str, weight: Optional[int] = None) -> Optional[str]:
    """Generate ILP formulation for FK computation."""
    from ..domain.braid.word import is_homogeneous_braid
    from ..domain.constraints.reduction import full_reduce
    from ..solver.ilp import ilp as generate_ilp

    braid_states = BraidStates(braid)

    if not is_homogeneous_braid(braid):
        # Fibered braid - load the provided inversion data
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()
        if not braid_states.validate():
            raise SignAssignmentError(
                f"The inversion data is not valid for braid {braid}."
            )
        braid_states.generate_position_assignments()

    relations = full_reduce(braid_states.get_state_relations())
    return generate_ilp(degree, relations, braid_states, outfile, weight=weight)
