"""
Core FK computation function.

This module provides the main fk() function for computing FK invariants.
"""

import json
import logging
import os
import time
from typing import Any, Dict, List, Optional, Union

from .presets import PRESETS
from .batch import fk_from_config
from ..domain.braid.states import BraidStates
from ..solver.inversion import get_sign_assignment_parallel
from ..solver.ilp import ilp
from ..infra.binary import binary_path, safe_unlink, run_fk_binary
from ..infra.config import load_inversion_file, load_ilp_file
from ..output.symbolic import format_symbolic_output, SYMPY_AVAILABLE

# Logger setup
logger = logging.getLogger("fk_logger")


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


def _assess_inversion(inversion: Dict[str, Any]) -> None:
    """Validate the structure of inversion data."""
    if "inversion_data" not in inversion:
        logger.error("Inversion data missing in inversion dictionary")
    if inversion.get("inversion_data") == "failure":
        logger.error("Inversion data invalid (failure)")


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
        Additional parameters (verbose, max_workers, etc.)

    Returns
    -------
    dict
        Dictionary containing comprehensive computation results.

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

    configure_logging(kwargs.get('verbose', False))

    # Build parameters with defaults and overrides
    params = {
        'verbose': False,
        'max_workers': 1,
        'chunk_size': 1 << 14,
        'include_flip': False,
        'max_shifts': None,
        'save_data': False,
        'save_dir': "data",
        'link_name': name,
        'symbolic': symbolic,
        'threads': threads if threads is not None else 1,
        'ilp': None,
        'ilp_file': None,
        'inversion': None,
        'inversion_file': None,
        'partial_signs': None,
        'weight': weight,
    }

    # Override with any provided kwargs
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
    """
    # Handle naming and save directories
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

    # Step 1: Inversion data
    if ilp is None and ilp_file is None:
        logger.debug("No ilp data provided, I need to calculate it")

        inversion_task = None
        if _progress_callback:
            inversion_task = _progress_callback.start_inversion(braid, degree)

        if inversion is not None and inversion_file is not None:
            logger.error("Both inversion data and inversion file are passed through. Pass one or the other")

        if inversion_file is not None and inversion is None:
            logger.debug("Loading inversion from file")
            inversion = load_inversion_file(inversion_file)

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

        if inversion is not None:
            _assess_inversion(inversion)

            if _progress_callback and inversion_task is not None:
                if 'metadata' in inversion and 'components' in inversion['metadata']:
                    _progress_callback.complete_inversion(inversion['metadata']['components'])
                else:
                    _progress_callback.complete_inversion(len(braid))

    if save_data:
        with open(link_name + "_inversion.json", "w") as f:
            json.dump(inversion, f)

    # Step 2: ILP data
    if ilp is not None and ilp_file is not None:
        logger.error("Both ilp data and ilp file are passed through. Pass one or the other")

    ilp_task = None
    if _progress_callback:
        ilp_task = _progress_callback.start_ilp_generation()

    if ilp is None and ilp_file is not None:
        logger.debug("Loading ILP from file")
        ilp = load_ilp_file(ilp_file)

    if ilp is None and ilp_file is None:
        logger.debug("Calculating ILP data")
        if inversion is None:
            raise ValueError("Inversion data is required for ILP calculation")
        ilp = _generate_ilp(
            braid,
            degree=degree,
            inversion_data=inversion["inversion_data"],
            outfile=link_name + "_ilp.csv",
            weight=weight,
        )
        logger.debug("ILP data calculated")

        if _progress_callback and ilp_task is not None:
            estimated_constraints = len(braid) * degree * 10
            _progress_callback.complete_ilp_generation(estimated_constraints)

    # Step 3: FK invariant computation
    if ilp_file is None:
        ilp_file_base = f"{link_name}_ilp"
    else:
        ilp_file_base = ilp_file[:ilp_file.index(".")]

    estimated_points = len(braid) * degree * 50

    fk_task = None
    if _progress_callback:
        fk_task = _progress_callback.start_fk_computation(threads, estimated_points)

    run_fk_binary(ilp_file_base, link_name, threads=threads, verbose=verbose)

    if _progress_callback and fk_task is not None:
        _progress_callback.update_fk_progress(estimated_points, estimated_points)
        try:
            with open(link_name + ".json", "r") as f:
                result_data = json.load(f)
                terms_count = len(result_data.get("terms", []))
                _progress_callback.complete_fk_computation(terms_count)
        except:
            _progress_callback.complete_fk_computation(0)

    with open(link_name + ".json", "r") as f:
        fk_result = json.load(f)

    # Clean up intermediate files unless explicitly saving
    if not save_data:
        safe_unlink(f"{link_name}_inversion.csv")
        safe_unlink(f"{link_name}_ilp.csv")
        safe_unlink(f"{link_name}.json")

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
            symbolic_repr = format_symbolic_output(fk_result, "pretty")
            if "metadata" in fk_result:
                fk_result["metadata"]["symbolic"] = symbolic_repr
            else:
                fk_result["symbolic"] = symbolic_repr
        except Exception as e:
            if verbose:
                print(f"Warning: Could not generate symbolic representation: {e}")
    elif symbolic and not SYMPY_AVAILABLE:
        if verbose:
            print("Warning: SymPy not available for symbolic output. Install with: pip install sympy")

    # Save updated metadata to file if save_data is enabled
    if save_data:
        with open(link_name + ".json", "w") as f:
            json.dump(fk_result, f, indent="\t")

    return fk_result


def _generate_ilp(braid: List[int], degree: int, inversion_data: Dict, outfile: str, weight: Optional[int] = None) -> Optional[str]:
    """Generate ILP formulation for FK computation."""
    from ..domain.braid.states import BraidStates
    from ..domain.braid.word import is_homogeneous_braid
    from ..domain.constraints.reduction import full_reduce
    from ..solver.ilp import ilp as generate_ilp

    braid_states = BraidStates(braid)

    if is_homogeneous_braid(braid):
        # Homogeneous braid - use default sign assignment
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
        return generate_ilp(degree, relations, braid_states, outfile, weight=weight)
    else:
        # Fibered braid - use provided inversion data
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()
        if braid_states.validate():
            braid_states.generate_position_assignments()
            all_relations = braid_states.get_state_relations()
            relations = full_reduce(all_relations)
            return generate_ilp(degree, relations, braid_states, outfile, weight=weight)
        else:
            raise Exception("The inversion data doesn't seem to be valid.")
