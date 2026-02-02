import logging
import sys
import os
import json
import time
import subprocess
import shutil
import pathlib
from typing import Optional, Dict, Any, List, Union
import argparse
from importlib import resources

HAS_YAML = False
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# Import symbolic output functionality
from .symbolic_output import (
    format_symbolic_output,
    print_symbolic_result,
    SYMPY_AVAILABLE,
)

# Import BraidStates for extracting component information
from .braidstates_links import BraidStates

from .sign_diagram_links_parallel_ansatz import get_sign_assignment_parallel
from .fklinks import FK as fk_ilp

# -------------------------------------------------------------------------
# Logger setup
# -------------------------------------------------------------------------
logger = logging.getLogger("fk_logger")


def configure_logging(verbose: bool) -> None:
    """Set logging level based on verbose flag."""
    handler = logging.StreamHandler()
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(fmt)

    # Remove old handlers so repeated runs don't duplicate output
    logger.handlers.clear()
    logger.addHandler(handler)

    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)  # Suppress INFO and DEBUG messages


# -------------------------------------------------------------------------
# Binary and file utilities
# -------------------------------------------------------------------------
def _binary_path(name: str) -> str:
    """Find the path to a required binary, checking PATH first then package."""
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
    """Safely remove a file, ignoring errors."""
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
    inversion["inversion_data"] = {
        int(k): v for k, v in inversion["inversion_data"].items()
    }
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


def _load_config_file(config_path: str) -> Dict[str, Any]:
    """Load configuration from JSON or YAML file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        if config_path.endswith((".yml", ".yaml")):
            if not HAS_YAML:
                raise ImportError(
                    "PyYAML is required for YAML config files. Install with: pip install PyYAML"
                )
            import yaml
            return yaml.safe_load(f)
        else:
            return json.load(f)


def _parse_bool(default_true: bool = False):
    """Legacy function for backward compatibility.
    
    Note: This function was used for argparse. With typer, boolean flags
    are handled directly by typer.Option, so this function is kept only
    for backward compatibility but should be replaced in new code.
    """
    # For typer, boolean flags are handled directly in function signatures
    # This function exists only for backward compatibility
    if default_true:
        return lambda parser, name, help_text: None
    else:
        return lambda parser, name, help_text: None


# -------------------------------------------------------------------------
# Preset configurations
# -------------------------------------------------------------------------
PRESETS = {
    "fast": {
        "max_workers": 1,
        "chunk_size": 1 << 12,  # Smaller chunks for faster start
        "include_flip": False,  # Disable flip symmetry
        "max_shifts": None,  # Don't limit shifts
        "verbose": False,
        "save_data": False,
        "threads": 1,
    },
    "accurate": {
        "max_workers": 4,
        "chunk_size": 1 << 16,  # Larger chunks for thoroughness
        "include_flip": False,  # Disable flip symmetry
        "max_shifts": None,  # No limit on shifts
        "verbose": True,
        "save_data": True,
        "threads": 1,
    },
    "parallel": {
        "max_workers": 8,  # High parallelism
        "chunk_size": 1 << 14,  # Balanced chunk size
        "include_flip": False,
        "max_shifts": None,
        "verbose": True,
        "save_data": False,
        "threads": 4,  # Use multiple threads for parallel preset
    },
}

# -------------------------------------------------------------------------
# Intelligent FK function that handles all interfaces
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Unified FK function - handles all interfaces intelligently
# -------------------------------------------------------------------------
def fk(
    braid_or_config: Union[List[int], str],
    degree: Optional[int] = None,
    config: Optional[str] = None,
    symbolic: bool = False,
    threads: Optional[int] = None,
    name: Optional[str] = None,
    **kwargs,
) -> Union[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """
    Unified FK invariant computation with intelligent interface detection.

    This function automatically detects the type of call and handles:
    1. Config file mode: fk("config.yaml") or fk(config="config.yaml")
    2. Simple mode: fk([1,-2,3], 2)

    Args:
        braid_or_config: Either a braid list [1,-2,3] OR config file path "config.yaml"
        degree: Computation degree (required unless using config file)
        config: Config file path (alternative to passing as first argument)
        symbolic: Generate symbolic polynomial representation using SymPy (default: False)
        threads: Number of threads for C++ FK computation (default: 1)
        name: Name for saved files when save_data is enabled (default: None)
        **kwargs: Additional parameters passed to _fk_compute (verbose, max_workers, etc.)

    Returns:
        Dictionary containing comprehensive computation results:
        {
            "braid": [1, -2, 3],                    # The braid used in computation
            "inversion_data": {...},                # Inversion assignment data
            "degree": 2,                            # Computation degree
            "components": 4,                        # Number of braid components/strands
            "fk": [[coeff_pairs], ...],             # FK coefficient matrix
            "symbolic": "polynomial_string"         # Symbolic representation (if requested)
        }

        For batch processing, returns a dictionary keyed by computation names:
        {
            "comp1": {"braid": [...], "inversion_data": {...}, "degree": ..., "components": ..., "fk": [...]},
            "comp2": {"braid": [...], "inversion_data": {...}, "degree": ..., "components": ..., "fk": [...]}
        }

    Examples:
        fk([1,-2,3], 2)                              # Simple mode
        fk([1,-2,3], 2, symbolic=True)              # With symbolic polynomial output
        fk([1,-2,3], 2, threads=4)                  # With 4 threads
        fk("config.yaml")                            # From configuration file

    Note:
        - Symbolic output requires SymPy (install with: pip install sympy)
        - The "components" field indicates the number of strands in the braid
        - FK coefficients are organized by powers of topological variables (x, y, etc.) and q
        - All advanced parameters (max_workers, verbose, etc.) can be passed as kwargs or set via config files
    """

    # ========== INTERFACE DETECTION ==========

    # 1. Config file mode detection
    if isinstance(braid_or_config, str) or config is not None:
        config_path = config or braid_or_config
        if not isinstance(config_path, str):
            raise TypeError("Config path must be a string")
        return _fk_from_config(config_path)

    # 2. Simple mode - just braid and degree with defaults
    braid = braid_or_config
    if degree is None:
        raise ValueError("degree is required when providing a braid list")

    configure_logging(kwargs.get('verbose', False))  # Quiet by default for simple mode

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
    }

    # Override with any provided kwargs
    params.update(kwargs)

    return _fk_compute(braid, degree, **params)


def _fk_from_config(
    config_path: str,
) -> Union[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """
    Handle FK computation from JSON/YAML configuration files.

    Supports both single computation and batch processing modes:
    - Single mode: Direct computation from config with braid, degree, and parameters
    - Batch mode: Multiple computations with global defaults and per-computation overrides

    Args:
        config_path: Path to JSON or YAML configuration file

    Returns:
        Single computation: Dict[str, Any] containing FK results
        Batch processing: Dict[str, Dict[str, Any]] keyed by computation names

    Configuration Format:
        Single computation:
        {
            "braid": [1, -2, 3],
            "degree": 2,
            "name": "my_knot",
            "preset": "accurate",
            "max_workers": 8,
            "save_data": true,
            "ilp_file": "path/to/precomputed.ilp",
            "inversion": {
                "0": [1, -1, 1],
                "1": [-1, 1]
            }
        }

        Batch processing:
        {
            "preset": "fast",
            "max_workers": 4,
            "save_data": true,
            "computations": [
                {
                    "name": "trefoil",
                    "braid": [1, 1, 1],
                    "degree": 2
                },
                {
                    "name": "figure_eight",
                    "braid": [1, -2, -1, -2],
                    "degree": 3,
                    "preset": "accurate"
                }
            ]
        }

    Note:
        - Global parameters are applied to all batch computations
        - Individual computations can override global settings
        - All presets disable flip symmetry by default
        - Batch mode provides progress tracking for multiple computations
        - Inversion data: just provide the component->signs dictionary; keys are auto-converted to integers
        - The 'name' parameter is used for file naming when 'save_data' is true
    """
    config_data = _load_config_file(config_path)

    # Check if this is a batch configuration (multiple braids)
    if "computations" in config_data:
        return _fk_batch_from_config(config_data, config_path)

    # Single computation mode (existing behavior)
    braid = config_data.get("braid")
    if not braid:
        raise ValueError("'braid' is required in config file")

    degree = config_data.get("degree")
    if degree is None:
        raise ValueError("'degree' is required in config file")

    # Process inversion data if present - convert to proper internal structure
    if "inversion" in config_data and config_data["inversion"] is not None:
        inversion_dict = config_data["inversion"]
        # Convert string keys to integers and wrap in proper structure
        inversion_data = {
            int(k): v for k, v in inversion_dict.items()
        }
        config_data["inversion"] = {
            "inversion_data": inversion_data,
            "braid": braid,
            "degree": degree,
        }

    # Use 'name' parameter as 'link_name' if provided
    name = config_data.get("name")
    if name and "link_name" not in config_data:
        config_data["link_name"] = name

    # Check if using preset in config
    preset = config_data.get("preset")
    if preset:
        filtered_config = {
            k: v
            for k, v in config_data.items()
            if k not in ["braid", "degree", "preset", "name"]
        }
        preset_config = PRESETS.get(preset, {}).copy()
        preset_config.update(filtered_config)
        verbose = preset_config.get("verbose", False)
        configure_logging(verbose)
        return _fk_compute(braid, degree, **preset_config)
    else:
        filtered_config = {
            k: v for k, v in config_data.items() if k not in ["braid", "degree", "name"]
        }
        verbose = filtered_config.get("verbose", False)
        configure_logging(verbose)
        return _fk_compute(braid, degree, **filtered_config)


def _fk_batch_from_config(
    config_data: Dict[str, Any], config_path: str
) -> Dict[str, Dict[str, Any]]:
    """
    Execute batch FK computations from configuration data with progress tracking.

    Processes multiple braid computations with global defaults and per-computation
    overrides, providing organized results keyed by computation names.

    Args:
        config_data: Dictionary containing batch configuration with 'computations' array
        config_path: Path to original config file (for error reporting)

    Returns:
        Dict[str, Dict[str, Any]]: Results keyed by computation names, each containing
        standard FK computation output (braid, inversion_data, degree, components, fk)

    Batch Processing Features:
        - Global parameter defaults applied to all computations
        - Individual computation parameter overrides
        - Named computations for organized output
        - Progress tracking during batch execution
        - Error handling continues processing remaining computations

    Note:
        - Computation names default to "computation_N" if not specified
        - Failed computations are logged but don't halt the batch
        - All computations use the same global logging configuration
    """
    computations = config_data.get("computations", [])
    if not computations:
        raise ValueError("'computations' array is empty in config file")

    # Global defaults from the config file
    global_defaults = {
        k: v for k, v in config_data.items() if k not in ["computations"]
    }

    results = {}
    total = len(computations)

    for i, computation in enumerate(computations, 1):
        # Get computation name (for results key)
        comp_name = computation.get("name", f"computation_{i}")

        # Required fields for each computation
        braid = computation.get("braid")
        if not braid:
            raise ValueError(f"'braid' is required for computation '{comp_name}'")

        degree = computation.get("degree")
        if degree is None:
            raise ValueError(f"'degree' is required for computation '{comp_name}'")

        # Merge global defaults with computation-specific parameters
        comp_config = global_defaults.copy()
        comp_config.update(
            {
                k: v
                for k, v in computation.items()
                if k not in ["name", "braid", "degree"]
            }
        )

        # Use computation name as link_name if not explicitly set
        comp_name_from_config = computation.get("name")
        if comp_name_from_config and "link_name" not in comp_config:
            comp_config["link_name"] = comp_name_from_config

        # Handle preset if specified
        preset = comp_config.get("preset")
        if preset:
            if preset not in PRESETS:
                raise ValueError(
                    f"Unknown preset '{preset}' for computation '{comp_name}'. Available: {list(PRESETS.keys())}"
                )

            # Start with preset, then apply overrides
            filtered_config = {k: v for k, v in comp_config.items() if k != "preset"}
            preset_config = PRESETS[preset].copy()
            preset_config.update(filtered_config)
            comp_config = preset_config

        # Configure logging (quiet for batch unless explicitly verbose)
        verbose = comp_config.get("verbose", False)
        if not verbose and total > 1:
            # For batch jobs, be quiet unless explicitly requested
            comp_config["verbose"] = False

        configure_logging(verbose)

        # Process inversion data if present - convert to proper internal structure
        if "inversion" in comp_config and comp_config["inversion"] is not None:
            inversion_dict = comp_config["inversion"]
            # Convert string keys to integers and wrap in proper structure
            inversion_data = {
                int(k): v for k, v in inversion_dict.items()
            }
            comp_config["inversion"] = {
                "inversion_data": inversion_data,
                "braid": braid,
                "degree": degree,
            }

        # Print progress for batch jobs
        if total > 1 and verbose:
            logger.info(f"Computing {comp_name} ({i}/{total})")
        """
        I need to figure out how to do this without messing up redirecting to file
        elif total > 1:
            print(f"Computing {comp_name} ({i}/{total})")
        """

        try:
            # Run the computation
            result = _fk_compute(
                braid,
                degree,
                **{
                    k: v for k, v in comp_config.items() if k not in ["braid", "degree"]
                },
            )
            results[comp_name] = result

        except Exception as e:
            error_msg = f"Failed to compute {comp_name}: {str(e)}"
            if verbose:
                logger.error(error_msg)
            results[comp_name] = {"error": str(e)}

    return results


# -------------------------------------------------------------------------
# Core computation function (internal)
# -------------------------------------------------------------------------
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
    _progress_callback: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Internal FK computation function performing the complete FK invariant calculation.

    This function executes the complete FK computation pipeline:
    1. Sign assignment calculation (inversion step) - can be parallelized
    2. ILP (Integer Linear Programming) formulation
    3. FK invariant computation using compiled helper binary
    4. Component count extraction from braid topology
    5. Optional symbolic polynomial conversion using SymPy

    Args:
        braid: List of integers representing braid word (e.g., [1, -2, 3])
        degree: Degree of the FK invariant to compute
        verbose: Enable verbose logging output (default: True)
        max_workers: Number of parallel workers for inversion calculation (default: 1)
        chunk_size: Chunk size for parallel processing tasks (default: 16384)
        include_flip: Include flip symmetry in inversion calculation (default: False)
        max_shifts: Maximum number of cyclic shifts to consider (default: None = unlimited)
        save_data: Save intermediate computation files (inversion/ILP/JSON) (default: False)
        save_dir: Directory path for saving intermediate files (default: "data")
        link_name: Base name for saved files, auto-generated if None (default: None)
        symbolic: Add symbolic polynomial representation using SymPy (default: False)
        threads: Number of threads for C++ FK computation (default: 1)
        ilp: Pre-computed ILP data string to skip ILP calculation (default: None)
        ilp_file: Path to file containing pre-computed ILP data (default: None)
        inversion: Pre-computed inversion data dictionary (default: None)
        inversion_file: Path to JSON file with pre-computed inversion data (default: None)
        partial_signs: Optional partial sign assignments for constrained inversion (default: None)

    Returns:
        Dict[str, Any]: Comprehensive computation results containing:
            - "braid": Original braid list used in computation
            - "inversion_data": Sign assignment data from inversion step
            - "degree": Computation degree used
            - "components": Number of components/strands in the braid topology
            - "fk": FK coefficient matrix as nested lists of [power, coefficient] pairs
            - "symbolic": Human-readable polynomial string (only if symbolic=True and SymPy available)

    Raises:
        ValueError: If braid list is empty or degree is invalid
        FileNotFoundError: If specified inversion_file or ilp_file doesn't exist
        ImportError: If symbolic=True but SymPy is not available

    Note:
        - This is an internal function - use the main fk() function for external calls
        - Flip symmetry is disabled by default for improved performance
        - The components count reflects braid topology, not just maximum absolute value
        - Symbolic output requires SymPy: pip install sympy
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
        
        # Start progress tracking if callback provided
        inversion_task = None
        if _progress_callback:
            inversion_task = _progress_callback.start_inversion(braid, degree)

        if inversion is not None and inversion_file is not None:
            logger.error(
                "Both inversion data and inversion file are passed through. Pass one or the other"
            )

        if inversion_file is not None and inversion is None:
            logger.debug("Loading inversion from file")
            inversion = _load_inversion_file(inversion_file)

        if inversion is None and inversion_file is None:
            logger.debug("Calculating inversion data")
            inversion = get_sign_assignment_parallel(
                braid,
                partial_signs=partial_signs,  # type: ignore
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
            
            # Complete inversion progress if tracking
            if _progress_callback and inversion_task is not None:
                if 'metadata' in inversion and 'components' in inversion['metadata']:
                    _progress_callback.complete_inversion(inversion['metadata']['components'])
                else:
                    _progress_callback.complete_inversion(len(braid))

    if save_data:
        with open(link_name + "_inversion.json", "w") as f:
            json.dump(inversion, f)

    # --- Step 2: ILP data ---
    if ilp is not None and ilp_file is not None:
        logger.error(
            "Both ilp data and ilp file are passed through. Pass one or the other"
        )

    # Start ILP progress if tracking
    ilp_task = None
    if _progress_callback:
        ilp_task = _progress_callback.start_ilp_generation()

    if ilp is None and ilp_file is not None:
        logger.debug("Loading ILP from file")
        ilp = _load_ilp_file(ilp_file)

    if ilp is None and ilp_file is None:
        logger.debug("Calculating ILP data")
        if inversion is None:
            raise ValueError("Inversion data is required for ILP calculation")
        ilp = fk_ilp(
            braid,
            degree=degree,
            inversion_data=inversion["inversion_data"],
            outfile=link_name + "_ilp.csv",
        )
        logger.debug("ILP data calculated")
        
        # Complete ILP progress if tracking
        if _progress_callback and ilp_task is not None:
            # Estimate constraints (this is approximate)
            estimated_constraints = len(braid) * degree * 10  # Rough estimate
            _progress_callback.complete_ilp_generation(estimated_constraints)

    # --- Step 3: FK invariant computation ---
    bin_path = _binary_path("fk_main")
    if ilp_file is None:
        ilp_file = f"{link_name}_ilp"
    else:
        ilp_file = ilp_file[:ilp_file.index(".")]
    
    # Estimate total points for progress tracking
    estimated_points = len(braid) * degree * 50  # Rough estimate
    
    # Start FK computation progress if tracking
    fk_task = None
    if _progress_callback:
        fk_task = _progress_callback.start_fk_computation(threads, estimated_points)
    
    cmd = [bin_path, ilp_file, f"{link_name}", "--threads", str(threads)]
    if verbose:
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    
    # Complete FK progress if tracking
    if _progress_callback and fk_task is not None:
        _progress_callback.update_fk_progress(estimated_points, estimated_points)  # Mark as complete
        # Get actual term count from result file for accurate completion
        try:
            with open(link_name + ".json", "r") as f:
                result_data = json.load(f)
                terms_count = len(result_data.get("terms", []))
                _progress_callback.complete_fk_computation(terms_count)
        except:
            _progress_callback.complete_fk_computation(0)  # Fallback

    with open(link_name + ".json", "r") as f:
        fk_result = json.load(f)

    # Clean up intermediate files unless explicitly saving
    if not save_data:
        _safe_unlink(f"{link_name}_inversion.csv")
        _safe_unlink(f"{link_name}_ilp.csv")
        _safe_unlink(f"{link_name}.json")

    # Get the braid that was actually used (may be modified during inversion computation)
    # If inversion was computed, the braid might have been canonicalized
    # If no inversion was computed (e.g., when ILP file provided), use original braid
    computed_braid = inversion["braid"] if inversion and "braid" in inversion else braid

    # Extract number of components (strands) from the braid
    braid_states = BraidStates(computed_braid)
    components = braid_states.n_components

    # Prepare result dictionary
    if "metadata" not in fk_result:
        fk_result["metadata"] = {}
    fk_result["metadata"]["braid"] = computed_braid
    fk_result["metadata"]["inversion"] = inversion["inversion_data"] if inversion else None
    fk_result["metadata"]["components"] = components

    # Add symbolic representation if requested and Sympy is available
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
            print(
                "Warning: SymPy not available for symbolic output. Install with: pip install sympy"
            )

    # Save updated metadata to file if save_data is enabled
    if save_data:
        with open(link_name + ".json", "w") as f:
            json.dump(fk_result, f, indent="\t")

    return fk_result 
