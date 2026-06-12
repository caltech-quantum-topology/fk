"""
Batch processing for FK computation.

This module provides functions for processing multiple FK computations
from configuration files.
"""

import logging
from typing import Any, Dict, Union

from .presets import PRESETS
from ..infra.config import load_config_file, parse_int_list

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


def fk_from_config(config_path: str) -> Union[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """
    Handle FK computation from JSON/YAML configuration files.

    Supports both single computation and batch processing modes.

    Parameters
    ----------
    config_path
        Path to JSON or YAML configuration file.

    Returns
    -------
    dict
        Single computation: Dict[str, Any] containing FK results
        Batch processing: Dict[str, Dict[str, Any]] keyed by computation names
    """
    config_data = load_config_file(config_path)

    # Check if this is a batch configuration
    if "computations" in config_data:
        return fk_batch_from_config(config_data, config_path)

    # Single computation mode
    braid = config_data.get("braid")
    if not braid:
        raise ValueError("'braid' is required in config file")

    degree = config_data.get("degree")
    if degree is None:
        raise ValueError("'degree' is required in config file")

    # Process inversion data if present
    if "inversion" in config_data and config_data["inversion"] is not None:
        inversion_dict = config_data["inversion"]
        inversion_data = {int(k): v for k, v in inversion_dict.items()}
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

    # Import compute function here to avoid circular imports
    from .compute import _fk_compute, configure_logging

    if preset:
        filtered_config = {
            k: v for k, v in config_data.items()
            if k not in ["braid", "degree", "preset", "name"]
        }
        preset_config = PRESETS.get(preset, {}).copy()
        preset_config.update(filtered_config)
        verbose = preset_config.get("verbose", False)
        configure_logging(verbose)
        return _fk_compute(braid, degree, **preset_config)
    else:
        filtered_config = {
            k: v for k, v in config_data.items()
            if k not in ["braid", "degree", "preset", "name"]
        }
        verbose = filtered_config.get("verbose", False)
        configure_logging(verbose)
        return _fk_compute(braid, degree, **filtered_config)


def fk_batch_from_config(
    config_data: Dict[str, Any],
    config_path: str
) -> Dict[str, Dict[str, Any]]:
    """
    Execute batch FK computations from configuration data.

    Parameters
    ----------
    config_data
        Dictionary containing batch configuration with 'computations' array.
    config_path
        Path to original config file (for error reporting).

    Returns
    -------
    dict
        Results keyed by computation names.
    """
    from .compute import _fk_compute, configure_logging

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
        comp_name = computation.get("name", f"computation_{i}")

        braid = computation.get("braid")
        if not braid:
            raise ValueError(f"'braid' is required for computation '{comp_name}'")

        degree = computation.get("degree")
        if degree is None:
            raise ValueError(f"'degree' is required for computation '{comp_name}'")

        # Merge global defaults with computation-specific parameters
        comp_config = global_defaults.copy()
        comp_config.update({
            k: v for k, v in computation.items()
            if k not in ["name", "braid", "degree"]
        })

        # Use computation name as link_name if not explicitly set
        comp_name_from_config = computation.get("name")
        if comp_name_from_config and "link_name" not in comp_config:
            comp_config["link_name"] = comp_name_from_config

        # Handle preset if specified
        preset = comp_config.get("preset")
        if preset:
            if preset not in PRESETS:
                raise ValueError(
                    f"Unknown preset '{preset}' for computation '{comp_name}'. "
                    f"Available: {list(PRESETS.keys())}"
                )
            filtered_config = {k: v for k, v in comp_config.items() if k != "preset"}
            preset_config = PRESETS[preset].copy()
            preset_config.update(filtered_config)
            comp_config = preset_config

        # Configure logging
        verbose = comp_config.get("verbose", False)
        if not verbose and total > 1:
            comp_config["verbose"] = False

        configure_logging(verbose)

        # Process inversion data if present
        if "inversion" in comp_config and comp_config["inversion"] is not None:
            inversion_dict = comp_config["inversion"]
            inversion_data = {int(k): v for k, v in inversion_dict.items()}
            comp_config["inversion"] = {
                "inversion_data": inversion_data,
                "braid": braid,
                "degree": degree,
            }

        if total > 1 and verbose:
            logger.info(f"Computing {comp_name} ({i}/{total})")

        try:
            result = _fk_compute(
                braid,
                degree,
                **{k: v for k, v in comp_config.items() if k not in ["braid", "degree", "preset"]}
            )
            results[comp_name] = result

        except Exception as e:
            error_msg = f"Failed to compute {comp_name}: {str(e)}"
            if verbose:
                logger.error(error_msg)
            results[comp_name] = {"error": str(e)}

    return results
