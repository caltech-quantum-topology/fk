"""
Batch processing for FK computation.

This module provides functions for processing multiple FK computations
from configuration files.
"""

import logging
from typing import Any, Dict, List, Union

from .presets import PRESETS
from ..infra.config import load_config_file

logger = logging.getLogger("fk_logger")


def _wrap_inversion(config: Dict[str, Any], braid: List[int], degree: int) -> None:
    """Convert a config-style inversion dict into the pipeline format, in place."""
    if config.get("inversion") is not None:
        inversion_data = {int(k): v for k, v in config["inversion"].items()}
        config["inversion"] = {
            "inversion_data": inversion_data,
            "braid": braid,
            "degree": degree,
        }


def _apply_preset(config: Dict[str, Any], context: str = "") -> Dict[str, Any]:
    """Merge a named preset (if any) under the explicit config values."""
    preset = config.get("preset")
    if not preset:
        return {k: v for k, v in config.items() if k != "preset"}
    if preset not in PRESETS:
        raise ValueError(
            f"Unknown preset '{preset}'{context}. Available: {list(PRESETS.keys())}"
        )
    merged = PRESETS[preset].copy()
    merged.update({k: v for k, v in config.items() if k != "preset"})
    return merged


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

    _wrap_inversion(config_data, braid, degree)

    # Use 'name' parameter as 'link_name' if provided
    name = config_data.get("name")
    if name and "link_name" not in config_data:
        config_data["link_name"] = name

    # Import compute function here to avoid circular imports
    from .compute import _fk_compute, configure_logging

    comp_config = _apply_preset(
        {k: v for k, v in config_data.items() if k not in ("braid", "degree", "name")}
    )
    configure_logging(comp_config.get("verbose", False))
    return _fk_compute(braid, degree, **comp_config)


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
        k: v for k, v in config_data.items() if k != "computations"
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
            if k not in ("name", "braid", "degree")
        })

        # Use computation name as link_name if not explicitly set
        if computation.get("name") and "link_name" not in comp_config:
            comp_config["link_name"] = computation["name"]

        comp_config = _apply_preset(comp_config, context=f" for computation '{comp_name}'")

        verbose = comp_config.get("verbose", False)
        configure_logging(verbose)

        _wrap_inversion(comp_config, braid, degree)

        if total > 1 and verbose:
            logger.info(f"Computing {comp_name} ({i}/{total})")

        try:
            result = _fk_compute(
                braid,
                degree,
                **{k: v for k, v in comp_config.items() if k not in ("braid", "degree")}
            )
            results[comp_name] = result

        except Exception as e:
            error_msg = f"Failed to compute {comp_name}: {str(e)}"
            if verbose:
                logger.error(error_msg)
            results[comp_name] = {"error": str(e)}

    return results
