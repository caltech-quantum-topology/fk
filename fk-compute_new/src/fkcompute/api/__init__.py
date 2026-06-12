"""
Application layer for FK computation.

This subpackage contains the main API:
- compute: Core fk() function
- batch: Batch processing
- presets: PRESETS dictionary
"""

from .compute import fk
from .batch import fk_batch_from_config, fk_from_config
from .presets import PRESETS

__all__ = [
    "fk",
    "fk_batch_from_config",
    "fk_from_config",
    "PRESETS",
]
