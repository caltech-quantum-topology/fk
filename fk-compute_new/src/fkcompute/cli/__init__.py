"""
Command-line interface for FK computation.

This subpackage contains:
- app: Typer app definition
- commands: All CLI commands
"""

from .app import app, main
from .commands import (
    interactive_command,
    simple_command,
    config_command,
    template_create_command,
)

__all__ = [
    "app",
    "main",
    "interactive_command",
    "simple_command",
    "config_command",
    "template_create_command",
]
