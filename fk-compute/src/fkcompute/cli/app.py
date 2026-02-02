"""
Typer application definition for FK computation CLI.

This module defines the main Typer app and its configuration.
"""

from __future__ import annotations

import sys
from typing import Optional, List

try:
    import typer
except ImportError:
    print("Error: typer is required. Please install it with: pip install typer")
    sys.exit(1)

# Create the main Typer app
app = typer.Typer(
    help="""
Compute the FK invariant for braids using inversion, ILP reduction, and compiled helper binary.

The FK invariant is a mathematical object used in knot theory to distinguish different
types of knots and links. This tool provides multiple interfaces:

INTERACTIVE MODE:
  fk                                  # Start interactive mode (prompted for all parameters)
  fk interactive                      # Same as above

SIMPLE USAGE:
  fk simple "[1,-2,3]" 2              # Quick computation with defaults

CONFIG FILE USAGE:
  fk config single.yaml               # Single computation from file
  fk config batch.json                # Multiple computations from file
  fk config c1.yaml c2.yaml c3.json   # Process multiple config files

TEMPLATE CREATION:
  fk template create my_config.yaml   # Create a blank configuration template

BRAID FORMATS:
  Braids can be specified in multiple formats:
  - JSON-style: "[1, -2, -3, 1]"
  - Comma-separated: "1,-2,-3,1"
  - Space-separated: "1 -2 -3 1"

For detailed documentation, run 'man fk' (after running 'fk-install-man').
""",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

# Sub-apps
template_app = typer.Typer(help="Template management commands")
history_app = typer.Typer(help="History and session management commands")

app.add_typer(template_app, name="template")
app.add_typer(history_app, name="history")


def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the fk CLI tool."""
    from .commands import handle_interactive

    if argv is None:
        argv = sys.argv

    # Handle legacy simple mode: fk "[1,-2,3]" 2 -> fk simple "[1,-2,3]" 2
    if len(argv) >= 3 and argv[1] not in [
        "simple", "config", "interactive", "template", "history", "-h", "--help", "--version"
    ]:
        try:
            int(argv[2])
            new_argv = [argv[0], "simple"] + argv[1:]
            original_argv = sys.argv[:]
            sys.argv = new_argv
            try:
                app()
            finally:
                sys.argv = original_argv
            return
        except (ValueError, IndexError):
            pass

    # Check if no arguments provided - default to interactive mode
    if len(argv) == 1:
        handle_interactive()
        return

    # Run typer normally
    app()
