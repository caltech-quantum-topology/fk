# src/fk-compute/cli.py
from __future__ import annotations
import sys
import argparse
import json
from typing import Optional, List

# Import all functionality from fk module
from .fk import (
    fk, PRESETS,
    _parse_int_list, _parse_bool,
    configure_logging
)

# Import symbolic functionality
from .symbolic_output import print_symbolic_result, SYMPY_AVAILABLE

# -------------------------------------------------------------------------
# CLI Parser
# -------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Build the main argument parser with subcommands."""
    p = argparse.ArgumentParser(
        description="""
Compute the FK invariant for braids using inversion, ILP reduction, and compiled helper binary.

The FK invariant is a mathematical object used in knot theory to distinguish different
types of knots and links. This tool provides two interfaces:

SIMPLE USAGE:
  fk simple "[1,-2,3]" 2              # Quick computation with defaults

CONFIG FILE USAGE:
  fk config single.yaml               # Single computation from file
  fk config batch.json                # Multiple computations from file

BRAID FORMATS:
  Braids can be specified in multiple formats:
  • JSON-style: "[1, -2, -3, 1]"
  • Comma-separated: "1,-2,-3,1"
  • Space-separated: "1 -2 -3 1"

For detailed documentation, run 'man fk' (after running 'fk-install-man').
""",
        prog="fk",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Simple - just the essentials
  fk simple "[1,-2,3]" 2

  # Config - single computation
  echo '{"braid": [1,-2,3], "degree": 2, "name": "my_knot", "save_data": true}' > single.json
  fk config single.json

  # Config - batch processing
  cat > batch.yaml << 'EOF'
  max_workers: 4
  save_data: true
  computations:
    - name: trefoil
      braid: [1, 1, 1]
      degree: 2
    - name: figure_eight
      braid: [1, -2, 1, -2]
      degree: 3
  EOF
  fk config batch.yaml

For more examples and detailed parameter descriptions, see 'man fk'.
"""
    )

    # Add subcommands
    subparsers = p.add_subparsers(dest="command", help="Available commands")

    # Simple compute command
    simple_parser = subparsers.add_parser(
        "simple",
        help="Simple FK computation with minimal options",
        description="Simple interface with sensible defaults. Uses quiet mode and standard parameters for quick computations.",
        epilog="Example: fk simple \"[1,-2,3]\" 2"
    )
    simple_parser.add_argument(
        "braid",
        type=str,
        help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"'
    )
    simple_parser.add_argument("degree", type=int, help="Computation degree")
    add_symbolic_simple = _parse_bool(default_true=False)
    add_symbolic_simple(simple_parser, "symbolic", "Print result in human-readable symbolic form using SymPy")

    # Config file command
    config_parser = subparsers.add_parser(
        "config",
        help="FK computation from configuration file",
        description="""
Load computation parameters from a JSON or YAML configuration file.
Supports both single computations and batch processing of multiple braids.

SINGLE COMPUTATION:
  {
    "braid": [1, -2, 3],
    "degree": 2,
    "name": "my_knot",
    "preset": "accurate",
    "max_workers": 8,
    "save_data": true,
    "ilp_file": "path/to/precomputed.ilp"
  }

BATCH PROCESSING:
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
        "braid": [1, -2, 1, -2],
        "degree": 3,
        "preset": "accurate"
      }
    ]
  }

Batch mode supports:
• Global defaults (applied to all computations)
• Individual computation overrides
• Named computations for organized results
• Progress tracking and error handling
• File naming based on computation 'name' when 'save_data' is enabled
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples: fk config single.yaml | fk config batch.json"
    )
    config_parser.add_argument(
        "config_path",
        type=str,
        help="Path to JSON or YAML configuration file"
    )

    return p


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------
def _print_result(result: dict, symbolic: bool = False) -> None:
    """Helper function to print result in requested format."""
    if symbolic:
        # Print symbolic representation if requested and SymPy is available
        if SYMPY_AVAILABLE:
            print_symbolic_result(result, format_type="pretty", show_matrix=False)
        else:
            print("Error: SymPy is required for symbolic output. Install with: pip install sympy")
            print("Falling back to JSON output:")
            print(json.dumps(result, indent=2, sort_keys=True))
    else:
        # Standard JSON output
        print(json.dumps(result, indent=2, sort_keys=True))


# -------------------------------------------------------------------------
# Command handlers
# -------------------------------------------------------------------------
def handle_simple(args) -> None:
    """Handle simple command."""
    braid = _parse_int_list(args.braid)
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

    result = fk(braid, args.degree, symbolic=getattr(args, 'symbolic', False))  # Simple mode auto-detected
    _print_result(result, getattr(args, 'symbolic', False))


def handle_config(args) -> None:
    """Handle config command."""
    result = fk(args.config_path)  # Config mode auto-detected
    _print_result(result, False)  # Config mode doesn't currently support CLI symbolic flag


# -------------------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the fk CLI tool."""
    if argv is None:
        argv = sys.argv

    parser = build_parser()
    args = parser.parse_args(argv[1:])

    # Route to appropriate handler based on subcommand
    if args.command == 'simple':
        handle_simple(args)
    elif args.command == 'config':
        handle_config(args)
    else:
        # No valid subcommand provided - show help
        parser.print_help()
