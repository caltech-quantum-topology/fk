# src/fk-compute/cli.py
from __future__ import annotations
import sys
import argparse
import json
from typing import Optional, List

# Import all functionality from fk module
from .fk import fk, PRESETS, _parse_int_list, _parse_bool, configure_logging

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
  • JSON-style: "[1, -2, -3, 1]"
  • Comma-separated: "1,-2,-3,1"
  • Space-separated: "1 -2 -3 1"

For detailed documentation, run 'man fk' (after running 'fk-install-man').
""",
        prog="fk",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Interactive - guided prompts
  fk
  fk interactive

  # Simple - just the essentials
  fk simple "[1,-2,3]" 2

  # Template - create configuration file
  fk template create my_config.yaml
  # Edit my_config.yaml, then run:
  fk config my_config.yaml

  # Config - single computation
  echo '{"braid": [1,-2,3], "degree": 2, "name": "my_knot", "save_data": true}' > single.json
  fk config single.json

  # Config - multiple files
  fk config config1.yaml config2.yaml config3.json

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
""",
    )

    # Add subcommands
    subparsers = p.add_subparsers(dest="command", help="Available commands")

    # Interactive mode command
    interactive_parser = subparsers.add_parser(
        "interactive",
        help="Interactive mode with guided prompts",
        description="Interactive mode that prompts you for all required parameters step by step.",
        epilog="Example: fk interactive",
    )

    # Simple compute command
    simple_parser = subparsers.add_parser(
        "simple",
        help="Simple FK computation with minimal options",
        description="Simple interface with sensible defaults. Uses quiet mode and standard parameters for quick computations.",
        epilog='Example: fk simple "[1,-2,3]" 2',
    )
    simple_parser.add_argument(
        "braid",
        type=str,
        help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"',
    )
    simple_parser.add_argument("degree", type=int, help="Computation degree")
    add_symbolic_simple = _parse_bool(default_true=False)
    add_symbolic_simple(
        simple_parser,
        "symbolic",
        "Print result in human-readable symbolic form using SymPy",
    )

    # Config file command
    config_parser = subparsers.add_parser(
        "config",
        help="FK computation from configuration file(s)",
        description="""
Load computation parameters from JSON or YAML configuration file(s).
Supports both single computations and batch processing of multiple braids.

MULTIPLE CONFIG FILES:
  Process multiple configuration files in sequence:
  fk config config1.yaml config2.yaml config3.json

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
        epilog="Examples: fk config single.yaml | fk config batch.json | fk config config1.yaml config2.yaml",
    )
    config_parser.add_argument(
        "config_paths", type=str, nargs="+", help="Path(s) to JSON or YAML configuration file(s)"
    )

    # Template command
    template_parser = subparsers.add_parser(
        "template",
        help="Template management commands",
        description="Create and manage FK computation configuration templates."
    )
    template_subparsers = template_parser.add_subparsers(dest="template_command", help="Template commands")

    # Template create subcommand
    template_create_parser = template_subparsers.add_parser(
        "create",
        help="Create a blank YAML template configuration file",
        description="""
Create a blank YAML template for FK computation configuration.

This generates a commented template file with examples for single
computation with all available options documented.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: fk template create my_config.yaml",
    )
    template_create_parser.add_argument(
        "output_path",
        type=str,
        nargs="?",
        default="fk_config.yaml",
        help="Output path for template file (default: fk_config.yaml)"
    )
    template_create_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing file if it exists"
    )

    return p


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------
def _prompt_interactive() -> dict:
    """Prompt user for parameters in interactive mode."""
    print("\n" + "=" * 60)
    print("FK Computation - Interactive Mode")
    print("=" * 60 + "\n")

    # Prompt for braid
    print("Enter braid word:")
    print("  Examples: [1,-2,3], 1,-2,3, or 1 -2 3")
    braid_input = input("Braid: ").strip()
    if not braid_input:
        raise ValueError("Braid cannot be empty")

    # Prompt for degree
    print("\nEnter computation degree:")
    print("  (positive integer)")
    while True:
        degree_input = input("Degree: ").strip()
        try:
            degree = int(degree_input)
            if degree > 0:
                break
            print("  Error: Degree must be a positive integer")
        except ValueError:
            print("  Error: Please enter a valid integer")

    # Prompt for threads (optional)
    print("\nEnter number of threads (optional):")
    print("  (press Enter for default)")
    threads_input = input("Threads: ").strip()
    threads = None
    if threads_input:
        try:
            threads = int(threads_input)
            if threads <= 0:
                print("  Warning: Invalid thread count, using default")
                threads = None
        except ValueError:
            print("  Warning: Invalid thread count, using default")
            threads = None

    # Prompt for name (optional)
    print("\nEnter computation name (optional):")
    print("  (used for output files if saving)")
    name = input("Name: ").strip()
    if not name:
        name = None

    # Prompt for symbolic output
    print("\nPrint result in symbolic form? (y/n)")
    symbolic_input = input("Symbolic: ").strip().lower()
    symbolic = symbolic_input in ["y", "yes", "true", "1"]

    # Prompt for save
    print("\nSave the output?")
    save_input = input("Save: ").strip().lower()
    save = save_input in ["y", "yes", "true", "1"]

    print("\n" + "=" * 60)
    print("Running computation...")
    print("=" * 60 + "\n")

    return {
        "braid": braid_input,
        "degree": degree,
        "threads": threads,
        "name": name,
        "symbolic": symbolic,
        "save_data": save
    }


def _print_result(result: dict, symbolic: bool = False) -> None:
    """Helper function to print result in requested format."""
    if symbolic:
        # Print symbolic representation if requested and SymPy is available
        if SYMPY_AVAILABLE:
            print_symbolic_result(result, format_type="pretty", show_matrix=False)
        else:
            print(
                "Error: SymPy is required for symbolic output. Install with: pip install sympy"
            )
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

    result = fk(
        braid, args.degree, symbolic=getattr(args, "symbolic", False)
    )  # Simple mode auto-detected
    _print_result(result, getattr(args, "symbolic", False))


def handle_config(args) -> None:
    """Handle config command."""
    config_paths = args.config_paths

    # If only one config file, process normally
    if len(config_paths) == 1:
        result = fk(config_paths[0])  # Config mode auto-detected
        _print_result(result, False)
        return

    # Multiple config files - process each one
    print(f"\nProcessing {len(config_paths)} configuration files...\n")
    print("=" * 60)

    results = {}
    for i, config_path in enumerate(config_paths, 1):
        print(f"\n[{i}/{len(config_paths)}] Processing: {config_path}")
        print("-" * 60)

        try:
            result = fk(config_path)
            results[config_path] = result

            # Print result for this config
            _print_result(result, False)
            print()

        except Exception as e:
            print(f"Error processing {config_path}: {e}")
            results[config_path] = {"error": str(e)}
            continue

    print("=" * 60)
    print(f"\nCompleted processing {len(config_paths)} configuration files.")

    # Summary
    successful = sum(1 for r in results.values() if "error" not in r)
    failed = len(results) - successful
    print(f"Successful: {successful}, Failed: {failed}")


def handle_interactive(args) -> None:
    """Handle interactive command."""
    params = _prompt_interactive()

    # Parse braid
    braid = _parse_int_list(params["braid"])
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

    # Call fk with collected parameters
    # Build kwargs based on what user provided
    kwargs = {
        "symbolic": params["symbolic"],
        "save_data": params["save_data"]
    }
    if params["threads"] is not None:
        kwargs["threads"] = params["threads"]
    if params["name"] is not None:
        kwargs["name"] = params["name"]

    result = fk(braid, params["degree"], **kwargs)
    _print_result(result, params["symbolic"])


def handle_template_create(args) -> None:
    """Handle template create command."""
    import os

    output_path = args.output_path

    # Check if file exists and overwrite not specified
    if os.path.exists(output_path) and not args.overwrite:
        print(f"Error: File '{output_path}' already exists. Use --overwrite to replace it.")
        sys.exit(1)

    # Generate template content
    template_content = """# FK Computation Configuration Template
# This is a YAML configuration file for computing FK invariants

# Required: Braid word as a list of integers
# Examples: [1, -2, 3], [1, 1, 1] (trefoil), [1, -2, 1, -2] (figure eight)
braid: [1, -2, 3]

# Required: Computation degree (positive integer)
degree: 2

# Optional: Name for this computation (used for output files if save_data is true)
# name: my_knot

# Optional: Preset configuration (fast, accurate, or parallel)
# Presets provide sensible defaults for different use cases
# preset: accurate

# Optional: Maximum number of worker processes for parallel computation
# max_workers: 4

# Optional: Chunk size for parallel processing
# chunk_size: 65536

# Optional: Number of threads for C++ computation
# threads: 1

# Optional: Include flip symmetry in computation
# include_flip: false

# Optional: Maximum number of shifts to consider
# max_shifts: null

# Optional: Enable verbose logging
# verbose: false

# Optional: Save computation data to files
# save_data: false

# Optional: Path to precomputed ILP file (advanced usage)
# ilp_file: path/to/precomputed.ilp

# Optional: Precomputed inversion data (advanced usage)
# Provide component->signs mapping
# inversion:
#   0: [1, -1, 1]
#   1: [-1, 1]
"""

    # Write template to file
    with open(output_path, 'w') as f:
        f.write(template_content)

    print(f"Template created: {output_path}")
    print(f"\nEdit the file to configure your computation, then run:")
    print(f"  fk config {output_path}")


# -------------------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the fk CLI tool."""
    if argv is None:
        argv = sys.argv

    parser = build_parser()

    # Check if no arguments provided - default to interactive mode
    if len(argv) == 1:
        argv = argv[:1] + ["interactive"]

    # Check if user is passing braid and degree without explicit "simple" subcommand
    # e.g., fk "[1,-2,3]" 2 instead of fk simple "[1,-2,3]" 2
    if len(argv) >= 3 and argv[1] not in [
        "simple",
        "config",
        "interactive",
        "-h",
        "--help",
    ]:
        # Try to detect if first arg looks like a braid and second like a degree
        try:
            # Check if second argument is an integer (degree)
            int(argv[2])
            # If so, insert "simple" as the subcommand
            argv = argv[:1] + ["simple"] + argv[1:]
        except (ValueError, IndexError):
            # Not a simple mode pattern, proceed normally
            pass

    args = parser.parse_args(argv[1:])

    # Route to appropriate handler based on subcommand
    if args.command == "simple":
        handle_simple(args)
    elif args.command == "config":
        handle_config(args)
    elif args.command == "interactive":
        handle_interactive(args)
    elif args.command == "template":
        if args.template_command == "create":
            handle_template_create(args)
        else:
            # Show template subcommand help
            parser.parse_args(["template", "-h"])
    else:
        # No valid subcommand provided - show help
        parser.print_help()
