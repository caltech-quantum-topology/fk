# src/fk-compute/cli.py
from __future__ import annotations
import sys
import json
from typing import Optional, List
from pathlib import Path

# Import all functionality from fk module
from .fk import fk, PRESETS, _parse_int_list, _parse_bool, configure_logging

# Import symbolic functionality
from .symbolic_output import print_symbolic_result, SYMPY_AVAILABLE

# Import typer
try:
    import typer
except ImportError:
    # Fallback to argparse if typer is not available
    print("Error: typer is required. Please install it with: pip install typer")
    sys.exit(1)


# -------------------------------------------------------------------------
# Typer App
# -------------------------------------------------------------------------
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
  • JSON-style: "[1, -2, -3, 1]"
  • Comma-separated: "1,-2,-3,1"
  • Space-separated: "1 -2 -3 1"

For detailed documentation, run 'man fk' (after running 'fk-install-man').
""",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

template_app = typer.Typer(help="Template management commands")
history_app = typer.Typer(help="History and session management commands")

app.add_typer(template_app, name="template")
app.add_typer(history_app, name="history")


# -------------------------------------------------------------------------
# Interactive mode command
# -------------------------------------------------------------------------
@app.command("interactive", deprecated=False)
def interactive_command(
    enhanced: bool = typer.Option(False, "--enhanced", help="Use enhanced interactive wizard with progress tracking"),
    quick: bool = typer.Option(False, "--quick", help="Use quick interactive mode with minimal prompts"),
) -> None:
    """Interactive mode with guided prompts."""
    try:
        if enhanced:
            # Import and run enhanced mode
            from .interactive import run_enhanced_interactive
            run_enhanced_interactive()
        elif quick:
            # Import and run quick mode
            from .interactive import run_quick_interactive
            run_quick_interactive()
        else:
            # Try enhanced mode by default, fallback to basic if import fails
            try:
                from .interactive import run_enhanced_interactive
                run_enhanced_interactive()
            except ImportError:
                handle_basic_interactive()
    except ImportError as e:
        # Fallback to basic interactive mode if enhanced components not available
        typer.echo("Enhanced interactive mode not available. Using basic mode.")
        typer.echo(f"Error: {e}")
        handle_basic_interactive()


# -------------------------------------------------------------------------
# Simple compute command
# -------------------------------------------------------------------------
@app.command("simple")
def simple_command(
    braid: str = typer.Argument(..., help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"'),
    degree: int = typer.Argument(..., help="Computation degree"),
    symbolic: bool = typer.Option(False, "--symbolic", help="Print result in human-readable symbolic form using SymPy"),
) -> None:
    """Simple FK computation with minimal options."""
    braid_list = _parse_int_list(braid)
    if not braid_list:
        raise typer.BadParameter("Could not parse braid into a non-empty list of integers")

    result = fk(braid_list, degree, symbolic=symbolic)  # Simple mode auto-detected
    _print_result(result, symbolic)


# -------------------------------------------------------------------------
# Config file command
# -------------------------------------------------------------------------
@app.command("config")
def config_command(
    config_paths: List[str] = typer.Argument(..., help="Path(s) to JSON or YAML configuration file(s)")
) -> None:
    """FK computation from configuration file(s).

    Supports both single computation and batch processing of multiple braids.

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
    """
    # If only one config file, process normally
    if len(config_paths) == 1:
        result = fk(config_paths[0])  # Config mode auto-detected
        _print_result(result, False)
        return

    # Multiple config files - process each one
    typer.echo(f"\nProcessing {len(config_paths)} configuration files...\n")
    typer.echo("=" * 60)

    results = {}
    for i, config_path in enumerate(config_paths, 1):
        typer.echo(f"\n[{i}/{len(config_paths)}] Processing: {config_path}")
        typer.echo("-" * 60)

        try:
            result = fk(config_path)
            results[config_path] = result

            # Print result for this config
            _print_result(result, False)
            typer.echo()

        except Exception as e:
            typer.echo(f"Error processing {config_path}: {e}", err=True)
            results[config_path] = {"error": str(e)}
            continue

    typer.echo("=" * 60)
    typer.echo(f"\nCompleted processing {len(config_paths)} configuration files.")

    # Summary
    successful = sum(1 for r in results.values() if "error" not in r)
    failed = len(results) - successful
    typer.echo(f"Successful: {successful}, Failed: {failed}")


# -------------------------------------------------------------------------
# Template create subcommand
# -------------------------------------------------------------------------
@template_app.command("create")
def template_create_command(
    output_path: str = typer.Argument(
        "fk_config.yaml",
        help="Output path for template file (default: fk_config.yaml)"
    ),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing file if it exists")
) -> None:
    """Create a blank YAML template configuration file.

    This generates a commented template file with examples for single
    computation with all available options documented.
    """
    output_file = Path(output_path)

    # Check if file exists and overwrite not specified
    if output_file.exists() and not overwrite:
        typer.echo(f"Error: File '{output_path}' already exists. Use --overwrite to replace it.", err=True)
        raise typer.Exit(1)

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
    output_file.write_text(template_content)

    typer.echo(f"Template created: {output_path}")
    typer.echo(f"\nEdit the file to configure your computation, then run:")
    typer.echo(f"  fk config {output_path}")


# -------------------------------------------------------------------------
# History management commands
# -------------------------------------------------------------------------
@history_app.command("show")
def history_show(
    limit: int = typer.Option(10, "--limit", "-l", help="Number of recent computations to show"),
    all: bool = typer.Option(False, "--all", "-a", help="Show all computations")
) -> None:
    """Show computation history."""
    try:
        from .interactive import ComputationHistory
        history = ComputationHistory()
        
        if all:
            limit = 100  # Show more than default
        elif limit <= 0:
            raise typer.BadParameter("Limit must be a positive integer")
        
        history.display_recent(limit=limit)
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("search")
def history_search(
    query: str = typer.Argument(..., help="Search query for computation history")
) -> None:
    """Search computation history."""
    try:
        from .interactive import ComputationHistory
        history = ComputationHistory()
        history.display_search_results(query)
        
        # Ask if user wants to reuse a computation
        if typer.confirm("Select a computation from results?", default=False):
            params = history.interactive_select()
            if params:
                typer.echo(f"Selected computation: {params}")
                # Run computation
                braid = _parse_int_list(params["braid"]) if isinstance(params["braid"], str) else params["braid"]
                if braid is not None:
                    result = fk(braid, params["degree"], **params)
                    _print_result(result, params.get("symbolic", False))
                else:
                    typer.echo("Error: Invalid braid data")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("clear")
def history_clear(
    confirm: bool = typer.Option(False, "--confirm", "-y", help="Skip confirmation prompt")
) -> None:
    """Clear computation history."""
    if not confirm:
        if not typer.confirm("Clear all computation history? This cannot be undone.", default=False):
            typer.echo("History clearing cancelled.")
            return
    
    try:
        from .interactive import ComputationHistory
        history = ComputationHistory()
        if history.clear_history():
            typer.echo("✅ History cleared successfully")
        else:
            typer.echo("❌ Failed to clear history")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("export")
def history_export(
    filepath: str = typer.Argument(..., help="File path to export history to")
) -> None:
    """Export computation history to a file."""
    try:
        from .interactive import ComputationHistory
        history = ComputationHistory()
        if history.export_history(filepath):
            typer.echo(f"✅ History exported to {filepath}")
        else:
            typer.echo("❌ Failed to export history")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("import")
def history_import(
    filepath: str = typer.Argument(..., help="File path to import history from")
) -> None:
    """Import computation history from a file."""
    try:
        from .interactive import ComputationHistory
        history = ComputationHistory()
        if history.import_history(filepath):
            typer.echo(f"✅ History imported from {filepath}")
        else:
            typer.echo("❌ Failed to import history")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


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
# Enhanced interactive mode functions
# -------------------------------------------------------------------------
def handle_interactive() -> None:
    """Handle enhanced interactive mode when called directly."""
    try:
        # Try to import enhanced interactive components
        from .interactive import run_enhanced_interactive
        run_enhanced_interactive()
    except ImportError as e:
        # Fallback to basic interactive mode if enhanced components not available
        typer.echo("Enhanced interactive mode not available. Using basic mode.")
        typer.echo(f"Error: {e}")
        handle_basic_interactive()


def handle_basic_interactive() -> None:
    """Handle basic interactive mode (fallback)."""
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


# -------------------------------------------------------------------------
# Main entry point with legacy argument handling
# -------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the fk CLI tool."""
    if argv is None:
        argv = sys.argv

    # Handle legacy simple mode: fk "[1,-2,3]" 2 -> fk simple "[1,-2,3]" 2
    if len(argv) >= 3 and argv[1] not in [
        "simple", "config", "interactive", "template", "-h", "--help", "--version"
    ]:
        try:
            # Check if second argument is an integer (degree)
            int(argv[2])
            # If so, inject "simple" as the subcommand
            new_argv = [argv[0], "simple"] + argv[1:]
            # Update sys.argv for typer
            original_argv = sys.argv[:]
            sys.argv = new_argv
            try:
                app()
            finally:
                sys.argv = original_argv
            return
        except (ValueError, IndexError):
            # Not a simple mode pattern, proceed normally
            pass

    # Check if no arguments provided - default to interactive mode
    if len(argv) == 1:
        handle_interactive()
        return

    # Run typer normally
    app()