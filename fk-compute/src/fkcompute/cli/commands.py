"""
CLI commands for FK computation.

This module defines all CLI commands for the fk tool.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

import typer

from .app import app, template_app, history_app
from ..api.compute import fk
from ..infra.config import parse_int_list
from ..output.symbolic import print_symbolic_result, SYMPY_AVAILABLE


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
            from ..interactive import run_enhanced_interactive
            run_enhanced_interactive()
        elif quick:
            from ..interactive import run_quick_interactive
            run_quick_interactive()
        else:
            try:
                from ..interactive import run_enhanced_interactive
                run_enhanced_interactive()
            except ImportError:
                handle_basic_interactive()
    except ImportError as e:
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
    braid_list = parse_int_list(braid)
    if not braid_list:
        raise typer.BadParameter("Could not parse braid into a non-empty list of integers")

    result = fk(braid_list, degree, symbolic=symbolic)
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
    """
    if len(config_paths) == 1:
        result = fk(config_paths[0])
        _print_result(result, False)
        return

    typer.echo(f"\nProcessing {len(config_paths)} configuration files...\n")
    typer.echo("=" * 60)

    results = {}
    for i, config_path in enumerate(config_paths, 1):
        typer.echo(f"\n[{i}/{len(config_paths)}] Processing: {config_path}")
        typer.echo("-" * 60)

        try:
            result = fk(config_path)
            results[config_path] = result
            _print_result(result, False)
            typer.echo()

        except Exception as e:
            typer.echo(f"Error processing {config_path}: {e}", err=True)
            results[config_path] = {"error": str(e)}
            continue

    typer.echo("=" * 60)
    typer.echo(f"\nCompleted processing {len(config_paths)} configuration files.")

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
    """Create a blank YAML template configuration file."""
    output_file = Path(output_path)

    if output_file.exists() and not overwrite:
        typer.echo(f"Error: File '{output_path}' already exists. Use --overwrite to replace it.", err=True)
        raise typer.Exit(1)

    template_content = """# FK Computation Configuration Template
# This is a YAML configuration file for computing FK invariants

# Required: Braid word as a list of integers
braid: [1, -2, 1, -2]

# Required: Computation degree (positive integer)
degree: 15 

# Optional: Name for this computation (used for output files if save_data is true)
# name: my_knot

# Optional: Inversion data
# inversion: {0: [1, 1, -1, 1, 1, 1, -1, 1]}

# Optional: Preset configuration (single thread, or parallel)
# preset: single thread

# Optional: Maximum number of worker processes for parallel computation
# max_workers: 4

# Optional: Number of threads for C++ computation
# threads: 1

# Optional: Enable verbose logging
# verbose: false

# Optional: Save computation data to files
# save_data: false
"""

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
    all_items: bool = typer.Option(False, "--all", "-a", help="Show all computations")
) -> None:
    """Show computation history."""
    try:
        from ..interactive import ComputationHistory
        history = ComputationHistory()

        if all_items:
            limit = 100
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
        from ..interactive import ComputationHistory
        history = ComputationHistory()
        history.display_search_results(query)

        if typer.confirm("Select a computation from results?", default=False):
            params = history.interactive_select()
            if params:
                typer.echo(f"Selected computation: {params}")
                braid = parse_int_list(params["braid"]) if isinstance(params["braid"], str) else params["braid"]
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
        from ..interactive import ComputationHistory
        history = ComputationHistory()
        if history.clear_history():
            typer.echo("History cleared successfully")
        else:
            typer.echo("Failed to clear history")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("export")
def history_export(
    filepath: str = typer.Argument(..., help="File path to export history to")
) -> None:
    """Export computation history to a file."""
    try:
        from ..interactive import ComputationHistory
        history = ComputationHistory()
        if history.export_history(filepath):
            typer.echo(f"History exported to {filepath}")
        else:
            typer.echo("Failed to export history")
    except ImportError as e:
        typer.echo(f"History management not available: {e}")


@history_app.command("import")
def history_import(
    filepath: str = typer.Argument(..., help="File path to import history from")
) -> None:
    """Import computation history from a file."""
    try:
        from ..interactive import ComputationHistory
        history = ComputationHistory()
        if history.import_history(filepath):
            typer.echo(f"History imported from {filepath}")
        else:
            typer.echo("Failed to import history")
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

    print("Enter braid word:")
    print("  Examples: [1,-2,3], 1,-2,3, or 1 -2 3")
    braid_input = input("Braid: ").strip()
    if not braid_input:
        raise ValueError("Braid cannot be empty")

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

    print("\nEnter computation name (optional):")
    print("  (used for output files if saving)")
    name = input("Name: ").strip()
    if not name:
        name = None

    print("\nPrint result in symbolic form? (y/n)")
    symbolic_input = input("Symbolic: ").strip().lower()
    symbolic = symbolic_input in ["y", "yes", "true", "1"]

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
        if SYMPY_AVAILABLE:
            print_symbolic_result(result, format_type="pretty", show_matrix=False)
        else:
            print("Error: SymPy is required for symbolic output. Install with: pip install sympy")
            print("Falling back to JSON output:")
            print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


def handle_interactive() -> None:
    """Handle enhanced interactive mode when called directly."""
    try:
        from ..interactive import run_enhanced_interactive
        run_enhanced_interactive()
    except ImportError as e:
        typer.echo("Enhanced interactive mode not available. Using basic mode.")
        typer.echo(f"Error: {e}")
        handle_basic_interactive()


def handle_basic_interactive() -> None:
    """Handle basic interactive mode (fallback)."""
    params = _prompt_interactive()

    braid = parse_int_list(params["braid"])
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

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
