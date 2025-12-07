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
types of knots and links. This tool provides multiple interfaces for different use cases:

SIMPLE USAGE:
  fk simple "[1,-2,3]" 2              # Quick computation with defaults

PRESET USAGE:
  fk preset "[1,-2,3]" 2 --preset fast      # Fast computation
  fk preset "[1,-2,3]" 2 --preset accurate  # Thorough computation
  fk preset "[1,-2,3]" 2 --preset parallel  # Multi-core optimized

CONFIG FILE USAGE:
  fk config single.yaml               # Single computation from file
  fk config batch.json                # Multiple computations from file

ADVANCED USAGE:
  fk advanced "[1,-2,3]" 2 --max-workers 4 --verbose  # Full control

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

  # Preset - predefined configurations
  fk preset "[1,-2,3]" 2 --preset accurate

  # Config - single computation
  echo '{"braid": [1,-2,3], "degree": 2, "name": "my_knot", "preset": "fast", "save_data": true}' > single.json
  fk config single.json

  # Config - batch processing
  cat > batch.yaml << 'EOF'
  preset: fast
  max_workers: 4
  save_data: true
  computations:
    - name: trefoil
      braid: [1, 1, 1]
      degree: 2
    - name: figure_eight
      braid: [1, -2, 1, -2]
      degree: 3
      preset: accurate
  EOF
  fk config batch.yaml

  # Advanced - full parameter control
  fk advanced "[1,-2,3]" 2 --max-workers 4 --verbose --save-data

  # Legacy mode (backward compatibility)
  fk "[1,-2,3]" 2 --verbose

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

    # Preset compute command
    preset_parser = subparsers.add_parser(
        "preset",
        help="FK computation using preset configurations",
        description="""
Use predefined configuration presets optimized for different scenarios:

• fast:     Single-threaded, limited shifts, no flip symmetry (fastest)
• accurate: Multi-core, unlimited shifts, no flip symmetry, saves data (most thorough)
• parallel: High parallelism, balanced settings, no flip symmetry (best for multi-core systems)
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: fk preset \"[1,-2,3]\" 2 --preset accurate"
    )
    preset_parser.add_argument(
        "braid",
        type=str,
        help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"'
    )
    preset_parser.add_argument("degree", type=int, help="Computation degree")
    preset_parser.add_argument(
        "--preset",
        choices=list(PRESETS.keys()),
        default="fast",
        help="Preset configuration to use"
    )
    add_symbolic_preset = _parse_bool(default_true=False)
    add_symbolic_preset(preset_parser, "symbolic", "Print result in human-readable symbolic form using SymPy")

    # Advanced compute command (full options)
    advanced_parser = subparsers.add_parser(
        "advanced",
        help="Advanced FK computation with all options",
        description="Full control over all computation parameters. Use this for fine-tuned computations when presets don't meet your needs.",
        epilog="Example: fk advanced \"[1,-2,3]\" 2 --max-workers 4 --verbose --save-data"
    )
    _add_advanced_arguments(advanced_parser)

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

    # Variables command - print symbolic relations
    variables_parser = subparsers.add_parser(
        "variables",
        help="Print symbolic relations for a braid",
        description="""
Print the symbolic relations and constraints for FK computation at a given degree.
This shows the mathematical structure of the problem in human-readable form,
including variable assignments, degree constraints, and relation inequalities.

For fibered braids, inversion data will be computed automatically.
For homogeneous braids, no additional data is needed.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example: fk variables "[1,-2,3]" 2'
    )
    variables_parser.add_argument(
        "braid",
        type=str,
        help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"'
    )
    variables_parser.add_argument("degree", type=int, help="Computation degree")
    variables_parser.add_argument(
        "-f", "--outfile",
        type=str,
        default=None,
        help="Output file to save symbolic relations (optional)"
    )

    # Empty config command - show parameters and generate template
    template_parser = subparsers.add_parser(
        "template",
        help="Show configuration parameters or generate template config file",
        description="""
Display all available configuration parameters and their meanings, or generate
a blank configuration file template that you can fill in.

MODES:
  fk template --list                 # List all parameters with descriptions
  fk template --generate config.yaml # Generate blank template file
  fk template --generate config.json # Generate JSON template file
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples: fk template --list | fk template --generate my_config.yaml"
    )
    template_group = template_parser.add_mutually_exclusive_group(required=True)
    template_group.add_argument(
        "--list",
        action="store_true",
        help="List all available configuration parameters with descriptions"
    )
    template_group.add_argument(
        "--generate",
        type=str,
        metavar="FILE",
        help="Generate a template configuration file (supports .yaml, .yml, .json)"
    )

    # Legacy compute command (backward compatibility)
    legacy_parser = subparsers.add_parser(
        "compute",
        help="Legacy compute command (same as advanced)",
        description="Legacy interface for backward compatibility. Equivalent to the 'advanced' subcommand.",
        epilog="Example: fk compute \"[1,-2,3]\" 2 --max-workers 4"
    )
    _add_advanced_arguments(legacy_parser)

    # If no subcommand provided, default to legacy behavior
    p.set_defaults(command="legacy")

    return p


def _add_advanced_arguments(parser: argparse.ArgumentParser) -> None:
    """Add all advanced arguments to a parser."""
    parser.add_argument(
        "braid",
        type=str,
        help='Braid word. Examples: "[1,-2,3]", "1,-2,3", or "1 -2 3"'
    )
    parser.add_argument("degree", type=int, help="Degree of the invariant computation")

    parser.add_argument("--ilp-file", type=str, default=None, help="Path to a precomputed ILP file")
    parser.add_argument("--ilp", type=str, default=None, help="ILP data as a raw string (advanced)")

    parser.add_argument("--inversion-file", type=str, default=None, help="Path to inversion JSON file")
    parser.add_argument(
        "--partial-signs",
        type=str,
        default=None,
        help='Optional partial sign assignments. Same formats as braid'
    )

    parser.add_argument("--max-workers", type=int, default=1, help="Parallel workers for inversion calculation")
    parser.add_argument("--chunk-size", type=int, default=(1 << 14), help="Chunk size for parallel tasks")
    parser.add_argument("--max-shifts", type=int, default=None, help="Maximum shifts considered in inversion")

    # Toggle flags
    add_verbose = _parse_bool(default_true=False)
    add_verbose(parser, "verbose", "Verbose logging")

    add_flip = _parse_bool(default_true=False)
    add_flip(parser, "include-flip", "Include flip symmetry in inversion")

    add_save = _parse_bool(default_true=False)
    add_save(parser, "save-data", "Save intermediate files (inversion/ILP/JSON)")

    add_print = _parse_bool(default_true=True)
    add_print(parser, "print-result", "Print the FK result to stdout")

    add_symbolic = _parse_bool(default_true=False)
    add_symbolic(parser, "symbolic", "Print result in human-readable symbolic form using SymPy")

    parser.add_argument("--save-dir", type=str, default="data", help="Directory to store data when --save-data is set")
    parser.add_argument("--link-name", type=str, default=None, help="Base name for output files (without extension)")


def build_legacy_parser() -> argparse.ArgumentParser:
    """Build legacy parser for backward compatibility when no subcommands used."""
    p = argparse.ArgumentParser(
        description="""
LEGACY MODE - Compute FK invariant from a braid using inversion and ILP stages.

This is the legacy interface maintained for backward compatibility.
For new usage, consider using the subcommand interface:

  fk simple "[1,-2,3]" 2                    # Quick computation
  fk preset "[1,-2,3]" 2 --preset accurate  # Preset-based
  fk advanced "[1,-2,3]" 2 --verbose        # Full control
  fk config myconfig.yaml                    # From config file

Run 'fk -h' to see all available subcommands and improved help.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: fk \"[1,-2,3]\" 2 --max-workers 4 --verbose"
    )
    _add_advanced_arguments(p)
    return p


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------
def _looks_like_negative_braid(s: str) -> bool:
    """Check if string looks like a braid with negative numbers."""
    import re
    # Pattern: starts with -, followed by digits/commas/more negative numbers
    # Examples: -1, -1,2, -1,-2,-3, -2,3,-1
    pattern = r'^-\d+(,-?\d+)*$'
    return bool(re.match(pattern, s))


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


def handle_preset(args) -> None:
    """Handle preset command."""
    braid = _parse_int_list(args.braid)
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

    result = fk(braid, args.degree, preset=args.preset, symbolic=getattr(args, 'symbolic', False))  # Preset mode
    _print_result(result, getattr(args, 'symbolic', False))


def handle_config(args) -> None:
    """Handle config command."""
    result = fk(args.config_path)  # Config mode auto-detected
    _print_result(result, False)  # Config mode doesn't currently support CLI symbolic flag


def handle_variables(args) -> None:
    """Handle variables command."""
    from .braidstates_links import BraidStates
    from .relations_links import full_reduce, print_symbolic_relations
    from .braids import is_homogeneous_braid

    braid = _parse_int_list(args.braid)
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

    # Determine if braid is homogeneous or fibered
    is_homogeneous = is_homogeneous_braid(braid)

    if is_homogeneous:
        # Homogeneous braid - no inversion data needed
        braid_states = BraidStates(braid)
        all_relations = braid_states.get_state_relations()
        relations = full_reduce(all_relations)
        print_symbolic_relations(args.degree, relations, braid_states, args.outfile)
    else:
        # Fibered braid - compute inversion data automatically
        print("Computing inversion data for fibered braid...")

        # Use the FK function to compute inversion data
        result = fk(braid, args.degree, verbose=False, save_data=False)

        if 'inversion_data' not in result:
            raise ValueError("Could not compute inversion data for this braid")

        # Extract inversion data and set up braid states
        inversion_data = result['inversion_data']
        braid_states = BraidStates(braid)
        braid_states.strand_signs = inversion_data
        braid_states.compute_matrices()

        if braid_states.validate():
            braid_states.generate_position_assignments()
            all_relations = braid_states.get_state_relations()
            relations = full_reduce(all_relations)
            print_symbolic_relations(args.degree, relations, braid_states, args.outfile)
        else:
            raise ValueError("Invalid inversion data computed for this braid")


def handle_template(args) -> None:
    """Handle template command."""
    if args.list:
        _print_config_parameters()
    elif args.generate:
        _generate_config_template(args.generate)


def _print_config_parameters() -> None:
    """Print all available configuration parameters with descriptions."""
    params = _get_config_parameter_docs()

    print("=" * 80)
    print("FK COMPUTATION CONFIGURATION PARAMETERS")
    print("=" * 80)
    print()

    for category, items in params.items():
        print(f"\n{category}")
        print("-" * len(category))
        for param_name, param_info in items.items():
            print(f"\n  {param_name}:")
            print(f"    Type:        {param_info['type']}")
            if 'default' in param_info:
                print(f"    Default:     {param_info['default']}")
            if 'options' in param_info:
                print(f"    Options:     {', '.join(str(o) for o in param_info['options'])}")
            print(f"    Description: {param_info['description']}")

    print("\n" + "=" * 80)
    print("USAGE EXAMPLES:")
    print("=" * 80)
    print("""
Minimal configuration (JSON):
  {
    "braid": [1, -2, 3],
    "degree": 2
  }

With optional parameters (YAML):
  braid: [1, -2, 3]
  degree: 2
  name: "my_knot"
  max_workers: 4
  save_data: true
  verbose: true
""")


def _generate_config_template(filepath: str) -> None:
    """Generate a blank configuration file template."""
    is_yaml = filepath.endswith(('.yaml', '.yml'))
    is_json = filepath.endswith('.json')

    if not is_yaml and not is_json:
        print(f"Error: File must have .yaml, .yml, or .json extension")
        sys.exit(1)

    if is_yaml:
        try:
            import yaml
        except ImportError:
            print("Error: PyYAML is required for YAML files. Install with: pip install PyYAML")
            sys.exit(1)

    template = _get_config_template(is_yaml)

    try:
        with open(filepath, 'w') as f:
            f.write(template)
        print(f"✓ Template configuration file created: {filepath}")
        print(f"\nEdit this file and run:")
        print(f"  fk config {filepath}")
    except Exception as e:
        print(f"Error writing template file: {e}")
        sys.exit(1)


def _get_config_parameter_docs() -> dict:
    """Return structured documentation for all configuration parameters."""
    return {
        "REQUIRED PARAMETERS": {
            "braid": {
                "type": "list[int]",
                "description": "Braid word as list of integers (e.g., [1, -2, 3])"
            },
            "degree": {
                "type": "int",
                "description": "Computation degree for FK invariant"
            }
        },
        "PRE-COMPUTED DATA": {
            "ilp": {
                "type": "string",
                "default": "None",
                "description": "Pre-computed ILP data as raw string (advanced)"
            },
            "ilp_file": {
                "type": "string",
                "default": "None",
                "description": "Path to pre-computed ILP file"
            },
            "inversion": {
                "type": "dict",
                "default": "None",
                "description": "Pre-computed inversion data as dict mapping component indices to sign lists"
            },
            "inversion_file": {
                "type": "string",
                "default": "None",
                "description": "Path to pre-computed inversion JSON file"
            }
        },
        "COMPUTATION PARAMETERS": {
            "partial_signs": {
                "type": "list[int]",
                "default": "None",
                "description": "Optional partial sign assignments for inversion calculation"
            },
            "max_workers": {
                "type": "int",
                "default": "1",
                "description": "Number of parallel workers for inversion calculation"
            },
            "chunk_size": {
                "type": "int",
                "default": "16384",
                "description": "Chunk size for parallel task processing"
            },
            "include_flip": {
                "type": "bool",
                "default": "false",
                "description": "Include flip symmetry in inversion computation"
            },
            "max_shifts": {
                "type": "int",
                "default": "None (unlimited)",
                "description": "Maximum number of cyclic shifts to consider"
            },
            "threads": {
                "type": "int",
                "default": "1",
                "description": "Number of threads for C++ FK computation engine"
            }
        },
        "OUTPUT AND LOGGING": {
            "verbose": {
                "type": "bool",
                "default": "false",
                "description": "Enable verbose logging output"
            },
            "save_data": {
                "type": "bool",
                "default": "false",
                "description": "Save intermediate files (inversion, ILP, results)"
            },
            "save_dir": {
                "type": "string",
                "default": "data",
                "description": "Directory path for saved data files"
            },
            "link_name": {
                "type": "string",
                "default": "auto-generated",
                "description": "Base name for output files (without extension)"
            },
            "symbolic": {
                "type": "bool",
                "default": "false",
                "description": "Generate symbolic polynomial representation using SymPy"
            },
            "name": {
                "type": "string",
                "default": "None",
                "description": "Name for this computation (used for file naming when save_data is enabled)"
            }
        }
    }


def _get_config_template(is_yaml: bool) -> str:
    """Generate a template configuration file content."""
    if is_yaml:
        return """# FK Computation Configuration File (YAML)
# ==========================================
# This is a template configuration file for FK invariant computation.
# Uncomment and modify the parameters you need.

# =============================================================================
# REQUIRED PARAMETERS
# =============================================================================
braid: [1, -2, 3]  # Braid word as list of integers
degree: 2           # Computation degree

# =============================================================================
# OPTIONAL PARAMETERS
# =============================================================================

# Name for this computation (used for file naming when save_data is true)
# name: "my_knot"

# Pre-computed data (provide file path or data directly)
# ilp_file: "path/to/precomputed.ilp"
# inversion_file: "path/to/inversion.json"
# inversion:
#   0: [1, -1, 1]
#   1: [-1, 1]

# Computation parameters
# max_workers: 1        # Number of parallel workers for inversion
# chunk_size: 16384     # Chunk size for parallel processing
# include_flip: false   # Include flip symmetry
# max_shifts: null      # Maximum cyclic shifts (null = unlimited)
# partial_signs: []     # Partial sign assignments
# threads: 1            # Number of threads for C++ computation

# Output and logging
# verbose: false        # Enable verbose logging
# save_data: false      # Save intermediate files
# save_dir: "data"      # Directory for saved files
# link_name: null       # Base name for output files
# symbolic: false       # Generate symbolic representation (requires SymPy)
"""
    else:  # JSON
        return """{
  "_comment": "FK Computation Configuration File (JSON)",
  "_comment2": "Remove all _comment fields before use",

  "braid": [1, -2, 3],
  "degree": 2,

  "_optional": "Optional parameters below - uncomment as needed",
  "name": "my_knot",

  "_precomputed": "Optional pre-computed data",
  "ilp_file": null,
  "inversion_file": null,
  "inversion": null,

  "_computation": "Computation parameters",
  "max_workers": 1,
  "chunk_size": 16384,
  "include_flip": false,
  "max_shifts": null,
  "partial_signs": null,
  "threads": 1,

  "_output": "Output and logging",
  "verbose": false,
  "save_data": false,
  "save_dir": "data",
  "link_name": null,
  "symbolic": false
}
"""


def handle_advanced(args) -> None:
    """Handle advanced/compute/legacy commands."""
    braid = _parse_int_list(args.braid)
    if not braid:
        raise ValueError("Could not parse braid into a non-empty list of integers")

    partial_signs = _parse_int_list(getattr(args, 'partial_signs', None))

    # Advanced mode - pass all parameters to unified fk function
    result = fk(
        braid,
        args.degree,
        ilp=getattr(args, 'ilp', None),
        ilp_file=getattr(args, 'ilp_file', None),
        inversion_file=getattr(args, 'inversion_file', None),
        partial_signs=partial_signs,
        max_workers=getattr(args, 'max_workers', 1),
        chunk_size=getattr(args, 'chunk_size', 1 << 14),
        include_flip=getattr(args, 'include_flip', False),
        max_shifts=getattr(args, 'max_shifts', None),
        verbose=getattr(args, 'verbose', False),
        save_data=getattr(args, 'save_data', False),
        save_dir=getattr(args, 'save_dir', 'data'),
        link_name=getattr(args, 'link_name', None),
        symbolic=getattr(args, 'symbolic', False),
    )

    # Print result in requested format
    if getattr(args, 'print_result', True):
        _print_result(result, getattr(args, 'symbolic', False))


# -------------------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv

    # Check if using help flag or subcommands
    has_help = len(argv) > 1 and argv[1] in ['-h', '--help']
    has_subcommand = len(argv) > 1 and argv[1] in ['simple', 'preset', 'config', 'advanced', 'compute', 'variables', 'template']

    if has_help or has_subcommand:
        parser = build_parser()
        args = parser.parse_args(argv[1:])

        if args.command == 'simple':
            handle_simple(args)
        elif args.command == 'preset':
            handle_preset(args)
        elif args.command == 'config':
            handle_config(args)
        elif args.command == 'variables':
            handle_variables(args)
        elif args.command == 'template':
            handle_template(args)
        elif args.command in ['advanced', 'compute']:
            handle_advanced(args)
        else:
            parser.print_help()
    else:
        # Legacy mode - parse as before for backward compatibility
        parser = build_legacy_parser()
        try:
            args = parser.parse_args(argv[1:])
            handle_advanced(args)
        except SystemExit as e:
            # Check if this might be due to negative numbers being treated as flags
            # Look for patterns like -1, -1,2, -1,-2,-3 etc.
            if len(argv) >= 3 and argv[1].startswith('-'):
                # Check if it looks like a braid string with negative numbers
                braid_str = argv[1]
                if _looks_like_negative_braid(braid_str):
                    # Try again with -- separator to treat negative numbers as positional args
                    try:
                        args = parser.parse_args(['--'] + argv[1:])
                        handle_advanced(args)
                    except SystemExit:
                        # If it still fails, show helpful error message
                        print(f"Error: Could not parse braid '{argv[1]}'. Try using quotes:", file=sys.stderr)
                        print(f"  fk \"{argv[1]}\" {' '.join(argv[2:])}", file=sys.stderr)
                        print("Or use the new subcommand interface:", file=sys.stderr)
                        print(f"  fk simple \"{argv[1]}\" {' '.join(argv[2:])}", file=sys.stderr)
                        raise
                else:
                    raise
            else:
                raise
