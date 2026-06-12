"""High-level prompts and menu functions for interactive mode."""

from typing import Dict, Any, Optional
from rich.console import Console
from rich.prompt import Prompt, Confirm
from rich.panel import Panel

from .ui import ValidatedInput, ComputationSummary

console = Console()


def show_main_menu() -> str:
    """Display main interactive menu."""
    console.print()

    menu_text = """
[bold cyan]FK Computation - Interactive Mode[/bold cyan]

[bold]Choose an option:[/bold]

1. 🧮 [cyan]New Computation[/cyan] - Start a new FK computation
2. ❓ [cyan]Help & Examples[/cyan] - Learn about FK invariants
3. 🚪 [cyan]Exit[/cyan] - Leave interactive mode
    """

    panel = Panel(
        menu_text.strip(), title="Main Menu", border_style="blue", padding=(1, 2)
    )

    console.print(panel)

    choice = Prompt.ask(
        "Select option", choices=["1", "2", "3"], default="1"
    )

    return choice


def get_computation_parameters() -> Optional[Dict[str, Any]]:
    """Collect all computation parameters interactively."""
    params = {}

    # Step 1: Get basic inputs
    console.print("\n[bold blue]Step 1: Basic Parameters[/bold blue]")
    console.print("Let's start with the essential information for your computation.\n")

    braid = ValidatedInput.get_braid()
    params["braid"] = braid

    degree = ValidatedInput.get_degree(braid=braid)
    params["degree"] = degree

    # Step 2: Preset selection
    console.print("\n[bold blue]Step 2: Computation Preset[/bold blue]")
    console.print("Choose a preset for optimized settings, or configure manually.\n")

    preset = ValidatedInput.get_preset()
    if preset:
        params["preset"] = preset
        params.update(_get_preset_parameters(preset))
    else:
        custom_params = ValidatedInput.get_custom_parameters()
        params.update(custom_params)

    # Step 3: Thread configuration
    console.print("\n[bold blue]Step 3: Thread Configuration[/bold blue]")
    console.print("Configure parallel computation for your system.\n")

    if not preset:
        threads = ValidatedInput.get_threads(preset)
    else:
        from ..api.presets import PRESETS
        threads = PRESETS[preset].get("threads", 1)

    params["threads"] = threads

    # Step 4: Optional settings
    console.print("\n[bold blue]Step 4: Optional Settings[/bold blue]")
    console.print("Configure output format and data saving.\n")

    name = ValidatedInput.get_name()
    if name:
        params["name"] = name

    symbolic_opts = ValidatedInput.get_symbolic()
    symbolic = symbolic_opts["symbolic"]
    params["symbolic"] = symbolic
    params["format_type"] = symbolic_opts["format_type"]

    save_data = ValidatedInput.get_save_preference()
    params["save_data"] = save_data

    # Step 5: Review and confirm
    console.print("\n[bold blue]Step 5: Review & Confirm[/bold blue]")
    console.print("Please review your computation settings before starting.\n")

    if ComputationSummary.show(params):
        return params
    else:
        console.print("Computation cancelled.", style="yellow")
        return None


def show_help_menu():
    """Display help and examples."""
    help_text = """
[bold cyan]FK Invariant Computation - Help[/bold cyan]

[bold]What is the FK Invariant?[/bold]
The FK invariant is a mathematical invariant of knots and links 
used in knot theory. It was introduced by Gukov and Manolescu in
arXiv:1904.060597. Here we compute it using the large color R-matrix formalism 
which was introduced by Park in arXiv:2004.02087. This program implements 
a search for inversion data, as well as a sum over R-matrix indices

[bold]Basic Concepts:[/bold]
• [green]Braid[/green]: A sequence of integers representing crossings
• [green]Degree[/green]: The degree in the x-variables
• [green]Components[/green]: The number of components in the link (one for knots)
• [green]Homogeneous[/green]: Braids with trivial sign diagrams
• [green]Fibered[/green]: Braids requiring additional inversion data

[bold]Common Examples:[/bold]
• Trefoil knot: [1, 1, 1]
• Figure-eight knot: [1, -2, 1, -2]
• Hopf link: [1, -1]

[bold]Preset Configurations:[/bold]
• [cyan]Fast[/cyan]: Quick computation, minimal resources (degree 1, 1 worker)
• [cyan]Accurate[/cyan]: Thorough computation, saves data (degree 2, 4 workers)
• [cyan]Parallel[/cyan]: High-performance (degree 3, 8 workers, 4 threads)

[bold]Tips:[/bold]
• Start with low degrees for testing
• Use homogeneous braids for faster computation
• Increase threads for parallel computations
• Always check the results against the q->1 limit and the MMR expansion

[bold]For more information:[/bold]
• Run 'fk --help' for command-line options
• Read the manual with 'man fk' (after 'fk-install-man')
• Visit the project documentation online
    """

    panel = Panel(
        help_text.strip(), title="Help & Examples", border_style="blue", padding=(1, 2)
    )

    console.print(panel)


def _get_preset_parameters(preset_name: str) -> Dict[str, Any]:
    """Get parameters for a specific preset."""
    from ..api.presets import PRESETS

    preset = PRESETS.get(preset_name, {})
    return {
        "max_workers": preset.get("max_workers", 1),
        "chunk_size": preset.get("chunk_size", 65536),
        "include_flip": preset.get("include_flip", False),
        "max_shifts": preset.get("max_shifts", None),
        "verbose": preset.get("verbose", False),
        "save_data": preset.get("save_data", False),
        # Note: intentionally NOT including 'degree' to avoid overriding user's choice
    }
