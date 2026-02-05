"""High-level prompts and menu functions for interactive mode."""

from typing import Dict, Any, Optional
from rich.console import Console
from rich.prompt import Prompt, Confirm
from rich.table import Table
from rich.panel import Panel


from .ui import ValidatedInput, ComputationSummary
from .history import ComputationHistory, SessionManager

console = Console()


def show_main_menu() -> str:
    """Display main interactive menu."""
    console.print()

    # Create nice menu panel
    menu_text = """
[bold cyan]FK Computation - Interactive Mode[/bold cyan]

[bold]Choose an option:[/bold]

1. 🧮 [cyan]New Computation[/cyan] - Start a new FK computation
2. 📋 [cyan]Recent Computations[/cyan] - View and reuse previous results
3. 🔍 [cyan]Search History[/cyan] - Search computation history
4. ⚙️  [cyan]Settings[/cyan] - Configure preferences
5. ❓ [cyan]Help & Examples[/cyan] - Learn about FK invariants
6. 🚪 [cyan]Exit[/cyan] - Leave interactive mode
    """

    panel = Panel(
        menu_text.strip(), title="Main Menu", border_style="blue", padding=(1, 2)
    )

    console.print(panel)

    choice = Prompt.ask(
        "Select option", choices=["1", "2", "3", "4", "5", "6"], default="1"
    )

    return choice


def get_computation_parameters(
    history: ComputationHistory, session: SessionManager
) -> Optional[Dict[str, Any]]:
    """Collect all computation parameters interactively."""
    params = {}

    # Step 1: Get basic inputs
    console.print("\n[bold blue]Step 1: Basic Parameters[/bold blue]")
    console.print("Let's start with the essential information for your computation.\n")

    # Get braid with enhanced validation
    braid = ValidatedInput.get_braid()
    params["braid"] = braid

    # Get degree with smart suggestions
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
        # Custom parameters
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
    console.print("Configure naming, output format, and data saving.\n")

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
        # Save to session preferences
        session.save_preference("last_preset", preset if preset else "custom")
        session.save_preference("last_threads", threads)
        session.save_preference("last_symbolic", symbolic)
        session.save_preference("last_save_data", save_data)

        return params
    else:
        # User cancelled
        console.print("Computation cancelled.", style="yellow")
        return None


def show_help_menu():
    """Display help and examples."""
    help_text = """
[bold cyan]FK Invariant Computation - Help[/bold cyan]

[bold]What is the FK Invariant?[/bold]
The FK invariant is a mathematical object used in knot theory to distinguish 
different types of knots and links. It's computed using inversion, 
Integer Linear Programming, and compiled helper binaries.

[bold]Basic Concepts:[/bold]
• [green]Braid[/green]: A sequence of integers representing crossings
• [green]Degree[/green]: The computational degree (higher = more detailed)
• [green]Components[/green]: Number of strands in the braid
• [green]Homogeneous[/green]: Braids that can be computed directly
• [green]Fibered[/green]: Braids requiring additional inversion data

[bold]Common Examples:[/bold]
• Trefoil knot: [1, 1, 1]
• Figure-eight knot: [1, -2, 1, -2] 
• Hopf link: [1, -1]
• Borromean rings: [1, 2, 1, 2, 1, 2]

[bold]Preset Configurations:[/bold]
• [cyan]Fast[/cyan]: Quick computation, minimal resources (degree 1, 1 worker)
• [cyan]Accurate[/cyan]: Thorough computation, saves data (degree 2, 4 workers)  
• [cyan]Parallel[/cyan]: High-performance (degree 3, 8 workers, 4 threads)

[bold]Tips:[/bold]
• Start with low degrees (1-2) for testing
• Use homogeneous braids for faster computation
• Increase threads for parallel computations
• Save computations you want to reference later
• Use the history search to find similar computations

[bold]For more information:[/bold]
• Run 'fk --help' for command-line options
• Read the manual with 'man fk' (after 'fk-install-man')
• Visit the project documentation online
    """

    panel = Panel(
        help_text.strip(), title="Help & Examples", border_style="blue", padding=(1, 2)
    )

    console.print(panel)


def show_settings_menu(session: SessionManager) -> bool:
    """Display and configure user settings."""
    while True:
        console.print()

        # Current settings table
        table = Table(show_header=False, box=None, width=50)
        table.add_column("Setting", style="bold")
        table.add_column("Value", style="cyan")

        table.add_row(
            "Default preset", session.get_preference("default_preset", "single thread")
        )
        table.add_row(
            "Default threads", str(session.get_preference("default_threads", 1))
        )
        table.add_row(
            "Auto-save history", str(session.get_preference("auto_save", False))
        )
        table.add_row(
            "Show advanced", str(session.get_preference("show_advanced", False))
        )
        table.add_row("Theme", session.get_preference("theme", "default"))

        console.print(
            Panel(table, title="Current Settings", border_style="blue", padding=(1, 2))
        )

        # Settings menu
        settings_menu = """
[bold]Settings Options:[/bold]

1. 🎯 Default preset - Change default computation preset
2. 🧵 Default threads - Change default thread count
3. 💾 Auto-save history - Toggle automatic history saving
4. 🎨 Theme - Change visual theme
5. 🗑️  Clear history - Remove all saved computations
6. ⬅️  Back - Return to main menu
        """

        console.print(
            Panel(
                settings_menu.strip(),
                title="Settings Menu",
                border_style="green",
                padding=(1, 2),
            )
        )

        choice = Prompt.ask(
            "Select setting", choices=["1", "2", "3", "4", "5", "6"], default="6"
        )

        if choice == "1":
            _change_default_preset(session)
        elif choice == "2":
            _change_default_threads(session)
        elif choice == "3":
            _toggle_auto_save(session)
        elif choice == "4":
            _change_theme(session)
        elif choice == "5":
            _clear_history_menu()
        elif choice == "6":
            return True  # Return to main menu


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


def _change_default_preset(session: SessionManager):
    """Change the default preset setting."""
    presets = ["single thread", "parallel"]
    current = session.get_preference("default_preset", "single thread")

    console.print(f"\nCurrent default preset: [cyan]{current}[/cyan]")

    new_preset = Prompt.ask("Select default preset", choices=presets, default=current)

    session.save_preference("default_preset", new_preset)
    console.print(f"✅ Default preset changed to [cyan]{new_preset}[/cyan]")


def _change_default_threads(session: SessionManager):
    """Change the default thread setting."""
    current = session.get_preference("default_threads", 1)

    console.print(f"\nCurrent default threads: [cyan]{current}[/cyan]")

    new_threads = int(Prompt.ask("Enter default thread count", default=str(current)))

    session.save_preference("default_threads", new_threads)
    console.print(f"✅ Default threads changed to [cyan]{new_threads}[/cyan]")


def _toggle_auto_save(session: SessionManager):
    """Toggle auto-save history setting."""
    current = session.get_preference("auto_save", False)

    console.print(f"\nCurrent auto-save: [cyan]{'On' if current else 'Off'}[/cyan]")

    new_value = Confirm.ask("Enable auto-save history?", default=current)

    session.save_preference("auto_save", new_value)
    console.print(
        f"✅ Auto-save changed to [cyan]{'On' if new_value else 'Off'}[/cyan]"
    )


def _change_theme(session: SessionManager):
    """Change the visual theme."""
    current = session.get_preference("theme", "default")
    themes = ["default", "dark", "light", "colorful"]

    console.print(f"\nCurrent theme: [cyan]{current}[/cyan]")

    new_theme = Prompt.ask("Select theme", choices=themes, default=current)

    session.save_preference("theme", new_theme)
    console.print(f"✅ Theme changed to [cyan]{new_theme}[/cyan]")


def _clear_history_menu():
    """Handle history clearing confirmation."""
    if Confirm.ask(
        "Are you sure you want to clear all computation history?\n"
        "This action cannot be undone.",
        default=False,
    ):
        history = ComputationHistory()
        history.clear_history()
        console.print("✅ History cleared successfully")
    else:
        console.print("History clearing cancelled")


def show_search_interface(history: ComputationHistory) -> Optional[Dict[str, Any]]:
    """Interactive search interface for computation history."""
    console.print("\n[bold cyan]Search Computation History[/bold cyan]")
    console.print("Search by braid pattern, name, or other criteria.\n")

    query = Prompt.ask("Enter search terms", default="").strip()

    if not query:
        console.print("No search terms provided.", style="yellow")
        return None

    # Show search results
    history.display_search_results(query)

    # Ask if user wants to select a result
    if Confirm.ask("Select a computation from results?", default=False):
        return history.interactive_select()

    return None
