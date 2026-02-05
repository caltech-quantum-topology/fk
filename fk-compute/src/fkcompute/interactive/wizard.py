"""Enhanced wizard for guided FK computations with progress tracking."""

import time
from typing import Dict, Any, Optional, List
from rich.console import Console
from rich.prompt import Confirm

from ..api.compute import fk
from ..infra.config import parse_int_list
from .ui import StatusMessage, ComputationSummary, ValidatedInput
from .progress import FKProgressTracker
from .history import ComputationHistory, SessionManager
from .prompts import show_main_menu, get_computation_parameters, show_help_menu, show_settings_menu, show_search_interface

console = Console()
quit_commands = ["exit", "Exit", "quit", "Quit", "q"]


class FKWizard:
    """Enhanced wizard for guided FK computations with progress tracking."""
    
    def __init__(self):
        self.history = ComputationHistory()
        self.session = SessionManager()
        self.running = True
    
    def run_wizard(self) -> bool:
        """Run main wizard interface."""
        console.print()
        console.print("[bold blue]🧮 FK Computation Wizard[/bold blue]")
        console.print("Enhanced interactive mode with progress tracking and history management")
        
        while self.running:
            try:
                choice = show_main_menu()
                self._handle_menu_choice(choice)
            except KeyboardInterrupt:
                console.print("\n\n⚠️  Interrupted by user", style="yellow")
                if Confirm.ask("Exit wizard?", default=True):
                    self.running = False
                continue
            except Exception as e:
                StatusMessage.error(f"Unexpected error: {e}")
                if Confirm.ask("Continue?", default=True):
                    continue
                else:
                    self.running = False
        
        console.print("\n👋 Thank you for using FK Computation Wizard!", style="bold blue")
        return True
    
    def _handle_menu_choice(self, choice: str):
        """Handle main menu choice."""
        if choice == "1":
            self._handle_new_computation()
        elif choice == "2":
            self._handle_recent_computations()
        elif choice == "3":
            self._handle_search_history()
        elif choice == "4":
            self._handle_settings()
        elif choice == "5":
            self._handle_help()
        elif choice == "6":
            self.running = False
        else:
            StatusMessage.warning("Invalid choice")
    
    def _handle_new_computation(self):
        """Handle new computation workflow."""
        try:
            # Get computation parameters
            params = get_computation_parameters(self.history, self.session)
            if not params:
                return  # User cancelled
            
            # Run computation with progress tracking
            self._run_computation_with_progress(params)
            
        except Exception as e:
            StatusMessage.error(f"Computation failed: {e}")
    
    def _handle_recent_computations(self):
        """Handle recent computations display and selection."""
        console.print("\n[bold cyan]Recent Computations[/bold cyan]")
        self.history.display_recent(limit=15)
        
        if Confirm.ask("Select a computation to reuse?", default=False):
            params = self.history.interactive_select()
            if params:
                console.print("\n[bold green]Reusing computation parameters[/bold green]")
                if ComputationSummary.show(params):
                    self._run_computation_with_progress(params)
    
    def _handle_search_history(self):
        """Handle history search."""
        from .prompts import show_search_interface
        params = show_search_interface(self.history)
        if params:
            console.print("\n[bold green]Using selected computation parameters[/bold green]")
            if ComputationSummary.show(params):
                self._run_computation_with_progress(params)
    
    def _handle_settings(self):
        """Handle settings menu."""
        from .prompts import show_settings_menu
        show_settings_menu(self.session)
    
    def _handle_help(self):
        """Handle help menu."""
        from .prompts import show_help_menu
        show_help_menu()
        console.print("\nPress Enter to continue...", style="dim")
        input()
    
    def _run_computation_with_progress(self, params: Dict[str, Any]):
        """Run FK computation with enhanced progress tracking."""
        start_time = time.time()
        
        try:
            # Parse braid
            if isinstance(params['braid'], str):
                braid = parse_int_list(params['braid'])
            else:
                braid = params['braid']
            
            if not braid:
                raise ValueError("Could not parse braid into a non-empty list of integers")
            
            # Run computation without debug output

            # Extract degree before removing from params to avoid passing it twice
            # Also strip format_type which is a display concern, not a compute param
            degree = params['degree']
            clean_params = {k: v for k, v in params.items() if k not in ('degree', 'braid', 'format_type')}

            # Start progress tracking
            with FKProgressTracker() as progress:
                # Hook into computation phases
                progress.start_inversion(braid, degree)

                # Add progress callback
                clean_params['_progress_callback'] = progress

                # Run the computation
                result = fk(braid_or_config=braid, degree=degree, **clean_params)
                
                # Complete progress tracking
                progress.complete_fk_computation(len(result.get('terms', [])))
                
                if params.get('symbolic'):
                    progress.start_symbolic_generation()
                    progress.complete_symbolic_generation()
                
                if params.get('save_data'):
                    progress.start_file_operations("Saving computation data")
                    progress.complete_file_operations()
            
            # Calculate elapsed time
            computation_time = time.time() - start_time
            
            # Display results
            self._display_results(result, params, computation_time)
            
            # Save to history
            self.history.save_computation(params, result, computation_time)
            
            StatusMessage.success(f"Computation completed in {computation_time:.1f} seconds")
            
        except Exception as e:
            StatusMessage.error(f"Computation failed: {e}")
            
            # Save failed attempt to history
            error_params = params.copy()
            error_result = {'error': str(e)}
            self.history.save_computation(error_params, error_result, time.time() - start_time)
    
    def _display_results(self, result: Dict[str, Any], params: Dict[str, Any], computation_time: float):
        """Display computation results in enhanced format."""
        console.print("\n" + "─" * 60)
        console.print("[bold green]🎉 Computation Results[/bold green]")
        console.print("─" * 60)
        
        # Show summary info
        from rich.table import Table
        table = Table(show_header=False, box=None)
        table.add_column("Property", style="bold cyan")
        table.add_column("Value")
        
        table.add_row("Braid", str(params['braid']))
        table.add_row("Degree", str(params['degree']))
        table.add_row("Computation time", f"{computation_time:.1f}s")
        table.add_row("Terms found", str(len(result.get('terms', []))))
        
        if 'metadata' in result:
            metadata = result['metadata']
            if 'components' in metadata:
                table.add_row("Components", str(metadata['components']))
        
        console.print(table)
        console.print()
        
        # Show actual results using existing print function
        from ..cli.commands import _print_result
        _print_result(
            result,
            params.get('symbolic', False),
            format_type=params.get('format_type', 'pretty'),
        )


class QuickWizard:
    """Simplified wizard for rapid computations."""
    
    @staticmethod
    def quick_computation():
        """Run a quick computation with minimal prompts."""
        from .ui import ValidatedInput
        console.print("\n[bold blue]⚡ Quick FK Computation[/bold blue]")
        console.print("Minimal prompts for rapid computation\n")
        
        try:
            # Get essential parameters only
            braid = ValidatedInput.get_braid()
            degree = ValidatedInput.get_degree(braid=braid)
            symbolic_opts = ValidatedInput.get_symbolic()
            symbolic = symbolic_opts["symbolic"]
            format_type = symbolic_opts["format_type"]

            # Basic parameters
            params = {
                'braid': braid,
                'degree': degree,
                'symbolic': symbolic,
                'format_type': format_type,
                'threads': 1,
                'save_data': False
            }

            # Show quick summary
            console.print("\n[bold]Quick Summary:[/bold]")
            console.print(f"  Braid: {braid}")
            console.print(f"  Degree: {degree}")
            if symbolic:
                console.print(f"  Symbolic: Yes ({format_type})")
            else:
                console.print(f"  Symbolic: No")
            console.print()

            if not Confirm.ask("Proceed?", default=True):
                console.print("Cancelled.")
                return

            # Run computation with basic progress
            start_time = time.time()
            result = fk(braid_or_config=braid, degree=degree, symbolic=symbolic)
            computation_time = time.time() - start_time

            # Display results
            console.print("\n[bold green]✅ Results:[/bold green]")
            from ..cli.commands import _print_result
            _print_result(result, symbolic, format_type=format_type)
            
            console.print(f"\n[dim]Completed in {computation_time:.1f} seconds[/dim]")
            
        except Exception as e:
            from .ui import StatusMessage
            StatusMessage.error(f"Quick computation failed: {e}")


# Export the main wizard class for backward compatibility
def run_enhanced_interactive() -> None:
    """Run the enhanced interactive wizard."""
    wizard = FKWizard()
    wizard.run_wizard()


def run_quick_interactive() -> None:
    """Run the quick interactive mode."""
    QuickWizard.quick_computation()
