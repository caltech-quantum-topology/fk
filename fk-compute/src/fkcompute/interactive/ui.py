"""Rich UI components for enhanced interactive experience."""

from typing import Any, List, Optional, Dict, Union
import json
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt, Confirm
from rich.text import Text
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.live import Live
from rich.layout import Layout
from rich.align import Align
from rich.rule import Rule

console = Console()


class BorderedSection:
    """Creates a bordered section with consistent styling."""
    
    def __init__(self, title: str, color: str = "blue", border_style: str = "blue"):
        self.title = title
        self.color = color
        self.border_style = border_style
        self.console = Console()
    
    def __enter__(self):
        # Display header with title
        title_text = Text(f" {self.title} ", style=f"bold {self.color}")
        console.print(Rule(title_text, style=self.border_style))
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # Display footer
        console.print(Rule(style=self.border_style))
        console.print()


class StatusMessage:
    """Displays colored status messages with appropriate icons."""
    
    @staticmethod
    def success(message: str):
        """Display success message with green checkmark."""
        console.print(f"✅ {message}", style="bold green")
    
    @staticmethod
    def error(message: str):
        """Display error message with red X."""
        console.print(f"❌ {message}", style="bold red")
    
    @staticmethod
    def warning(message: str):
        """Display warning message with yellow warning sign."""
        console.print(f"⚠️  {message}", style="bold yellow")
    
    @staticmethod
    def info(message: str):
        """Display info message with blue info icon."""
        console.print(f"ℹ️  {message}", style="bold blue")
    
    @staticmethod
    def step(step_num: int, total: int, message: str):
        """Display step indicator."""
        console.print(f"[{step_num}/{total}] {message}", style="cyan")


class ValidatedInput:
    """Enhanced input handling with real-time validation."""
    
    @staticmethod
    def get_braid(prompt: str = "Enter braid word") -> List[int]:
        """Get and validate braid input with examples and real-time feedback."""
        while True:
            with BorderedSection("Braid Input", color="cyan"):
                # Show examples
                console.print("Examples:", style="dim")
                console.print("  • [1, -2, 3]  - JSON-style list")
                console.print("  • 1,-2,3      - Comma-separated")
                console.print("  • 1 -2 3      - Space-separated")
                console.print()
                
                # Common braid examples
                console.print("Common braids:", style="dim")
                console.print("  • [1,1,1]    - Trefoil knot (3 crossings)")
                console.print("  • [1,-2,1,-2] - Figure-eight knot (4 crossings)")
                console.print()
                
                braid_input = console.input(f"[bold]{prompt}:[/] ").strip()
                
                if not braid_input:
                    StatusMessage.error("Braid cannot be empty")
                    continue
                
                # Parse and validate braid
                try:
                    braid = ValidatedInput._parse_braid(braid_input)
                    
                    # Show braid analysis
                    ValidatedInput._analyze_braid(braid)
                    
                    # Confirm if user wants to proceed
                    if Confirm.ask("Continue with this braid?"):
                        return braid
                        
                except ValueError as e:
                    StatusMessage.error(f"Invalid braid: {e}")
    
    @staticmethod
    def get_degree(prompt: str = "Enter computation degree", braid: Optional[List[int]] = None) -> int:
        """Get and validate degree with smart suggestions."""
        while True:
            with BorderedSection("Degree Selection", color="green"):
                if braid:
                    # Provide degree suggestions based on braid complexity
                    crossings = len(braid)
                    suggestions = []
                    if crossings <= 3:
                        suggestions = [(1, "Quick test"), (2, "Recommended")]
                    elif crossings <= 6:
                        suggestions = [(1, "Basic"), (2, "Recommended"), (3, "Thorough")]
                    else:
                        suggestions = [(1, "Basic"), (2, "Moderate"), (3, "Thorough"), (4, "Comprehensive")]
                    
                    console.print(f"Based on {crossings} crossings, suggested degrees:", style="dim")
                    for deg, desc in suggestions:
                        console.print(f"  • {deg} - {desc}")
                    console.print()
                
                degree_input = console.input(f"[bold]{prompt}:[/] ").strip()
                
                try:
                    degree = int(degree_input)
                    if degree <= 0:
                        StatusMessage.error("Degree must be a positive integer")
                        continue
                    
                    # Warn about high degrees
                    if degree > 4:
                        StatusMessage.warning(f"Degree {degree} may take significant time to compute")
                        if not Confirm.ask("Continue anyway?"):
                            continue
                    
                    return degree
                    
                except ValueError:
                    StatusMessage.error("Please enter a valid integer")
    
    @staticmethod
    def get_preset() -> Optional[str]:
        """Let user select from available presets."""
        from ..fk import PRESETS
        
        presets = list(PRESETS.keys())
        
        with BorderedSection("Preset Selection", color="magenta"):
            console.print("Choose a computation preset for optimized settings:", style="dim")
            console.print()
            
            # Create table for preset comparison
            table = Table(show_header=True, header_style="bold blue")
            table.add_column("Preset", style="cyan", width=10)
            table.add_column("Workers", justify="center", width=8)
            table.add_column("Threads", justify="center", width=8)
            table.add_column("Description", style="dim")
            
            for preset_name in presets:
                preset = PRESETS[preset_name]
                table.add_row(
                    preset_name,
                    str(preset.get('max_workers', 'N/A')),
                    str(preset.get('threads', 'N/A')),
                    {
                        'fast': 'Quick computation, minimal resources',
                        'accurate': 'Thorough computation, saves data',
                        'parallel': 'High-performance parallel processing'
                    }.get(preset_name, 'Custom settings')
                )
            
            console.print(table)
            console.print()
            
            # Add custom option
            console.print("0. [dim]Custom settings[/dim] - Configure all parameters manually")
            console.print()
            
            choice = Prompt.ask(
                "Select preset",
                choices=presets + ["0", "custom"],
                default="accurate"
            )
            
            if choice in ["0", "custom"]:
                return None
            return choice
    
    @staticmethod
    def get_custom_parameters() -> Dict[str, Any]:
        """Get custom computation parameters interactively."""
        params = {}
        
        with BorderedSection("Custom Parameters", color="yellow"):
            console.print("Configure computation parameters manually:", style="dim")
            console.print()
            
            # Max workers
            params['max_workers'] = int(Prompt.ask(
                "Max workers",
                default=4,
                show_default=True
            ))
            
            # Chunk size
            params['chunk_size'] = int(Prompt.ask(
                "Chunk size",
                default=65536,
                show_default=True
            ))
            
            # Include flip
            params['include_flip'] = Confirm.ask(
                "Include flip symmetry",
                default=False
            )
            
            # Max shifts
            if Confirm.ask("Limit shifts"):
                params['max_shifts'] = int(Prompt.ask("Maximum shifts"))
            else:
                params['max_shifts'] = None
            
            # Verbose
            params['verbose'] = Confirm.ask(
                "Enable verbose logging",
                default=True
            )
            
            # Save data
            params['save_data'] = Confirm.ask(
                "Save computation data",
                default=True
            )
            
            return params
    
    @staticmethod
    def get_threads(preset: Optional[str] = None) -> int:
        """Get thread count with preset-based suggestions."""
        with BorderedSection("Thread Configuration", color="green"):
            preset_threads = 1
            if preset:
                from ..fk import PRESETS
                preset_threads = PRESETS[preset].get('threads', 1)
                console.print(f"Preset '{preset}' suggests {preset_threads} threads", style="dim")
                console.print()
            
            threads = int(Prompt.ask(
                "Number of threads",
                default=preset_threads,
                show_default=True
            ))
            
            if threads > 8:
                StatusMessage.warning(f"{threads} threads may not provide additional benefit due to hardware limitations")
            
            return threads
    
    @staticmethod
    def get_name() -> Optional[str]:
        """Get optional computation name."""
        with BorderedSection("Computation Naming", color="blue"):
            console.print("Name is optional and used for saved files", style="dim")
            console.print()
            
            name = Prompt.ask("Computation name (press Enter to skip)", default="").strip()
            return name if name else None
    
    @staticmethod
    def get_symbolic() -> bool:
        """Get symbolic output preference."""
        with BorderedSection("Output Format", color="purple"):
            console.print("Choose output format:", style="dim")
            console.print("  • JSON - Machine-readable, standard format")
            console.print("  • Symbolic - Human-readable polynomial format (requires SymPy)")
            console.print()
            
            return Confirm.ask(
                "Generate symbolic output",
                default=False
            )
    
    @staticmethod
    def get_save_preference() -> bool:
        """Get save data preference."""
        return Confirm.ask(
            "Save computation data to files",
            default=False
        )
    
    @staticmethod
    def _parse_braid(braid_input: str) -> List[int]:
        """Parse braid string into list of integers."""
        # Remove surrounding brackets if present
        braid_input = braid_input.strip()
        if braid_input.startswith('[') and braid_input.endswith(']'):
            braid_input = braid_input[1:-1]
        
        # Try JSON first
        try:
            val = json.loads(f"[{braid_input}]")
            if isinstance(val, list):
                return [int(x) for x in val]
        except Exception:
            pass
        
        # Try splitting on commas and spaces
        parts = []
        for part in braid_input.replace(',', ' ').split():
            if part.strip():
                parts.append(int(part.strip()))
        
        if not parts:
            raise ValueError("Could not parse braid - no numbers found")
        
        return parts
    
    @staticmethod
    def _analyze_braid(braid: List[int]):
        """Analyze and display information about the braid."""
        try:
            from ..braids import is_homogeneous_braid
            from ..braidstates_links import BraidStates
            
            # Basic statistics
            crossings = len(braid)
            max_strand = max(abs(x) for x in braid)
            min_strand = min(abs(x) for x in braid)
            components = max_strand
            
            console.print(Rule("Braid Analysis", style="cyan"))
            
            # Show statistics
            stats_table = Table(show_header=False, box=None)
            stats_table.add_column("Property", style="bold")
            stats_table.add_column("Value")
            
            stats_table.add_row("Crossings", str(crossings))
            stats_table.add_row("Strands", str(components))
            stats_table.add_row("Max absolute value", str(max_strand))
            
            # Type classification
            if is_homogeneous_braid(braid):
                braid_type = "Homogeneous (standard computation)"
                type_style = "green"
            else:
                braid_type = "Fibered (requires inversion data)"
                type_style = "yellow"
            
            stats_table.add_row("Type", f"[{type_style}]{braid_type}[/{type_style}]")
            
            console.print(stats_table)
            
            # Complexity estimate
            if crossings <= 3:
                complexity = "Simple"
                color = "green"
            elif crossings <= 6:
                complexity = "Moderate" 
                color = "yellow"
            else:
                complexity = "Complex"
                color = "red"
            
            console.print(f"\nComplexity: [{color}]{complexity}[/{color}]", style="bold")
            
        except Exception as e:
            console.print(f"Analysis failed: {e}", style="dim red")


class ComputationSummary:
    """Display a summary of computation parameters before execution."""
    
    @staticmethod
    def show(params: Dict[str, Any]) -> bool:
        """Show parameter summary and ask for confirmation."""
        with BorderedSection("Computation Summary", color="blue"):
            # Create summary table
            table = Table(show_header=False, box=None)
            table.add_column("Parameter", style="bold cyan", width=15)
            table.add_column("Value")
            
            # Braid (show abbreviated if long)
            braid = params['braid']
            if len(braid) > 8:
                braid_str = f"[{', '.join(map(str, braid[:3]))}, ..., {', '.join(map(str, braid[-3:]))}]"
            else:
                braid_str = str(braid)
            table.add_row("Braid", braid_str)
            table.add_row("Degree", str(params['degree']))
            
            if params.get('preset'):
                table.add_row("Preset", params['preset'])
            
            if params.get('threads'):
                table.add_row("Threads", str(params['threads']))
            
            if params.get('name'):
                table.add_row("Name", params['name'])
            
            table.add_row("Symbolic", "Yes" if params.get('symbolic') else "No")
            table.add_row("Save data", "Yes" if params.get('save_data') else "No")
            
            console.print(table)
            console.print()
            
            return Confirm.ask("Proceed with computation?", default=True)