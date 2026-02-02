"""Progress tracking for FK computations."""

from typing import Optional, Dict, Any, Union
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn
from rich.live import Live
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
import time
import threading

console = Console()


class FKProgressTracker:
    """Tracks and displays progress during FK computation phases."""
    
    def __init__(self):
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False
        )
        self.tasks = {}
        self.current_phase = None
        self.start_time = None
        self.live = None
        self.lock = threading.Lock()
    
    def __enter__(self):
        """Start progress tracking context."""
        self.start_time = time.time()
        self.live = Live(self.progress, console=console, refresh_per_second=10)
        self.live.start()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop progress tracking."""
        if self.live:
            self.live.stop()
    
    def start_inversion(self, braid: list, degree: int) -> Optional[Any]:
        """Start inversion calculation phase."""
        with self.lock:
            task_id = self.progress.add_task(
                f"🔍 Computing inversion (degree {degree})", 
                total=None
            )
            self.tasks['inversion'] = task_id
            self.current_phase = "inversion"
            
            StatusMessage.info(f"Starting inversion computation for {len(braid)}-crossing braid")
            return task_id
    
    def update_inversion_progress(self, current: int, total: int, description: str = ""):
        """Update inversion progress."""
        with self.lock:
            if 'inversion' in self.tasks:
                if description:
                    self.progress.update(
                        self.tasks['inversion'], 
                        description=f"🔍 {description}",
                        completed=current,
                        total=total
                    )
                else:
                    self.progress.update(
                        self.tasks['inversion'],
                        completed=current,
                        total=total
                    )
    
    def complete_inversion(self, components: int):
        """Mark inversion as complete."""
        with self.lock:
            if 'inversion' in self.tasks:
                self.progress.update(
                    self.tasks['inversion'],
                    description=f"✅ Inversion complete ({components} components)",
                    completed=100,
                    total=100
                )
                time.sleep(0.5)  # Brief pause to show completion
    
    def start_ilp_generation(self) -> Optional[Any]:
        """Start ILP generation phase."""
        with self.lock:
            task_id = self.progress.add_task(
                "📋 Generating ILP formulation",
                total=None
            )
            self.tasks['ilp'] = task_id
            self.current_phase = "ilp"
            StatusMessage.info("Generating Integer Linear Programming formulation")
            return task_id
    
    def complete_ilp_generation(self, constraints: int):
        """Mark ILP generation as complete."""
        with self.lock:
            if 'ilp' in self.tasks:
                self.progress.update(
                    self.tasks['ilp'],
                    description=f"✅ ILP complete ({constraints} constraints)",
                    completed=100,
                    total=100
                )
                time.sleep(0.5)
    
    def start_fk_computation(self, threads: int, total_points: int) -> Optional[Any]:
        """Start FK computation phase."""
        with self.lock:
            task_id = self.progress.add_task(
                f"⚡ Computing FK invariant ({threads} threads)",
                total=total_points
            )
            self.tasks['fk'] = task_id
            self.current_phase = "fk"
            StatusMessage.info(f"Computing FK invariant for {total_points} points")
            return task_id
    
    def update_fk_progress(self, current: int, total: int, point_info: Optional[Dict] = None):
        """Update FK computation progress."""
        with self.lock:
            if 'fk' in self.tasks:
                if point_info:
                    description = f"⚡ Computing point {current}/{total} (deg {point_info.get('degree', '?')})"
                else:
                    description = f"⚡ Computing FK invariant ({current}/{total})"
                
                self.progress.update(
                    self.tasks['fk'],
                    description=description,
                    completed=current,
                    total=total
                )
    
    def complete_fk_computation(self, terms_count: int):
        """Mark FK computation as complete."""
        with self.lock:
            if 'fk' in self.tasks:
                self.progress.update(
                    self.tasks['fk'],
                    description=f"✅ FK complete ({terms_count} terms)",
                    completed=100,
                    total=100
                )
                time.sleep(0.5)
    
    def start_symbolic_generation(self) -> Optional[Any]:
        """Start symbolic output generation."""
        with self.lock:
            task_id = self.progress.add_task(
                "🎨 Generating symbolic representation",
                total=None
            )
            self.tasks['symbolic'] = task_id
            self.current_phase = "symbolic"
            StatusMessage.info("Converting to symbolic polynomial format")
            return task_id
    
    def complete_symbolic_generation(self):
        """Mark symbolic generation as complete."""
        with self.lock:
            if 'symbolic' in self.tasks:
                self.progress.update(
                    self.tasks['symbolic'],
                    description="✅ Symbolic generation complete",
                    completed=100,
                    total=100
                )
                time.sleep(0.5)
    
    def start_file_operations(self, operation: str) -> Optional[Any]:
        """Start file operations (saving/loading)."""
        with self.lock:
            task_id = self.progress.add_task(
                f"💾 {operation}",
                total=None
            )
            self.tasks['file_ops'] = task_id
            return task_id
    
    def complete_file_operations(self):
        """Mark file operations as complete."""
        with self.lock:
            if 'file_ops' in self.tasks:
                self.progress.update(
                    self.tasks['file_ops'],
                    description="✅ File operations complete",
                    completed=100,
                    total=100
                )
                time.sleep(0.3)
    
    def show_phase_summary(self):
        """Display a summary of current computation phase."""
        if not self.start_time:
            return
        
        elapsed = time.time() - self.start_time
        
        # Create summary panel
        summary_text = Text()
        summary_text.append("Computation Status\n", style="bold blue")
        summary_text.append(f"Elapsed time: {elapsed:.1f}s\n", style="dim")
        
        if self.current_phase:
            phase_descriptions = {
                'inversion': 'Calculating sign assignments',
                'ilp': 'Generating ILP constraints',
                'fk': 'Computing FK invariant',
                'symbolic': 'Creating symbolic output',
                'file_ops': 'Saving/loading files'
            }
            
            current_desc = phase_descriptions.get(self.current_phase, 'Processing')
            summary_text.append(f"Current phase: {current_desc}", style="cyan")
        
        panel = Panel(
            summary_text,
            title="Status",
            border_style="blue",
            padding=(0, 2)
        )
        
        console.print(panel)


class StatusMessage:
    """Display status messages during computation."""
    
    @staticmethod
    def info(message: str):
        console.print(f"ℹ️  {message}", style="bold blue")
    
    @staticmethod
    def success(message: str):
        console.print(f"✅ {message}", style="bold green")
    
    @staticmethod
    def warning(message: str):
        console.print(f"⚠️  {message}", style="bold yellow")
    
    @staticmethod
    def error(message: str):
        console.print(f"❌ {message}", style="bold red")


class ResourceMonitor:
    """Monitor and display system resources during computation."""
    
    def __init__(self):
        self.monitoring = False
        self.start_time = None
    
    def start_monitoring(self):
        """Start resource monitoring."""
        self.start_time = time.time()
        self.monitoring = True
    
    def stop_monitoring(self):
        """Stop resource monitoring."""
        self.monitoring = False
    
    def get_status(self) -> Dict[str, Any]:
        """Get current resource status."""
        if not self.start_time:
            return {}
        
        elapsed = time.time() - self.start_time
        
        # Basic status (in a real implementation, you might use psutil for detailed monitoring)
        status = {
            'elapsed': elapsed,
            'memory_usage': 'N/A',  # Would require psutil
            'cpu_usage': 'N/A'     # Would require psutil
        }
        
        return status
    
    def display_status(self):
        """Display current resource status."""
        if not self.monitoring:
            return
        
        status = self.get_status()
        
        # Create status table
        table = Table(show_header=False, box=None, width=40)
        table.add_column("Metric", style="dim")
        table.add_column("Value", style="bold")
        
        table.add_row("Elapsed", f"{status['elapsed']:.1f}s")
        table.add_row("Memory", status['memory_usage'])
        table.add_row("CPU", status['cpu_usage'])
        
        panel = Panel(
            table,
            title="Resources",
            border_style="yellow",
            padding=(0, 1)
        )
        
        console.print(panel)