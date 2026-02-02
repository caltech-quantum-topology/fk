"""Computation history and session management."""

import json
import os
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt, Confirm
from rich.panel import Panel
from rich.text import Text

console = Console()


class ComputationHistory:
    """Manages computation history and sessions."""
    
    def __init__(self, history_file: Optional[str] = None):
        self.history_file = history_file or self._get_default_history_file()
        self.history_dir = Path(self.history_file).parent
        self.history_dir.mkdir(parents=True, exist_ok=True)
    
    def _get_default_history_file(self) -> str:
        """Get default history file path."""
        home = Path.home()
        return str(home / ".fkcompute" / "history.json")
    
    def save_computation(self, params: Dict[str, Any], result: Dict[str, Any], 
                      computation_time: float = 0.0) -> str:
        """Save a computation to history."""
        # Ensure history file exists
        if not os.path.exists(self.history_file):
            self._create_empty_history()
        
        # Load existing history
        history = self._load_history()
        
        # Create computation entry
        computation_id = self._generate_computation_id()
        entry = {
            'id': computation_id,
            'timestamp': datetime.now().isoformat(),
            'params': params,
            'result_summary': self._summarize_result(result),
            'computation_time': computation_time,
            'success': 'error' not in result
        }
        
        # Add to history
        history['computations'].append(entry)
        
        # Limit history size (keep last 100 computations)
        if len(history['computations']) > 100:
            history['computations'] = history['computations'][-100:]
        
        # Save updated history
        self._save_history(history)
        
        return computation_id
    
    def get_recent_computations(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent computations."""
        history = self._load_history()
        return history['computations'][-limit:]
    
    def search_computations(self, query: str) -> List[Dict[str, Any]]:
        """Search computations by braid, name, or other criteria."""
        history = self._load_history()
        results = []
        
        query = query.lower()
        
        for comp in history['computations']:
            # Search in various fields
            searchable_text = " ".join([
                comp.get('params', {}).get('name', ''),
                str(comp.get('params', {}).get('braid', [])),
                comp.get('params', {}).get('preset', ''),
                comp.get('result_summary', {}).get('type', '')
            ]).lower()
            
            if query in searchable_text:
                results.append(comp)
        
        return results
    
    def load_configuration(self, computation_id: str) -> Optional[Dict[str, Any]]:
        """Load parameters from a previous computation."""
        history = self._load_history()
        
        for comp in history['computations']:
            if comp['id'] == computation_id:
                return comp['params'].copy()
        
        return None
    
    def delete_computation(self, computation_id: str) -> bool:
        """Delete a computation from history."""
        history = self._load_history()
        
        original_length = len(history['computations'])
        history['computations'] = [
            comp for comp in history['computations'] 
            if comp['id'] != computation_id
        ]
        
        if len(history['computations']) < original_length:
            self._save_history(history)
            return True
        
        return False
    
    def clear_history(self) -> bool:
        """Clear all computation history."""
        if Confirm.ask("Clear all computation history? This cannot be undone.", default=False):
            self._create_empty_history()
            return True
        return False
    
    def export_history(self, filepath: str) -> bool:
        """Export history to a file."""
        try:
            history = self._load_history()
            with open(filepath, 'w') as f:
                json.dump(history, f, indent=2)
            return True
        except Exception as e:
            console.print(f"Failed to export history: {e}", style="red")
            return False
    
    def import_history(self, filepath: str) -> bool:
        """Import history from a file."""
        try:
            with open(filepath, 'r') as f:
                imported_history = json.load(f)
            
            # Validate structure
            if 'computations' not in imported_history:
                console.print("Invalid history file format", style="red")
                return False
            
            # Merge with existing history
            current_history = self._load_history()
            current_history['computations'].extend(imported_history['computations'])
            
            # Limit size and save
            if len(current_history['computations']) > 100:
                current_history['computations'] = current_history['computations'][-100:]
            
            self._save_history(current_history)
            return True
            
        except Exception as e:
            console.print(f"Failed to import history: {e}", style="red")
            return False
    
    def display_recent(self, limit: int = 10):
        """Display recent computations in a nice table."""
        computations = self.get_recent_computations(limit)
        
        if not computations:
            console.print("No computation history found", style="dim")
            return
        
        console.print(f"\n[bold]Recent Computations (Last {len(computations)})[/bold]")
        
        table = Table(show_header=True, header_style="bold blue")
        table.add_column("ID", width=8)
        table.add_column("Time", width=16)
        table.add_column("Braid", width=20)
        table.add_column("Degree", width=6)
        table.add_column("Type", width=10)
        table.add_column("Time", width=8)
        table.add_column("Status", width=8)
        
        for comp in reversed(computations):
            # Format timestamp
            try:
                dt = datetime.fromisoformat(comp['timestamp'])
                time_str = dt.strftime("%m/%d %H:%M")
            except:
                time_str = "Unknown"
            
            # Format braid (abbreviate if long)
            braid = comp.get('params', {}).get('braid', [])
            if len(braid) > 6:
                braid_str = f"[{', '.join(map(str, braid[:3]))},..., {braid[-1]}]"
            else:
                braid_str = str(braid)
            
            # Get other info
            degree = comp.get('params', {}).get('degree', '?')
            comp_type = comp.get('result_summary', {}).get('type', 'Unknown')
            comp_time = comp.get('computation_time', 0)
            status = "✅ Success" if comp.get('success') else "❌ Failed"
            status_style = "green" if comp.get('success') else "red"
            
            table.add_row(
                comp['id'][:8],
                time_str,
                braid_str,
                str(degree),
                comp_type,
                f"{comp_time:.1f}s",
                f"[{status_style}]{status}[/{status_style}]"
            )
        
        console.print(table)
    
    def display_search_results(self, query: str):
        """Display search results."""
        results = self.search_computations(query)
        
        if not results:
            console.print(f"No results found for: {query}", style="dim")
            return
        
        console.print(f"\n[bold]Search Results for: '{query}'[/bold]")
        
        table = Table(show_header=True, header_style="bold blue")
        table.add_column("ID", width=8)
        table.add_column("Time", width=16)
        table.add_column("Braid", width=20)
        table.add_column("Degree", width=6)
        table.add_column("Name", width=15)
        table.add_column("Status", width=8)
        
        for comp in results:
            # Format timestamp
            try:
                dt = datetime.fromisoformat(comp['timestamp'])
                time_str = dt.strftime("%m/%d %H:%M")
            except:
                time_str = "Unknown"
            
            # Get info
            braid = comp.get('params', {}).get('braid', [])
            degree = comp.get('params', {}).get('degree', '?')
            name = comp.get('params', {}).get('name', '')[:15]
            status = "✅ Success" if comp.get('success') else "❌ Failed"
            status_style = "green" if comp.get('success') else "red"
            
            # Format braid
            if len(braid) > 6:
                braid_str = f"[{', '.join(map(str, braid[:3]))},..., {braid[-1]}]"
            else:
                braid_str = str(braid)
            
            table.add_row(
                comp['id'][:8],
                time_str,
                braid_str,
                str(degree),
                name,
                f"[{status_style}]{status}[/{status_style}]"
            )
        
        console.print(table)
    
    def interactive_select(self) -> Optional[Dict[str, Any]]:
        """Let user interactively select a previous computation."""
        computations = self.get_recent_computations(20)
        
        if not computations:
            console.print("No computation history found", style="dim")
            return None
        
        console.print("\n[bold]Select a computation to reuse:[/bold]")
        
        # Create numbered list
        for i, comp in enumerate(reversed(computations), 1):
            braid = comp.get('params', {}).get('braid', [])
            degree = comp.get('params', {}).get('degree', '?')
            name = comp.get('params', {}).get('name', '')
            time_str = comp.get('timestamp', '')[:16]  # Just date time part
            
            display_name = f"{name} - " if name else ""
            console.print(f"{i:2d}. {display_name}Braid {braid}, Degree {degree} ({time_str})")
        
        console.print(" 0. Cancel")
        
        choice = Prompt.ask("Select computation", choices=[str(i) for i in range(len(computations) + 1)])
        
        if choice == "0":
            return None
        
        try:
            index = int(choice) - 1
            if 0 <= index < len(computations):
                selected_comp = reversed(computations).__next__() if index == 0 else list(reversed(computations))[index]
                return selected_comp['params'].copy()
        except (ValueError, IndexError):
            pass
        
        return None
    
    def _load_history(self) -> Dict[str, Any]:
        """Load history from file."""
        try:
            with open(self.history_file, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return self._create_empty_history()
    
    def _save_history(self, history: Dict[str, Any]):
        """Save history to file."""
        with open(self.history_file, 'w') as f:
            json.dump(history, f, indent=2)
    
    def _create_empty_history(self) -> Dict[str, Any]:
        """Create empty history file."""
        empty_history = {
            'version': '1.0',
            'created': datetime.now().isoformat(),
            'computations': []
        }
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(self.history_file), exist_ok=True)
        
        self._save_history(empty_history)
        return empty_history
    
    def _generate_computation_id(self) -> str:
        """Generate unique computation ID."""
        import uuid
        return str(uuid.uuid4())[:8]
    
    def _summarize_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Create a summary of the computation result."""
        summary = {}
        
        if 'error' in result:
            summary['type'] = 'Error'
            summary['error'] = result['error']
        else:
            # Extract key information from result
            if 'metadata' in result:
                metadata = result['metadata']
                summary['components'] = metadata.get('components', 'Unknown')
                summary['terms_count'] = len(result.get('terms', []))
                summary['type'] = f"{summary['components']} components, {summary['terms_count']} terms"
            else:
                summary['type'] = 'Unknown format'
        
        return summary


class SessionManager:
    """Manages user session and preferences."""
    
    def __init__(self, config_file: Optional[str] = None):
        self.config_file = config_file or self._get_default_config_file()
        self.preferences = self._load_preferences()
    
    def _get_default_config_file(self) -> str:
        """Get default config file path."""
        home = Path.home()
        return str(home / ".fkcompute" / "config.json")
    
    def save_preference(self, key: str, value: Any):
        """Save a user preference."""
        self.preferences[key] = value
        self._save_preferences()
    
    def get_preference(self, key: str, default: Any = None) -> Any:
        """Get a user preference."""
        return self.preferences.get(key, default)
    
    def _load_preferences(self) -> Dict[str, Any]:
        """Load preferences from file."""
        try:
            with open(self.config_file, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            # Create default preferences
            defaults = {
                'default_preset': 'accurate',
                'default_threads': 1,
                'auto_save': False,
                'theme': 'default',
                'recent_presets': [],
                'show_advanced': False
            }
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.config_file), exist_ok=True)
            
            with open(self.config_file, 'w') as f:
                json.dump(defaults, f, indent=2)
            
            return defaults
    
    def _save_preferences(self):
        """Save preferences to file."""
        try:
            with open(self.config_file, 'w') as f:
                json.dump(self.preferences, f, indent=2)
        except Exception as e:
            console.print(f"Failed to save preferences: {e}", style="red")