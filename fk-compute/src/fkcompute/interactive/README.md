# Interactive Package

This package provides a modern, user-friendly interactive interface for FK computation with rich formatting, progress tracking, and session management. It uses the [Rich](https://rich.readthedocs.io/) library for beautiful terminal output.

## Modules

### `wizard.py`
Main wizard classes for guided FK computations.

**Classes:**

| Class | Description |
|-------|-------------|
| `FKWizard` | Full-featured wizard with menus, history, and settings |
| `QuickWizard` | Simplified wizard for rapid computations |

**Functions:**

| Function | Description |
|----------|-------------|
| `run_enhanced_interactive()` | Launch the full wizard interface |
| `run_quick_interactive()` | Launch quick computation mode |

```python
from fkcompute.interactive import run_enhanced_interactive, run_quick_interactive

# Full wizard with menus
run_enhanced_interactive()

# Quick mode with minimal prompts
run_quick_interactive()
```

### `ui.py`
Rich UI components for enhanced interactive experience.

**Classes:**

| Class | Description |
|-------|-------------|
| `BorderedSection` | Context manager for bordered UI sections |
| `StatusMessage` | Static methods for colored status messages |
| `ValidatedInput` | Enhanced input with validation and feedback |
| `ComputationSummary` | Display parameter summary before execution |

```python
from fkcompute.interactive import BorderedSection, ValidatedInput, StatusMessage

# Create a bordered section
with BorderedSection("Input Section", color="cyan"):
    # Content displayed within borders
    pass

# Display status messages
StatusMessage.success("Computation complete")
StatusMessage.error("Invalid input")
StatusMessage.warning("Large braid detected")
StatusMessage.info("Processing...")

# Get validated braid input
braid = ValidatedInput.get_braid()
degree = ValidatedInput.get_degree(braid=braid)
preset = ValidatedInput.get_preset()
```

### `progress.py`
Progress tracking for FK computations with live updates.

**Classes:**

| Class | Description |
|-------|-------------|
| `FKProgressTracker` | Tracks and displays computation phase progress |
| `StatusMessage` | Status message display utilities |
| `ResourceMonitor` | System resource monitoring (optional) |

```python
from fkcompute.interactive import FKProgressTracker

# Use as context manager for automatic start/stop
with FKProgressTracker() as progress:
    # Inversion phase
    progress.start_inversion(braid, degree)
    # ... computation ...
    progress.complete_inversion(components=2)

    # ILP phase
    progress.start_ilp_generation()
    # ... computation ...
    progress.complete_ilp_generation(constraints=100)

    # FK computation phase
    progress.start_fk_computation(threads=4, total_points=1000)
    for i in range(1000):
        progress.update_fk_progress(i, 1000)
    progress.complete_fk_computation(terms_count=50)
```

### `history.py`
Computation history and session management.

**Classes:**

| Class | Description |
|-------|-------------|
| `ComputationHistory` | Manages computation history with search/export |
| `SessionManager` | Manages user preferences and settings |

```python
from fkcompute.interactive import ComputationHistory, SessionManager

# Computation history
history = ComputationHistory()

# Save a computation
history.save_computation(params, result, computation_time=1.5)

# Get recent computations
recent = history.get_recent_computations(limit=10)

# Search history
results = history.search_computations("trefoil")

# Display in table format
history.display_recent(limit=10)

# Export/import
history.export_history("backup.json")
history.import_history("backup.json")

# Session preferences
session = SessionManager()
session.save_preference("default_preset", "parallel")
preset = session.get_preference("default_preset", "single thread")
```

### `prompts.py`
High-level prompts and menu functions.

**Functions:**

| Function | Description |
|----------|-------------|
| `show_main_menu()` | Display main interactive menu |
| `get_computation_parameters(history, session)` | Collect all parameters interactively |
| `show_help_menu()` | Display help and examples |
| `show_settings_menu(session)` | Display and configure settings |
| `show_search_interface(history)` | Interactive history search |

```python
from fkcompute.interactive import show_main_menu, get_computation_parameters

# Show main menu and get choice
choice = show_main_menu()

# Collect parameters through guided prompts
params = get_computation_parameters(history, session)
if params:
    # User confirmed, proceed with computation
    result = fk(**params)
```

## Wizard Features

### Main Menu Options

1. **New Computation** - Start a fresh FK computation with guided input
2. **Recent Computations** - View and reuse previous computations
3. **Search History** - Search through computation history
4. **Settings** - Configure preferences
5. **Help & Examples** - Learn about FK invariants
6. **Exit** - Leave interactive mode

### Guided Input Flow

The wizard guides users through:

1. **Braid Input** - With validation, examples, and braid analysis
2. **Degree Selection** - With suggestions based on braid complexity
3. **Preset Selection** - Choose optimized configuration or custom
4. **Thread Configuration** - Configure parallelism
5. **Optional Settings** - Name, symbolic output, save preferences
6. **Review & Confirm** - Summary before execution

### Progress Tracking Phases

| Phase | Icon | Description |
|-------|------|-------------|
| Inversion | 🔍 | Computing sign assignments |
| ILP | 📋 | Generating ILP formulation |
| FK Computation | ⚡ | Computing FK invariant |
| Symbolic | 🎨 | Generating symbolic output |
| File Operations | 💾 | Saving/loading data |

## Data Storage

### History File
Location: `~/.fkcompute/history.json`

Stores up to 100 recent computations with:
- Timestamp
- Input parameters
- Result summary
- Computation time
- Success/failure status

### Config File
Location: `~/.fkcompute/config.json`

Stores user preferences:
- Default preset
- Default thread count
- Auto-save setting
- Theme preference

## Usage Examples

### CLI Integration

```bash
# Start enhanced interactive mode
fk interactive --enhanced

# Start quick mode
fk interactive --quick

# Default (when no args)
fk
```

### Programmatic Usage

```python
from fkcompute.interactive import FKWizard

# Create and run wizard
wizard = FKWizard()
wizard.run_wizard()

# Access wizard components
wizard.history.display_recent()
wizard.session.save_preference("theme", "dark")
```

### Custom Progress Integration

```python
from fkcompute import fk
from fkcompute.interactive import FKProgressTracker

# Integrate progress tracker with computation
with FKProgressTracker() as progress:
    progress.start_inversion(braid, degree)
    result = fk(
        braid,
        degree,
        _progress_callback=progress  # Hook into computation
    )
    progress.complete_fk_computation(len(result['terms']))
```

## Dependencies

Requires the `rich` library (optional dependency):

```bash
pip install fkcompute[interactive]
# or
pip install rich
```

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+C` | Cancel current operation |
| `Enter` | Confirm/proceed |
| Numbers | Select menu options |

## Error Handling

The wizard gracefully handles:
- Invalid braid input (with retry prompts)
- Computation failures (with error display)
- Keyboard interrupts (with exit confirmation)
- Missing dependencies (with install instructions)

Failed computations are also saved to history for debugging.
