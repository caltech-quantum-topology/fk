# CLI Package

This package provides the command-line interface for FK computation using Typer.

## Modules

### `app.py`
Typer application definition and main entry point.

**Objects:**
- `app`: Main Typer application
- `template_app`: Subcommand group for templates
- `history_app`: Subcommand group for history

**Functions:**
- `main(argv)`: CLI entry point

### `commands.py`
All CLI command implementations.

**Commands:**
- `interactive_command`: Interactive wizard mode
- `simple_command`: Quick computation
- `config_command`: Config file processing
- `template_create_command`: Create config templates
- `history_*`: History management commands

## Usage

### Interactive Mode

```bash
# Start interactive wizard (default when no args)
fk

# Explicitly call interactive
fk interactive

# Enhanced mode with progress tracking
fk interactive --enhanced

# Quick mode with minimal prompts
fk interactive --quick
```

### Simple Mode

```bash
# Basic computation
fk simple "[1,1,1]" 2

# With symbolic output
fk simple "[1,-2,3]" 2 --symbolic

# Various braid formats work
fk simple "1,1,1" 2
fk simple "1 1 1" 2
```

### Config File Mode

```bash
# Single config file
fk config myconfig.yaml

# Multiple config files
fk config config1.yaml config2.yaml config3.json

# Batch processing (config with "computations" array)
fk config batch.yaml
```

### Template Management

```bash
# Create default template
fk template create

# Create with custom name
fk template create my_config.yaml

# Overwrite existing
fk template create existing.yaml --overwrite
```

### History Management

```bash
# Show recent computations
fk history show
fk history show --limit 20
fk history show --all

# Search history
fk history search "trefoil"

# Clear history
fk history clear
fk history clear --confirm  # Skip prompt

# Export/import
fk history export backup.json
fk history import backup.json
```

## Braid Format

The CLI accepts braids in multiple formats:

| Format | Example | Description |
|--------|---------|-------------|
| JSON | `"[1, -2, 3]"` | Standard JSON array |
| Comma | `"1,-2,3"` | Comma-separated |
| Space | `"1 -2 3"` | Space-separated |

**Important:** Quote the braid string to prevent shell interpretation of negative signs.

```bash
# Correct
fk simple "[1,-2,3]" 2
fk simple "1,-2,3" 2

# May fail (shell interprets -2 as flag)
fk simple 1,-2,3 2
```

## Output

### JSON Output (Default)

```bash
fk simple "[1,1,1]" 2
```

```json
{
  "terms": [...],
  "metadata": {
    "braid": [1, 1, 1],
    "components": 1,
    "inversion": {...},
    ...
  }
}
```

### Symbolic Output

```bash
fk simple "[1,1,1]" 2 --symbolic
```

```
     2
-q⋅x  + q
```

## Configuration File Format

### Single Computation

```yaml
braid: [1, -2, 3]
degree: 2
name: my_knot
preset: parallel
verbose: true
```

### Batch Processing

```yaml
preset: single thread

computations:
  - name: trefoil
    braid: [1, 1, 1]
    degree: 2

  - name: figure_eight
    braid: [1, -2, 1, -2]
    degree: 3
    preset: parallel
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, computation failure, etc.) |

## Environment

The CLI respects these patterns:

```bash
# Pipe output
fk simple "[1,1,1]" 2 > result.json

# Combine with jq
fk simple "[1,1,1]" 2 | jq '.metadata.components'

# Batch processing
for f in configs/*.yaml; do
  fk config "$f"
done
```

## Help

```bash
# General help
fk --help

# Command-specific help
fk simple --help
fk config --help
fk template create --help
fk history --help
```

## Programmatic Access

For Python scripts, prefer the API directly:

```python
# Instead of subprocess call
from fkcompute import fk

result = fk([1, 1, 1], 2)
```

The CLI is best for:
- Interactive exploration
- Shell scripts
- Quick one-off computations
- Batch processing via config files
