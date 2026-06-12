# CLI

The CLI entrypoint is `fk` (Typer app in `src/fkcompute/cli/app.py`).

Help:

```bash
fk --help
fk simple --help
fk config --help
```

## Commands

### Interactive

```bash
fk
fk interactive
fk interactive --enhanced
fk interactive --quick
```

Behavior:

- Enhanced mode uses `rich` and lives under `src/fkcompute/interactive/`.
- If `rich` is not installed, the CLI falls back to a basic stdin/stdout prompt loop.

### Simple

```bash
fk simple "[1,-2,3]" 2
```

Options:

- `--symbolic`: print a symbolic polynomial (requires SymPy)
- `--format {pretty,inline,latex,mathematica}`: choose symbolic format (auto-enables symbolic if not `pretty`)

### Config

Run one or more config files:

```bash
fk config one.yaml
fk config a.yaml b.yaml c.json
```

If the config file contains a `computations:` list, it is treated as a batch file (see `docs/config.md`).

### Template

Create a starter YAML file:

```bash
fk template create
fk template create my_config.yaml
fk template create my_config.yaml --overwrite
```

Note: the current template is intentionally minimal; `docs/config.md` contains a more complete option reference.

### Print-As

Format a saved JSON result file as a symbolic expression:

```bash
fk print-as data/trefoil.json --format latex
```

This is useful when you ran a compute with `save_data: true` and want to reformat without recomputing.

### History

History is available when the optional interactive dependencies are installed.

```bash
fk history show
fk history search trefoil
fk history export history.json
fk history import history.json
fk history clear
```

History files live under `~/.fkcompute/`.

## Braid Input Formats

The CLI accepts braids in several string formats (parsed by `src/fkcompute/infra/config.py:parse_int_list`):

- JSON-style: `"[1, -2, 3]"`
- Comma-separated: `"1,-2,3"`
- Space-separated: `"1 -2 3"`

Always quote the argument so `-2` is not parsed as a flag:

```bash
fk simple "1 -2 3" 2
```

## Legacy Shortcut

The CLI supports a legacy invocation:

```bash
fk "[1,-2,3]" 2
```

This is rewritten to `fk simple "[1,-2,3]" 2` by `src/fkcompute/cli/app.py`.

## Output

- Default output is JSON (printed to stdout).
- With `--symbolic` / `--format`, the CLI prints a symbolic expression.

See `docs/output.md` for the exact result schema.

## Man Page

Install and open:

```bash
fk-install-man
man fk
```
