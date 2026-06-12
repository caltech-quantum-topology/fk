# Mathematica / Wolfram Language (`FkCompute` Paclet)

This repo includes a Wolfram Language wrapper (Paclet) under `mathematica/FkCompute/`.

It calls into Python via a small JSON bridge implemented in `src/fkcompute/mathematica_bridge.py`.

## Requirements

- Mathematica 12+
- Python 3.9+ with `fkcompute` installed
- SymPy if you want symbolic expressions returned (recommended)

## Install The Paclet

From Mathematica:

```wl
PacletInstall[Directory["/path/to/fk-compute/mathematica/FkCompute"]]
Needs["FkCompute`"]
```

## Basic Usage

Compute from a braid + degree:

```wl
res = FkCompute[{1, 1, 1}, 2];
```

Compute from a YAML/JSON config file (same format as the `fk config` CLI):

```wl
res = FkCompute["config.yaml"];
```

## Options

Options are defined in `mathematica/FkCompute/Kernel/FkCompute.wl`.

Common ones:

- "PythonExecutable" -> Automatic | "/path/to/python"
- "Symbolic" -> True|False
- "Threads" -> Automatic|int
- "Verbose" -> True|False
- "Preset" -> "single thread" | "parallel"
- "MaxWorkers" / "ChunkSize" / "IncludeFlip" / "MaxShifts"
- "SaveData" / "SaveDir" / "Name"

Example:

```wl
res = FkCompute[{1, -2, 1, -2}, 3,
  "Symbolic" -> True,
  "Threads" -> 8,
  "Preset" -> "parallel",
  "Verbose" -> True
];
```

## What Gets Returned?

Compute calls (`FkCompute[braid, degree, ...]`):

- If SymPy is available in the Python environment, the bridge returns a string formatted in Mathematica syntax, and the Paclet converts it to a Wolfram Language expression.
- If SymPy is not available, the bridge returns the full JSON result dict, and the Paclet returns it as an Association.

Config calls (`FkCompute["config.yaml"]`) always return the full result dict (Association).

Internally, the bridge supports returning the full dict for compute calls via a `full` flag in the JSON request; the Paclet currently exposes the default behavior.

## Debugging

- Confirm Python executable:

```wl
FkComputeVersion[]
```

You can also set a default Python executable (see `mathematica/FkCompute/Kernel/FkCompute.wl`):

```wl
$FkComputePythonExecutable = "/usr/bin/python3";
```

- If you see Python bridge failures, ensure the selected Python has `fkcompute` installed and can run:

```bash
python -m fkcompute.mathematica_bridge <<<'{"mode":"version"}'
```

- If symbolic formatting fails, install SymPy in the same environment:

```bash
pip install "fkcompute[symbolic]"
```
