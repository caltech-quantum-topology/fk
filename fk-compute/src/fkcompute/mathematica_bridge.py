"""Mathematica/Wolfram Language bridge.

This module implements a tiny JSON-in/JSON-out interface intended to be
called from Mathematica via RunProcess.

Example (shell):
  python -m fkcompute.mathematica_bridge <<<'{"mode":"version"}'
"""

from __future__ import annotations

import json
import sys
import traceback
from dataclasses import asdict, dataclass
from typing import Any, Dict, Optional


@dataclass
class _Error:
    type: str
    message: str
    traceback: str


def _json_default(obj: Any) -> Any:
    # Best-effort conversion for occasional non-JSON scalars.
    try:
        if hasattr(obj, "item"):
            return obj.item()
    except Exception:
        pass
    return str(obj)


def _read_request() -> Dict[str, Any]:
    raw = sys.stdin.read()
    if not raw.strip():
        return {"mode": "version"}
    return json.loads(raw)


def _write_response(payload: Dict[str, Any]) -> None:
    sys.stdout.write(json.dumps(payload, default=_json_default))
    sys.stdout.flush()


def _apply_preset(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    preset = kwargs.pop("preset", None)
    if not preset:
        return kwargs

    from fkcompute.api.presets import PRESETS

    if preset not in PRESETS:
        raise ValueError(f"Unknown preset '{preset}'. Available: {list(PRESETS.keys())}")

    merged = PRESETS[preset].copy()
    merged.update(kwargs)
    return merged


def _handle_compute(req: Dict[str, Any]) -> Any:
    from fkcompute.api.compute import fk
    from fkcompute.output.symbolic import format_symbolic_output, SYMPY_AVAILABLE

    braid = req.get("braid")
    degree = req.get("degree")
    if braid is None or degree is None:
        raise ValueError("'braid' (list[int]) and 'degree' (int) are required")

    kwargs = req.get("kwargs") or {}
    if not isinstance(kwargs, dict):
        raise TypeError("'kwargs' must be an object/dict")

    kwargs = _apply_preset(dict(kwargs))

    # Extract bridge-specific options before calling fk
    format_type = kwargs.pop("format_type", "mathematica")
    full = kwargs.pop("full", False)

    result = fk(braid, int(degree), **kwargs)

    # By default return only the Mathematica expression string.
    # Pass "full": true in the request to get the complete result dict instead.
    if not full and SYMPY_AVAILABLE:
        return format_symbolic_output(result, format_type=format_type)

    # Full mode: attach the expression and return everything
    if SYMPY_AVAILABLE:
        result["expression"] = format_symbolic_output(result, format_type=format_type)

    return result


def _handle_config(req: Dict[str, Any]) -> Any:
    from fkcompute.api.batch import fk_from_config

    config_path = req.get("config")
    if not config_path or not isinstance(config_path, str):
        raise ValueError("'config' (string path) is required")
    return fk_from_config(config_path)


def _handle_version() -> Dict[str, Any]:
    import platform

    try:
        from fkcompute import __version__ as fk_version
    except Exception:
        fk_version = None

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "fkcompute": fk_version,
    }


def main() -> None:
    try:
        req = _read_request()
        mode = req.get("mode", "compute")

        if mode == "compute":
            result = _handle_compute(req)
        elif mode == "config":
            result = _handle_config(req)
        elif mode == "version":
            result = _handle_version()
        else:
            raise ValueError("Unknown mode. Expected one of: compute, config, version")

        _write_response({"ok": True, "result": result})

    except Exception as e:
        err = _Error(
            type=type(e).__name__,
            message=str(e),
            traceback=traceback.format_exc(),
        )
        _write_response({"ok": False, "error": asdict(err)})
        # Still exit non-zero so callers can detect bridge failures.
        raise


if __name__ == "__main__":
    main()
