# fkcompute

Compute the **FK invariant** for braids using inversion, ILP reduction, and a compiled helper binary.

This package bundles both Python logic and C++ executables to make the computation portable and easy to use from the command line or within Python code.

---

## Features

- Parse braids and compute their FK invariant.
- Uses inversion assignments and ILP reduction internally.
- Calls an optimized C++ helper binary (`fk_segments_links`) compiled at install time.
- Provides both:
  - **CLI** command `fk`
  - **Python API** via `fkcompute.fk`
- From CLI run `fk --help` for further instructions
- From python shell run `fkcompute.fk` for further instructions

---

## Installation

You’ll need a working C++ toolchain (GCC/Clang on Linux/macOS, MSVC on Windows).  
Then install directly:

```bash
pip install .
