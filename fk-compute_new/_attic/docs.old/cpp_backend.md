# C++ Backend (`fk_main`)

The heavy FK polynomial computation is implemented in C++ and exposed via the `fk_main` executable.

Python calls this binary in Phase 3 after it has written an input file (custom CSV) in Phase 2.

## Where It Lives

- Main executable source: `cpp/main.cpp`
- Core implementation: `cpp/src/fk_computation.cpp` and friends
- Headers: `cpp/include/fk/*`

## How Python Finds The Binary

Resolver: `src/fkcompute/infra/binary.py:binary_path`.

Resolution order:

1. Use `fk_main` from `PATH` (if present)
2. Otherwise use the packaged binary at `fkcompute/_bin/fk_main`

This is useful if you are developing the C++ backend and want to override the packaged binary.

## Invocation Contract

The C++ program expects base paths without extensions:

```bash
fk_main <input_base> <output_base> [--threads N] [--verbose]
```

It reads:

- `<input_base>.csv`

and writes:

- `<output_base>.json`

Python wrapper:

- `src/fkcompute/infra/binary.py:run_fk_binary(ilp_file_base, output_base, threads=...)`

## Building

### Build Via CMake (recommended for packaging)

From the repo root:

```bash
cmake -B build -S .
cmake --build build
```

The root `CMakeLists.txt` copies the resulting binary to `src/fkcompute/_bin/fk_main`.

### Build Via Makefile (useful for C++ iteration)

```bash
make -C cpp main
```

## OpenMP / Threads

The backend is designed to run with multiple compute threads.

- The Python API parameter `threads` is passed through to `fk_main --threads`.
- OpenMP is detected at build time (see `CMakeLists.txt` and `cpp/Makefile`).

If OpenMP is not available, the code will still compile, but parallel speedups may be limited.

## Manual Smoke Test

If you have an ILP file (for example one generated with `save_data: true`):

```bash
./src/fkcompute/_bin/fk_main data/trefoil_d2_ilp data/trefoil_d2 --threads 4 --verbose
```

This should write `data/trefoil_d2.json`.
