# Performance Notes

This project has two main sources of parallelism:

1. Python multiprocessing during inversion search (Phase 1)
2. OpenMP-style threading inside the C++ backend (Phase 3)

The best settings depend on your braid type and whether you can skip phases using precomputed data.

## Key Knobs

Inversion search (Python):

- `max_workers`: number of worker processes
- `chunk_size`: index range per worker task
- `include_flip`: try mirrored braid variants (increases search)
- `max_shifts`: limit cyclic shifts (reduces search)
- `partial_signs`: constrain some signs to shrink the space

C++ compute:

- `threads`: passed to `fk_main --threads`

I/O (avoid recompute):

- `save_data: true` to persist inversion + ILP + final JSON
- `inversion` / `inversion_file` to skip Phase 1
- `ilp_file` to skip Phase 1 and Phase 2

## Presets

If you are using config files, presets are the easiest way to start:

- `single thread`: low overhead, good for small problems
- `parallel`: uses `cpu_count-1` for both workers and threads

Presets live in `src/fkcompute/api/presets.py`.

## Practical Guidance

- If inversion search dominates:
  - Provide `inversion` when you can (or reuse `*_inversion.json` from a prior run).
  - Try `include_flip: false` unless you know you need flips.
  - Limit `max_shifts` for very long braids.
  - Use `partial_signs` to fix known segments.

- If C++ compute dominates:
  - Increase `threads`.
  - Consider skipping Python phases by reusing an `ilp_file`.

- If you see CPU oversubscription:
  - Reduce either `max_workers` or `threads`.
  - Be aware that BLAS implementations (OpenBLAS) may also use threads; setting `OPENBLAS_NUM_THREADS=1` can help when the C++ backend is already parallel.
