# _attic

Files moved here during the 2026-06 refactor instead of being deleted,
because this directory is not under version control.

- `fk_old/` — pre-package legacy scripts (superseded by `src/fkcompute/`)
- `docs.old/` — old markdown docs (superseded by `docs/` Sphinx sources)
- `src-backup-pre-refactor/` — full snapshot of `src/` taken immediately
  before the refactor, for diffing/rollback
- `test_find_sign_assignment_full.py` — pre-refactor copy of the one test
  file that was modified
- `history.py` — interactive session-history module that was never wired
  into the CLI
- `example_symbolic.py` — unused example module

If everything works for a while, this whole directory can be deleted.

## 2026-06-12: C++ backend refactor

- `cpp-backup-pre-refactor/` — full snapshot of `cpp/` sources taken before
  the C++ refactor. Removed from the live tree afterwards:
  - `src/fk_segments_links.cpp` + `include/fk/fk.hpp` — legacy pre-FKComputation main
  - `src/inequality_solver.cpp`, `src/solution_pool_1a_double_links.cpp`,
    `include/fk/btree.hpp` — superseded by logic inside `fk_computation.cpp`
  - `src/{multivariable_polynomial,bmpoly,hmpoly,zmpoly}.cpp` + headers —
    alternative polynomial backends that no longer compiled against the
    fk_main pipeline (missing exportToJson/fmpz methods); FMPoly is the
    only working backend
  - bit-rotted examples/tests that referenced removed APIs
