# Task: Implement `BMPoly` (Basic Multivariable Polynomial)

## High-level goal

Implement a new C++ class **`BMPoly`** (Basic Multivariable Polynomial) that:

1. Has **the same public API** as `MultivariablePolynomial`, so that it can be swapped in for that class in `FKComputation` with minimal or no changes.
2. Uses an **internal vector-based representation** for coefficients, similar to how multivariable polynomials are handled in `src/fk_segments_links.cpp`, instead of `std::unordered_map` (or similar) as in `MultivariablePolynomial`.
3. Supports **negative powers of the variables** by tracking a **ground degree**, similar to how `FMPoly` behaves.
4. Has its own **header and source file**.
5. Is covered by **unit tests**, and the whole project **still compiles and tests pass**.

You will:
- Inspect the existing polynomial types (`MultivariablePolynomial`, `FMPoly`, and the polynomials in `src/fk_segments_links.cpp`).
- Design and implement `BMPoly`.
- Add tests.
- Integrate it into the build system so everything compiles.
