
# Background & Problem

* **Current model:** Dense in (x), sparse in (q). Coefficients live in a flat `vector<bilvector<int>>` indexed by a row-major mapping of non-negative (x)-exponents ((a_1,\dots,a_n)) with fixed per-variable caps (0\le a_i\le d_i).
* **Limitations:**

  1. **No negative (x)-exponents** → cannot represent Laurent monomials (x_i^{-k}).
  2. **Memory blow-up** when (n) or the degree caps grow: you allocate (\prod_i (d_i+1)) slots even if most are zero.
  3. **Multiplication work scales with dense grid**: tries all slot pairs, many of which are zero.
  4. **Rigid degree window**: terms exceeding caps are silently dropped; not ideal when the natural support is small but unbounded.

# Goal

* Support **negative (x)-exponents** and avoid dense storage/iteration over empty (x)-slots, while preserving the public API used by callers (constructor shape, `get/set/addToCoefficient`, `getQPolynomial`, `+=/-=/*=`, `print`, `exportToJson`).

# New Data Structure (Solution)

* **Key idea:** Make the representation **sparse in (x)**.
* **Storage:**
  `unordered_map<vector<int>, bilvector<int>, VectorHash> coeffs_`

  * **Key:** full (x)-exponent vector (a=(a_1,\dots,a_n)), allowing **negative integers**.
  * **Value:** the existing `bilvector<int>` holding the sparse (q)-polynomial for that monomial.
* **Hash:** custom `VectorHash` that mixes signed ints robustly; iteration order is unspecified (use sort for deterministic output).

# How This Solves the Problems

* **Negative exponents:** directly supported by allowing negative components in the (x)-key.
* **Space efficiency:** memory ∝ **number of occupied (x)-monomials**, not the dense grid size.
* **Speed:** arithmetic loops only over **present keys**; multiplication visits pairs of occupied monomials and convolves their (q)-parts.
* **No hard degree window:** terms are kept whenever produced; optional pruning can be added by policy (see “Extensions”).

# Algorithms (unchanged logic, different loops)

* **Addition/Subtraction:** for each ((a, Q_a(q))) in RHS, accumulate into `coeffs_[a]` (create-on-write), then prune if the resulting (q)-polynomial is all zeros.
* **Multiplication:** nested loop over occupied keys ((a, Q_a)) and ((b, Q_b)); destination key (c=a+b). In (q), do the usual discrete convolution (i+j) with `bilvector<int>`, skipping zero entries.
* **Zero pruning:** after updates, drop any (a) whose (q)-polynomial became identically zero.

# API Impact (caller-visible)

* **Kept intact:** class name `MultivariablePolynomial`; constructors; `getCoefficient`, `setCoefficient`, `addToCoefficient`, `getQPolynomial`, `+=/-=/*=`, `isZero`, `clear`, `print`, `exportToJson`.
* **Compatibility fields:** `getMaxXDegrees()` and `getBlockSizes()` still exist **as advisory** (no longer used for indexing/caps). If callers relied on them to enforce bounds, that behavior changes (see Risks).
* **Removal:** external dense container access (if anyone reached into the old `vector<bilvector<...>>`)—replace with public API or a new `getCoefficientMap()` accessor (read-only).

# File/Build Changes

* Replace the old header with the new one (class name preserved).
* **Delete** or **stop compiling** the old `.cpp` implementation; the replacement is header-only.

# Serialization / Printing

* **JSON change:** instead of implicit dense ordering, we now emit explicit triplets `{ "x": [a1,...,an], "q": j, "c": c }`. This is robust to sparsity and negative exponents.
* **Pretty print:** order is hash-map dependent; for deterministic ordering, gather keys and sort before printing (optional).

# Complexity (qualitative)

* **Before:**

  * Add/Sub: (O(X\cdot Q)) with (X=\prod_i(d_i+1)).
  * Mul: up to (O(X^2\cdot \bar{Q}^2)).
* **After:**

  * Add/Sub: (O(S\cdot Q)), (S=) number of occupied (x)-monomials.
  * Mul: (O(S_A \cdot S_B \cdot \bar{Q}^2)).
  * Typically **much smaller** when support is sparse.

# Risks / Behavioral Differences

* **No implicit truncation:** products are not dropped for exceeding a degree cap (since none is defined). If you need truncation, add an optional **TermFilter** (see Extensions).
* **Ordering:** iteration is not stable; tests that depend on a printed order should sort.
* **Metadata semantics:** anything that expected `maxXDegrees`/`blockSizes` to be authoritative must be updated (or leave them as advisory only).
* **Performance knobs:** `unordered_map` has overhead; if hot, consider Abseil’s `flat_hash_map` or reserve() heuristics.

# Migration Steps

1. Swap in the new header and remove the `.cpp`.
2. Run tests.

   * If some tests depended on dense print order, adjust them to compare sets of terms or sort output.
   * If logic assumed truncation at `maxXDegrees`, either:

     * add a `TermFilter` (below), or
     * clamp after operations.
3. (Optional) Update any internal utilities that read the old JSON to accept the new triplet format.

# Testing Checklist

* Construct with previous constructors; set/get coefficients at positive and **negative** (x)-exponents.
* `+=, -=` with overlapping and disjoint supports.
* `*=` correctness on mixed positive/negative exponents and negative/positive (q)-powers.
* Zero pruning after cancelling terms.
* JSON round-trip (if you have an importer) or at least JSON schema sanity.
* Large sparse case (e.g., (n=8), few nonzeros) to confirm footprint/time improvements.

# Optional Extensions (small patches if you need them)

* **TermFilter (truncation policy):**
  Plug a predicate `keep(a, q, c)` into `addToCoefficient`/`*=` to drop terms outside desired regions (e.g., total degree, per-var min/max, (|q|) bound).
* **Templated coefficients:**
  `template<typename T>` instead of `int` to support bigints, rationals, or mod (p).
* **Deterministic printing:**
  Provide `printSorted(order)` that sorts by lexicographic (a) then by (q).
* **Performance:**
  Reserve map capacity based on expected support; consider thread-local accumulation during `*=` and merge.

---

**One-liner summary:**
We replace the dense (x)-grid with a sparse hash map keyed by the full (possibly negative) (x)-exponent vector, retaining the existing (q)-sparse `bilvector`. This enables negative exponents, slashes memory and runtime on sparse polynomials, and preserves your public API with only advisory metadata changing in meaning.
