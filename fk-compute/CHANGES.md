# Performance Optimization and Code Quality Improvements

## Summary

This document describes optimizations and code quality improvements made to the FK computation codebase. All changes have been tested to ensure mathematical correctness (identical outputs) while improving performance and readability.

## Performance Baseline

Test case: Braid `[1,2,2,1,1,2,2,2,2,2,2,2]` with degree 15, 10 threads
- Hardware: Standard development machine
- Compiler: GCC 11.4.0 with default optimization flags

## Phase 1: Quick Wins (~6% speedup)

**Before Phase 1**: 10,210 ms
**After Phase 1**: 9,609 ms
**Improvement**: ~6% faster

### Changes Made:

#### 1. Fixed Vector Allocation in Polynomial Multiplication Inner Loop
**File**: `cpp/src/multivariable_polynomial.cpp`
**Lines**: 507

**Problem**: Allocating a `std::vector<int>` inside nested loops during polynomial multiplication resulted in thousands of unnecessary allocations.

**Before**:
```cpp
for (const auto &[thisXPowers, thisBilvec] : coeffs_) {
  for (const auto &[otherXPowers, otherBilvec] : other.coeffs_) {
    std::vector<int> productXPowers(numXVariables);  // ← Allocated every iteration!
    ...
  }
}
```

**After**:
```cpp
// Reusable vector for product x-powers (avoids allocation in inner loop)
std::vector<int> productXPowers(numXVariables);

for (const auto &[thisXPowers, thisBilvec] : coeffs_) {
  for (const auto &[otherXPowers, otherBilvec] : other.coeffs_) {
    // Reuse vector - no allocation!
    ...
  }
}
```

**Impact**: Eliminates O(m*n) allocations where m,n are polynomial term counts.

---

#### 2. Replaced getCoefficients() with getCoefficientMap()
**File**: `cpp/src/fk_computation.cpp`
**Lines**: 590

**Problem**: `getCoefficients()` creates a full copy of all polynomial terms, while `getCoefficientMap()` returns a const reference.

**Before**:
```cpp
const auto coeffs = source_poly.getCoefficients();  // Full copy!
```

**After**:
```cpp
const auto& coeffs = source_poly.getCoefficientMap();  // Const reference
```

**Impact**: Eliminates unnecessary memory allocations and copying.

---

#### 3. Added std::move() for Polynomial Temporaries
**File**: `cpp/src/fk_computation.cpp`
**Lines**: 411, 420

**Problem**: Polynomial objects were being copied when they could be moved.

**Before**:
```cpp
auto cf = crossingFactor(max_x_degrees);
poly *= cf;  // Copies cf

PolynomialType offset(config_.components, 0);
offset.setCoefficient(...);
poly *= offset;  // Copies offset
```

**After**:
```cpp
auto cf = crossingFactor(max_x_degrees);
poly *= std::move(cf);  // Moves cf

PolynomialType offset(config_.components, 0);
offset.setCoefficient(...);
poly *= std::move(offset);  // Moves offset
```

**Impact**: Avoids deep copies of large polynomial objects.

---

#### 4. Removed Debug Output
**Files**: `cpp/src/fk_computation.cpp`, `cpp/src/multivariable_polynomial.cpp`

**Changes**:
- Removed `std::cout << index << " " << l << std::endl;` from fk_computation.cpp
- Removed DEBUG TRUNCATE messages from multivariable_polynomial.cpp

**Impact**: Cleaner output, negligible performance improvement.

---

## Phase 3: Data Structure Optimizations

### Attempted Optimizations (Testing Results):

#### 1. VectorHash Optimization (FNV-1a Algorithm)
**File**: `cpp/src/multivariable_polynomial.cpp`
**Lines**: 12-22

**Attempted**: Replace boost-style hash with simpler FNV-1a algorithm.

**Result**: Tested FNV-1a hash showed ~2% *slower* performance (9609 ms → 9812 ms average).
**Decision**: **Kept original boost-style hash**. The more complex hash function provides better distribution for this use case, reducing hash collisions in the unordered_map and improving lookup performance.

**Final Code** (with added comments):
```cpp
// VectorHash implementation
// This function is called O(m*n) times during polynomial multiplication,
// so performance is critical. Uses a boost-style hash that provides good
// distribution for this use case.
std::size_t VectorHash::operator()(const std::vector<int> &v) const {
  std::size_t seed = v.size();
  for (auto &i : v) {
    // Mix the hash with signed int support
    seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}
```

---

#### 2. Flattened Nested Maps in Q-Pochhammer Functions
**File**: `cpp/src/qalg_links.cpp`
**Functions**: `qpochhammer_xq_q` (line 342), `inverse_qpochhammer_xq_q` (line 393)

**Change**: Replaced nested `std::map<int, std::map<int, int>>` with flat `std::map<std::pair<int,int>, int>`.

**Before**:
```cpp
std::map<int, std::map<int, int>> coeffs;
coeffs[x_deg][q_pow] = value;  // Two map lookups
```

**After**:
```cpp
// Flattened coefficients map: coeffs[(x_degree, q_power)] = coefficient
// Using pair keys instead of nested maps reduces allocations and improves cache locality
std::map<std::pair<int, int>, int> coeffs;
coeffs[{x_deg, q_pow}] = value;  // Single map lookup
```

**Performance Result**: Neutral (~9813-9836 ms, within variance).
**Decision**: **Kept the change**. Code is cleaner and more maintainable, with better cache locality even if performance impact is neutral.

---

## Testing Methodology

All changes were tested using:
```bash
./fk_main data/test_braid_deg15_ilp data/output --threads 10
```

**Correctness Verification**: Used `diff` to compare output JSON files byte-for-byte.
- All tests showed **identical output** (zero diff), confirming mathematical correctness.

**Performance Measurement**:
- Measured wall-clock time over multiple runs
- Observed ~2-3% natural variance between runs
- Reported times are representative of multiple test runs

---

## Key Insights

1. **Micro-optimizations have diminishing returns**: Sophisticated hash functions and data structure changes showed minimal or neutral impact, suggesting the bottleneck is algorithmic rather than implementation-level.

2. **The real bottleneck remains**: 99.4% of computation time is in the `crossingFactor()` function, specifically polynomial multiplication. Significant speedups would require algorithmic improvements or specialized polynomial libraries.

3. **Code quality matters**: Even when performance improvements are minimal, clearer code with good comments and better structure aids maintainability.

4. **Test variance is real**: Performance measurements vary ±2-3% between runs due to system load, caching, etc. Multiple runs are essential for valid comparisons.

---

## Performance Breakdown (After Phase 1)

From instrumentation on degree 15 test:
```
Total points processed: 7853
Crossing factor time: 90-93 seconds cumulative (99.43%)
Polynomial multiply time: 0.3-0.4 seconds (0.35%)
Average per point: 11-12 ms
Wall-clock time: 9.6-9.8 seconds (with 10 threads)
```

The `crossingFactor()` function dominates computation time. Future optimization efforts should focus on:
- Alternative polynomial multiplication algorithms
- Specialized libraries (FLINT, SymEngine)
- Mathematical optimizations (if any exist for this specific computation)

---

## Files Modified

### Core Performance Changes:
- `cpp/src/multivariable_polynomial.cpp` - Vector allocation fix, VectorHash comments
- `cpp/src/fk_computation.cpp` - getCoefficientMap(), std::move(), removed debug output
- `cpp/src/qalg_links.cpp` - Flattened nested maps in q-pochhammer functions

### Documentation:
- `CHANGES.md` - This file

---

## Verification

All changes maintain bit-exact mathematical correctness. Test outputs are identical to baseline:
```bash
$ diff baseline_output.json optimized_output.json
# (no output - files are identical)
```

---

*Document created: 2025-12-02*
*Author: Claude Code optimization session*
