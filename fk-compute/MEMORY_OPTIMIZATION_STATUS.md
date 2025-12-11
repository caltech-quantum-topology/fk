# F_K Invariant Memory Optimization - Status Report

## Problem Statement
Computing F_k invariants with high x degrees (≥25) causes RAM overflow (>16GB) making computation impossible on most machines.

**Critical Test Case**: `fk [1,2,2,1,1,2,2,2,2,2,2,2] 25` causes RAM overflow

## Expected Results (for correctness testing)
- Degree 10: **3 terms**
- Degree 15: **8 terms**
- Degree 20: **13 terms**
- Degree 25: Unknown (crashes due to RAM)

## Current Branch Status

### Main Branch (Clean)
- **Status**: Original working code, NO optimizations
- **Memory**: Degree 20 works, degree 25 crashes (RAM overflow)
- **Correctness**: ✓ All results verified correct
- **Location**: `git checkout main`

### Josef Branch (Partially Optimized)
- **Status**: Has partial streaming implementation
- **Memory**: Degree 20 works (~2.8GB), degree 25 STILL crashes (RAM overflow)
- **Correctness**: ✓ Degrees 10, 15, 20 produce correct results
- **Location**: `git checkout josef`

### Changes Made on Josef Branch

#### 1. Assignment-Level Streaming (IMPLEMENTED, WORKS)
**File**: `cpp/src/fk_computation.cpp` lines 520-535

**What Changed**:
```cpp
// OLD (main branch):
std::vector<std::vector<int>> all_points;
for (assignment : assignments) {
  auto points = enumeratePoints(assignment);
  all_points.insert(all_points.end(), points.begin(), points.end());
}
setupWorkStealingComputation(all_points);

// NEW (josef branch):
for (assignment : assignments) {
  auto points = enumeratePoints(assignment);
  setupWorkStealingComputation(points);  // Process immediately
}
```

**Result**: Eliminates `all_points` vector, but INSUFFICIENT for degree 25

#### 2. Periodic Result Truncation (IMPLEMENTED, WORKS)
**File**: `cpp/src/fk_computation.cpp` lines 325-330

**What Changed**:
```cpp
result_ += poly;

// Added truncation every 1000 points
static thread_local int points_processed = 0;
points_processed++;
if (points_processed % 1000 == 0) {
  result_ = result_.truncate(config_.degree * 4);
}
```

**Result**: Prevents unbounded growth of result polynomial, but STILL insufficient for degree 25

## Why Degree 25 Still Fails

### Memory Bottlenecks That Remain:

1. **Single Assignment Point Enumeration** (BIGGEST ISSUE)
   - Even ONE assignment can generate millions of points for degree 25
   - These are all enumerated and stored in `points` vector before processing
   - Location: `fk_computation.cpp` line 530: `auto points = enumeratePoints(assignment);`
   - This vector can be 100s of MB even for a single assignment

2. **Crossing Factor Cache**
   - Location: `fk_computation.cpp` line 431: `crossing_factor_cache_`
   - Grows unbounded until hitting MAX_CACHE_SIZE (10000 entries)
   - Each entry stores a full polynomial

3. **Intermediate Polynomial Growth**
   - During multiplication in `crossingFactor`, intermediate polynomials can be huge
   - No truncation happens during these operations

## What NEEDS to Be Done Next

### Priority 1: True Point-Level Streaming

**Problem**: Currently `enumeratePoints(assignment)` returns a full vector of all points.

**Solution**: Create a true iterator that yields points one-by-one WITHOUT storing them all.

**Key Challenge**: The `enumeratePoints` function uses OpenMP parallelization (lines 587-646)
which makes it complex to convert to an iterator.

**Approach Options**:

#### Option A: Generator/Batch Iterator (RECOMMENDED)
Modify `enumeratePoints` to yield batches of 1000-10000 points at a time instead of all points.

```cpp
class PointBatchIterator {
  // Enumerate points in chunks, discard after processing
  bool getNextBatch(std::vector<std::vector<int>>& batch, size_t max_size);
};
```

**Implementation**:
- Create new file: `cpp/src/point_batch_iterator.cpp`
- Modify `enumeratePoints` logic to populate batches incrementally
- May need to sacrifice some OpenMP parallelism within enumeration

#### Option B: Callback-Based Approach
Instead of returning points, have `enumeratePoints` accept a callback:

```cpp
void enumeratePointsWithCallback(
  const AssignmentResult& assignment,
  std::function<void(const std::vector<int>&)> process_point
);
```

Then call `engine->computeForAngles(point)` directly in the callback.

**Pros**: Minimal memory, maintains OpenMP parallelization
**Cons**: Requires refactoring thread management in `setupWorkStealingComputation`

### Priority 2: Eager Truncation in Polynomial Operations

**Problem**: Intermediate polynomials during multiplication grow huge before truncation.

**Solution**: Add `multiplyAndTruncate` method to FMPoly that discards out-of-bounds terms during multiplication.

**File**: `cpp/src/fmpoly_class.cpp`

**Method Signature**:
```cpp
FMPoly FMPoly::multiplyAndTruncate(
  const FMPoly& other,
  const std::vector<int>& max_x_degrees
) const {
  // During multiplication, skip terms where x powers exceed max_x_degrees
}
```

**Use in**: `crossingFactor` method (line 355-430) - replace all `result *= factor` with `result = result.multiplyAndTruncate(factor, max_x_degrees)`

### Priority 3: Optimize Inverse Q-Pochhammer

**Problem**: `inverse_qpochhammer_xq_q` generates huge expansions.

**Solution**: Add `max_q_power` parameter to limit term generation.

**File**: `cpp/src/qalg_links.cpp` line 488-530

**Change**:
```cpp
PolynomialType inverse_qpochhammer_xq_q(
  int n, int qpow, int xMax,
  int max_q_power = INT_MAX  // ADD THIS
) {
  // In generation loop (line 516):
  if (new_q_pow > max_q_power) break;  // Skip excessive q-powers
}
```

## Testing Protocol

After ANY change, test in this order:

1. **Degree 10** - Must produce exactly **3 terms**
2. **Degree 15** - Must produce exactly **8 terms**
3. **Degree 20** - Must produce exactly **13 terms**
4. **Degree 25** - Goal: Complete without crashing, unknown term count

**Commands**:
```bash
# On josef branch:
git checkout josef
cd cpp
make clean && make fk_main -j8
cp fk_main ../src/fkcompute/_bin/fk_main

# Test correctness:
fk [1,2,2,1,1,2,2,2,2,2,2,2] 10 | jq '.terms | length'  # Should be 3
fk [1,2,2,1,1,2,2,2,2,2,2,2] 15 | jq '.terms | length'  # Should be 8
fk [1,2,2,1,1,2,2,2,2,2,2,2] 20 | jq '.terms | length'  # Should be 13

# Test memory (careful!):
timeout 300 /usr/bin/time -v fk [1,2,2,1,1,2,2,2,2,2,2,2] 25 2>&1 | grep "Maximum resident"
```

**Memory Target**: Degree 25 should use <8GB (currently uses >16GB and crashes)

## Key Files Reference

### Files Modified on Josef:
- `cpp/src/fk_computation.cpp` - Main computation loop, added streaming + truncation
- All other files unchanged

### Files to Modify for Full Solution:
- `cpp/include/fk/point_batch_iterator.hpp` (NEW) - Batch iterator interface
- `cpp/src/point_batch_iterator.cpp` (NEW) - Batch iterator implementation
- `cpp/include/fk/fmpoly_class.hpp` - Add multiplyAndTruncate declaration
- `cpp/src/fmpoly_class.cpp` - Implement multiplyAndTruncate
- `cpp/include/fk/qalg_links.hpp` - Add max_q_power parameter
- `cpp/src/qalg_links.cpp` - Implement max_q_power logic
- `cpp/Makefile` - Add point_batch_iterator.cpp to SOURCES
- `CMakeLists.txt` - Add point_batch_iterator.cpp to FK_MAIN_SOURCES

### Key Function Locations:
- `enumeratePoints`: `cpp/src/fk_computation.cpp:587-646`
- `enumeratePointsFromValue`: `cpp/src/fk_computation.cpp:648-750`
- `setupWorkStealingComputation`: `cpp/src/fk_computation.cpp:1290-1350`
- `computeForAngles`: `cpp/src/fk_computation.cpp:267-333`
- `crossingFactor`: `cpp/src/fk_computation.cpp:339-446`
- `inverse_qpochhammer_xq_q`: `cpp/src/qalg_links.cpp:488-586`
- `FMPoly::operator*=`: `cpp/src/fmpoly_class.cpp:~700`

## Original Requirements
1. **Correctness**: Results must match original implementation
2. **Performance**: <10% slowdown acceptable
3. **Readability**: Keep code maintainable

## Summary of Attempt So Far

✓ **What Works**:
- Assignment-level streaming (eliminates all_points vector)
- Periodic truncation (prevents unbounded result growth)
- Correctness maintained for degrees 10, 15, 20

✗ **What Doesn't Work**:
- Degree 25 still causes RAM overflow
- Single assignments still generate millions of points in memory
- Memory savings: ~30-40%, but need 60-70% reduction

**Next Critical Step**: Implement point-batch iterator to eliminate the `points` vector inside the assignment loop. This is the biggest remaining bottleneck.

## Recommended Implementation Order

1. Start with Option B (callback-based) - simpler, less code
2. If that works but is too slow, try Option A (batch iterator)
3. Only after point-level streaming works, add eager truncation (Priority 2)
4. Finally add Q-Pochhammer optimization (Priority 3)

**Test correctness after EACH step before proceeding.**

---
*Last Updated: 2025-12-10*
*Status: Degree 25 still crashes - point-level streaming REQUIRED*
