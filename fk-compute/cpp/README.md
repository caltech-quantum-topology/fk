# FK Computation Project

A C++ implementation for computing FK with multivariable polynomial support and arithmetic operations.

## Project Structure

```
├── include/fk/           # Header files (.hpp)
│   ├── bilvector.hpp            # Bidirectional vector template class
│   ├── btree.hpp                # Binary tree template class
│   ├── linalg.hpp               # Linear algebra operations
│   ├── multivariable_polynomial.hpp  # Multivariable polynomial class
│   ├── qalg_links.hpp           # Q-algebra and binomial operations
│   ├── solution_pool_1a_double_links.hpp  # ILP solution pooling
│   └── string_to_int.hpp        # String parsing utilities
├── src/                  # Source files (.cpp)
│   ├── fk_segments_links.cpp    # Main FK computation implementation
│   ├── linalg.cpp               # Linear algebra implementations
│   ├── multivariable_polynomial.cpp  # Polynomial operations
│   ├── qalg_links.cpp           # Q-algebra implementations
│   ├── solution_pool_1a_double_links.cpp  # ILP solver
│   └── string_to_int.cpp        # String parsing implementations
├── tests/                # Test files
│   ├── arithmetic_test.cpp      # Polynomial arithmetic tests
│   └── testing.cpp              # Q-algebra function tests
├── examples/             # Example programs
│   └── polynomial_example.cpp   # Polynomial usage demonstration
├── build/                # Generated object files (created automatically)
├── docs/                 # Documentation (for future use)
├── Makefile              # Build configuration
└── README.md             # This file
```

## Building the Project

### Prerequisites
- C++17 compatible compiler (g++ recommended)
- Make utility

### Quick Start
```bash
# Build main executable and example
make

# Build and run polynomial example
make run-example

# Build and run arithmetic tests
make run-arithmetic

# Build main FK computation
make fk_segments_links

# Clean all generated files
make clean
```

### Available Make Targets

**Build Targets:**
- `all` - Build main executable and example (default)
- `fk_segments_links` - Build main FK computation executable
- `polynomial_example` - Build polynomial example
- `arithmetic_test` - Build arithmetic operations test
- `testing` - Build testing executable

**Debug Builds:**
- `debug` - Build main executable with debug flags
- `debug-example` - Build example with debug flags

**Run Targets:**
- `run` - Build and run main executable
- `run-example` - Build and run polynomial example
- `run-arithmetic` - Build and run arithmetic test
- `run-test` - Build and run testing executable

**Utility Targets:**
- `clean` - Remove all generated files
- `clean-obj` - Remove only object files
- `check` - Verify all files compile
- `format` - Format code with clang-format (if available)
- `structure` - Show project directory structure
- `help` - Show all available targets

## Features

### Multivariable Polynomial Support
- Polynomials in variables q, x₁, x₂, ..., xₙ
- Arbitrary positive/negative q powers using bilvector storage
- Configurable maximum degrees for each x variable
- Sparse representation for efficiency

### Arithmetic Operations
- Addition (`+`, `+=`)
- Subtraction (`-`, `-=`)
- Multiplication (`*`, `*=`)
- Compatibility checking for polynomial operations
- Support for negative coefficients

### Core Components
- **bilvector**: Template class for bidirectional indexing
- **btree**: Binary tree for efficient substring search
- **FK Class**: Main FK computation
- **Q-algebra**: Q-binomial and Pochhammer symbol calculations
- **Linear Algebra**: Matrix operations and transformations
- **ILP Solver**: Integer linear programming for solution pooling

### JSON Export
- Export polynomials to JSON format
- Metadata preservation (variable count, max degrees)
- Compatible with existing data processing workflows

## Example Usage

```cpp
#include "fk/multivariable_polynomial.hpp"

// Create polynomial in q, x1, x2 with max degree 3
MultivariablePolynomial poly(2, 3);

// Set coefficients: 5*q^2*x1^1*x2^2
poly.setCoefficient(2, {1, 2}, 5);

// Arithmetic operations
MultivariablePolynomial poly2(2, 3);
poly2.setCoefficient(1, {0, 1}, 3);

auto sum = poly + poly2;  // Addition
auto product = poly * poly2;  // Multiplication

// Export to JSON
poly.exportToJson("output");
```

## Development

### Code Style
- C++17 standard
- Header-only template classes for performance
- Clear separation of interface (.hpp) and implementation (.cpp)
- Comprehensive error handling with exceptions

### Testing
The project includes comprehensive tests:
- Polynomial arithmetic verification
- Q-algebra function validation
- Error handling and edge cases
- Performance benchmarks

### Contributing
1. Follow the existing code style
2. Add tests for new features
3. Update documentation as needed
4. Use `make format` to maintain consistent formatting