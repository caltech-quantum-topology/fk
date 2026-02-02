# Output Package

This package handles output formatting for FK computation results, primarily converting numerical coefficient data into human-readable symbolic polynomials using SymPy.

## Modules

### `symbolic.py`
SymPy-based symbolic output formatting.

**Functions:**

| Function | Description |
|----------|-------------|
| `format_symbolic_output(fk_result, format_type)` | Format as string |
| `print_symbolic_result(fk_result, format_type)` | Print formatted result |
| `matrix_to_polynomial(fk_result)` | Convert to SymPy expression |
| `list_to_q_polynomial(q_polyL)` | Convert q-terms to polynomial |
| `collect_by_x_powers(polynomial)` | Reorganize by x powers |

**Constants:**
- `SYMPY_AVAILABLE`: Boolean indicating if SymPy is installed

```python
from fkcompute.output.symbolic import (
    format_symbolic_output,
    print_symbolic_result,
    matrix_to_polynomial,
    SYMPY_AVAILABLE,
)

# Check if SymPy is available
if not SYMPY_AVAILABLE:
    print("Install SymPy: pip install sympy")

# Format result as string
fk_result = fk([1, 1, 1], 2)
symbolic_str = format_symbolic_output(fk_result, "pretty")
print(symbolic_str)

# Different format types
latex_str = format_symbolic_output(fk_result, "latex")
plain_str = format_symbolic_output(fk_result, "str")

# Print directly
print_symbolic_result(fk_result, format_type="pretty")
```

## Format Types

| Type | Description | Use Case |
|------|-------------|----------|
| `"pretty"` | Unicode pretty-print | Terminal display |
| `"latex"` | LaTeX markup | Documentation, papers |
| `"str"` | Plain string | Programmatic use |

**Examples:**

```python
# Pretty format (terminal)
    2
-q⋅x  + q

# LaTeX format
- q x^{2} + q

# String format
-q*x**2 + q
```

## Variable Naming

Variables are chosen based on the number of link components:

| Components | Variables |
|------------|-----------|
| 1 | `x` |
| 2 | `x, y` |
| 3+ | `a, b, c, ...` (skipping `q`) |

The quantum parameter `q` is always reserved.

## Usage Examples

### Basic Symbolic Output

```python
from fkcompute import fk
from fkcompute.output.symbolic import format_symbolic_output

result = fk([1, 1, 1], 2, symbolic=True)

# Access via metadata
print(result["metadata"]["symbolic"])

# Or format manually
symbolic = format_symbolic_output(result, "pretty")
print(symbolic)
```

### Direct Polynomial Access

```python
from fkcompute.output.symbolic import matrix_to_polynomial

result = fk([1, 1, 1], 2)
poly = matrix_to_polynomial(result)

# poly is a SymPy expression - can manipulate it
from sympy import symbols, expand, simplify

q = symbols('q')
expanded = expand(poly)
evaluated = poly.subs(q, 2)  # Evaluate at q=2
```

### Custom Formatting

```python
from fkcompute.output.symbolic import (
    matrix_to_polynomial,
    collect_by_x_powers,
)
import sympy as sp

result = fk([1, 1], 2)  # Hopf link (2 components)
poly = matrix_to_polynomial(result)

# Collect terms by x powers
collected = collect_by_x_powers(poly)

# Custom LaTeX with options
latex = sp.latex(collected, order='grlex')

# Factor if possible
factored = sp.factor(poly)
```

### For Documentation

```python
from fkcompute import fk
from fkcompute.output.symbolic import format_symbolic_output

# Generate LaTeX for paper
result = fk([1, 1, 1], 2)
latex = format_symbolic_output(result, "latex")

# Insert into document
doc_text = f"""
The FK polynomial for the trefoil is:
$$
{latex}
$$
"""
```

## FK Result Structure

The symbolic formatter expects this structure:

```python
{
    "terms": [
        {
            "x": [power1, power2, ...],  # Powers of x variables
            "q_terms": [
                {"q": power, "c": coefficient},
                ...
            ]
        },
        ...
    ],
    "metadata": {
        "num_x_variables": int,
        "max_x_degrees": [int, ...],
        ...
    }
}
```

## Error Handling

```python
from fkcompute.output.symbolic import format_symbolic_output

# Missing SymPy
if not SYMPY_AVAILABLE:
    # Returns error message string instead of raising
    result = format_symbolic_output(fk_result, "pretty")
    # "Error formatting symbolic output: SymPy is required..."

# Large coefficients (>10000 digits)
# Raises ValueError with descriptive message
```

## Performance Notes

- Symbolic conversion is relatively fast for typical results
- Very large polynomials (many terms) may take longer
- Consider caching symbolic output if reusing
- The `collect_by_x_powers` function improves readability but adds processing time
