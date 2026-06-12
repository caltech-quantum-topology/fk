"""

"Symbolic Abelian Group"

This module provides a lightweight symbolic algebra system for manipulating
linear sums of multiple variables. It supports basic arithmetic
operations, substitution, comparison, and algebraic manipulation needed for ILP
constraint solving in the invariant computation pipeline.

Key features:
- Multi-variate polynomial representation using coefficient arrays
- Arithmetic operations: addition, subtraction, scalar multiplication/division
- Comparison operations for constant expressions
- Variable substitution and symbolic manipulation
- String representation with algebraic notation
"""

import copy
import numpy as np
from string import ascii_lowercase, ascii_uppercase

class Symbol:
    """
    Represents a polynomial expression with multiple variables.
    
    A Symbol stores a polynomial as a coefficient array where:
    - var[0] = constant term
    - var[i] = coefficient of variable x_i for i > 0
    
    Example:
        Symbol(2) creates variable x_2 (represented as [0, 0, 1])
        3*Symbol(1) + 5 creates 3*x_1 + 5 (represented as [5, 3])
    """
    
    def __init__(self, index=0):
        """
        Initialize a Symbol.
        
        Args:
            index: Variable index (0 for constant, i>0 for variable x_i)
                  Creates the elementary symbol x_index with coefficient 1
        """
        self.var = np.zeros(index + 1)  # Coefficient array
        self.var[-1] = 1  # Set coefficient of x_index to 1
    def __add__(self, b):
        """
        Add two symbolic expressions or add a scalar to this expression.
        
        Args:
            b: Symbol object or scalar (int/float) to add
            
        Returns:
            New Symbol representing the sum
            
        Implementation handles coefficient array size differences by zero-padding.
        """
        if isinstance(b, int) or isinstance(b, float):
            # Scalar addition: add to constant term
            sym = copy.deepcopy(self)
            sym.var[0] += b
            return sym
        
        # Symbol addition: align coefficient arrays and add element-wise
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            # Extend self.var to match b.var size
            var = np.concatenate((self.var, np.zeros(s2 - s1)))
            new_var = var + b.var
        elif s2 < s1:
            # Extend b.var to match self.var size
            var = np.concatenate((b.var, np.zeros(s1 - s2)))
            new_var = self.var + var
        else:
            # Same size: direct addition
            new_var = self.var + b.var
        
        # Create result symbol
        sym = Symbol()
        sym.var = new_var
        return sym
    def __radd__(self, b):
        """Right addition: scalar + Symbol."""
        sym = copy.deepcopy(self)
        sym.var[0] += b
        return sym
        
    def __sub__(self, b):
        """
        Subtract symbolic expressions or subtract a scalar from this expression.
        
        Args:
            b: Symbol object or scalar (int/float) to subtract
            
        Returns:
            New Symbol representing the difference
        """
        if isinstance(b, int) or isinstance(b, float):
            # Scalar subtraction: subtract from constant term
            sym = copy.deepcopy(self)
            sym.var[0] -= b
            return sym
        
        # Symbol subtraction: align coefficient arrays and subtract element-wise
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            # Extend self.var to match b.var size
            var = np.concatenate((self.var, np.zeros(s2 - s1)))
            new_var = var - b.var
        elif s2 < s1:
            # Extend b.var to match self.var size
            var = np.concatenate((b.var, np.zeros(s1 - s2)))
            new_var = self.var - var
        else:
            # Same size: direct subtraction
            new_var = self.var - b.var
        
        # Create result symbol
        sym = Symbol()
        sym.var = new_var
        return sym
        
    def __rsub__(self, b):
        """Right subtraction: scalar - Symbol = -(Symbol - scalar)."""
        sym = copy.deepcopy(self)
        sym.var *= -1  # Negate all coefficients
        sym.var[0] += b  # Add scalar to constant term
        return sym
    def __mul__(self, b):
        """
        Multiply Symbol by a scalar (polynomial multiplication not supported).
        
        Args:
            b: Scalar (int/float) to multiply by
            
        Returns:
            New Symbol with all coefficients scaled by b
            
        Raises:
            TypeError: If b is not a scalar (polynomial*polynomial not implemented)
        """
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var *= b  # Scale all coefficients
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")
        
    def __rmul__(self, b):
        """Right multiplication: scalar * Symbol."""
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var *= b
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")
        
    def __truediv__(self, b):
        """
        Divide Symbol by a scalar.
        
        Args:
            b: Scalar (int/float) to divide by
            
        Returns:
            New Symbol with all coefficients divided by b
            
        Raises:
            TypeError: If b is not a scalar
        """
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var /= b  # Scale all coefficients
            return sym
        raise TypeError("Division of a Symbol object by anything other than int or float is not supported!")
    def _validate_comparison(self, b):
        """
        Helper method to validate and extract values for comparison operations.
        
        Args:
            b: Value to compare against (int or Symbol)
            
        Returns:
            Tuple of (self_value, b_value) for comparison
            
        Raises:
            Exception: If comparison is not valid (non-constant symbols, wrong types)
        """
        if not self.is_constant():
            raise Exception("Can only compare constant Symbols!")
        
        a = self.constant()
        
        if isinstance(b, int):
            return a, b
        elif isinstance(b, Symbol):
            if b.is_constant():
                return a, b.constant()
            else:
                raise Exception("Can only compare constant Symbols!")
        else:
            raise Exception('Can only compare Symbol objects to other Symbol objects or integers!')

    def __gt__(self, b):
        """Greater than comparison (only for constant symbols)."""
        a, b_val = self._validate_comparison(b)
        return a > b_val
        
    def __lt__(self, b):
        """Less than comparison (only for constant symbols)."""
        a, b_val = self._validate_comparison(b)
        return a < b_val
        
    def __ge__(self, b):
        """Greater than or equal comparison (only for constant symbols)."""
        a, b_val = self._validate_comparison(b)
        return a >= b_val
        
    def __le__(self, b):
        """Less than or equal comparison (only for constant symbols)."""
        a, b_val = self._validate_comparison(b)
        return a <= b_val
    def subs(self, dictionary):
        """
        Substitute values for variables in the symbolic expression.
        
        Args:
            dictionary: Dict mapping variables to their substitution values
                       Keys can be int indices or Symbol objects
                       Values can be int constants or Symbol expressions
        
        Returns:
            New Symbol with substitutions applied
            
        Raises:
            RuntimeError: If attempting to substitute for a constant
            ValueError: If substitution creates circular dependency
            TypeError: If key is not int or Symbol
        """
        new_var = copy.deepcopy(self.var)
        
        for (key, value) in dictionary.items():
            # Convert key to index
            if isinstance(key, Symbol):
                key_index = key.index()
            elif isinstance(key, int):
                key_index = key
            else:
                raise TypeError('index is neither an int or Symbol!')
                
            # Skip if variable doesn't exist in this expression
            if key_index >= len(self.var):
                continue
                
            # Prevent substitution for constant term
            if key_index == 0:  # Note: should be 'key_index == 0' not 'key == one'
                raise RuntimeError("Attempted to substitute some value for a constant!")
                
            # Check for circular substitution
            if isinstance(value, Symbol) and len(value.var) > key_index and value.var[key_index] != 0:
                raise ValueError('Tautological Substitution!')
            
            # Perform substitution
            if isinstance(value, int):
                # Scalar substitution: add coefficient * scalar to constant term
                new_var[0] = new_var[0] + new_var[key_index] * value
            else:
                # Symbol substitution: add coefficient * value_expression
                s1 = new_var.size
                s2 = value.var.size
                
                if s1 < s2:
                    # Extend new_var to accommodate value.var
                    var = np.concatenate((new_var, np.zeros(s2 - s1)))
                    new_var = var + var[key_index] * value.var
                elif s2 < s1:
                    # Extend value.var to match new_var
                    var = np.concatenate((value.var, np.zeros(s1 - s2)))
                    new_var = new_var + new_var[key_index] * var
                else:
                    # Same size: direct operation
                    new_var = new_var + new_var[key_index] * value.var
                    
            # Remove the substituted variable
            new_var[key_index] = 0
            
        # Create result symbol
        sym = Symbol()
        sym.var = new_var
        return sym
    def as_coefficients_dict(self):
        """
        Get coefficient dictionary representation.
        
        Returns:
            Dict mapping Symbol variables to their coefficients (non-zero only)
        """
        return {Symbol(key): value for key, value in enumerate(self.var) if value != 0}
        
    def constant(self):
        """Get the constant term of the expression."""
        return self.var[0]
        
    def is_constant(self):
        """
        Check if this symbol represents a constant (no variables).
        
        Returns:
            True if expression is constant, False if it contains variables
        """
        for index in range(1, len(self.var)):
            if self.var[index] != 0:
                return False
        return True
        
    def free_symbols(self):
        """
        Get list of free variables in this expression.
        
        Returns:
            List of Symbol objects representing variables with non-zero coefficients
        """
        return [Symbol(index) for index in range(1, len(self.var)) if self.var[index] != 0]
    def __getitem__(self, index):
        """Get coefficient at given index (note: should be __getitem__, not __get_item__)."""
        return self.var[index]
        
    def __eq__(self, value) -> bool:
        """
        Check equality with another Symbol.
        
        Args:
            value: Object to compare with
            
        Returns:
            True if both represent the same polynomial expression
        """
        if isinstance(value, Symbol):
            # Compare coefficient arrays, handling different sizes
            s1 = self.var.size
            s2 = value.var.size
            
            if s1 < s2:
                # Extend self.var with zeros for comparison
                var = np.concatenate((self.var, np.zeros(s2 - s1)))
                for index in range(s2):
                    if var[index] != value.var[index]:
                        return False
            elif s2 < s1:
                # Extend value.var with zeros for comparison
                var = np.concatenate((value.var, np.zeros(s1 - s2)))
                for index in range(s1):
                    if self.var[index] != var[index]:
                        return False
            else:
                # Same size: direct comparison
                for index in range(s1):
                    if self.var[index] != value.var[index]:
                        return False
            return True
        else:
            return False
            
    def __hash__(self):
        """
        Hash function for Symbol objects (enables use in sets/dicts).
        
        Returns:
            Hash based on non-zero trailing coefficients
        """
        # Remove trailing zeros to ensure equivalent polynomials hash the same
        to_hash = list(self.var)
        for index in reversed(range(len(to_hash))):
            if to_hash[index] == 0:
                to_hash.pop(index)
            else:
                break
        return hash(tuple(to_hash))
    def _format_polynomial(self):
        """
        Helper method to format polynomial as algebraic string.
        
        Returns:
            String representation like "3a + 2b - c + 5"
        """
        syms = [''] + list(ascii_lowercase) + list(ascii_uppercase)  # Variable names
        string = ''
        
        for index in range(len(self.var)):
            if self.var[index] != 0:
                try:
                    to_print = int(self.var[index])
                except:
                    to_print = self.var[index]  # Keep as float if conversion fails
                    
                sign = 2 * (to_print > 0) - 1  # Extract sign (+1 or -1)
                
                # Add appropriate connector (+/-)
                if string != '':
                    if sign == 1:
                        string += ' + '
                    elif sign == -1:
                        string += ' - '
                    else:
                        raise Exception('Sign should only be 1 or -1!')
                
                if index == 0:
                    # Constant term
                    string += str(abs(to_print) if string != '' else to_print)
                else:
                    # Variable term
                    if sign < 0 and string == '':
                        string = '-'  # Leading negative sign
                    
                    if abs(to_print) != 1:
                        # Non-unit coefficient: show number + variable
                        string += str(abs(to_print)) + syms[index]
                    else:
                        # Unit coefficient: show only variable
                        string += syms[index]
                        
        return string if string != '' else '0'

    def __repr__(self):
        """String representation for debugging/development."""
        return self._format_polynomial()
        
    def __str__(self):
        """String representation for display."""
        return self._format_polynomial()
    def index(self, a=None):
        """
        Extract variable index from elementary symbol.
        
        Args:
            a: Symbol to get index from (if None, use self)
            
        Returns:
            Integer index of the elementary variable
            
        Raises:
            Exception: If symbol is not elementary (has multiple variables or coefficients ≠ 1)
            
        Note: Elementary symbols have exactly one coefficient equal to 1, all others 0.
        """
        if a is None:
            # Get index from self
            found = False
            for index_ in range(len(self.var)):
                val = self.var[index_]
                if val == 1:
                    if not found:
                        found = True
                        ind = index_
                    else:
                        raise Exception('Tried to index non-elementary symbol!')
                elif val != 0:
                    raise Exception('Tried to index non-elementary symbol!')
            if not found:
                raise Exception('No elementary variable found!')
            return ind
        else:
            # Get index from parameter a
            if isinstance(a, int):
                return a
            
            found = False
            for index_ in range(len(a.var)):
                val = a.var[index_]
                if val == 1:
                    if not found:
                        found = True
                        ind = index_
                    else:
                        raise Exception('Tried to index non-elementary symbol!')
                elif val != 0:
                    raise Exception('Tried to index non-elementary symbol!')
            if not found:
                raise Exception('No elementary variable found!')
            return ind

def solve(symbol, index):
    """
    Solve for a variable in a linear equation: symbol = 0.
    
    Rearranges the equation to isolate the specified variable:
    If symbol = c₀ + c₁x₁ + c₂x₂ + ... + cᵢxᵢ + ... = 0
    Then xᵢ = -(c₀ + c₁x₁ + c₂x₂ + ... + cᵢ₋₁xᵢ₋₁ + cᵢ₊₁xᵢ₊₁ + ...) / cᵢ
    
    Args:
        symbol: Symbol equation to solve (assumed = 0)
        index: Variable index to solve for (int or Symbol)
        
    Returns:
        List containing one Symbol representing the solution expression
        
    Raises:
        TypeError: If index is not int or Symbol
        ZeroDivisionError: If coefficient of variable is zero
        
    Example:
        solve(2*x + 3*y - 6, x) returns [3 - 1.5*y] (representing x = 3 - 1.5*y)
    """
    # Convert Symbol index to integer
    if not isinstance(index, int):
        if isinstance(index, Symbol):
            index = index.index()
        else:
            raise TypeError('index is neither an int or Symbol!')
    
    new_symbol = copy.deepcopy(symbol)
    
    # Check if variable has non-zero coefficient
    if new_symbol.var[index] != 0:
        # Solve: divide all coefficients by -coefficient of target variable
        # This isolates the variable: x_i = (-(other terms)) / coefficient
        new_symbol.var /= (-new_symbol.var[index])
    else:
        raise ZeroDivisionError("Cannot solve for variable with zero coefficient!")
    
    # Remove the target variable from the expression
    new_symbol.var[index] = 0
    
    return [new_symbol]

def symbols(n: int):
    """
    Create a list of symbolic variables.
    
    Args:
        n: Number of symbols to create
        
    Returns:
        List of Symbol objects representing variables x_1, x_2, ..., x_n
        
    Example:
        symbols(3) creates [Symbol(1), Symbol(2), Symbol(3)] representing x_1, x_2, x_3
    """
    return [Symbol(j) for j in range(1, n + 1)]

# Global symbolic constants used throughout the system
one = Symbol()        # Represents constant 1: [1, 0, 0, ...]
zero = one - 1        # Represents constant 0: [0, 0, 0, ...]  
nunity = zero - 1     # Represents constant -1: [-1, 0, 0, ...]

# Creates a symabelgroup expression from a coefficient dictionary, recalling that keys are Symbol objects
def expr_from_dict (dict_):
    expression = 0
    for (key, value) in zip(dict_.keys(), dict_.values()):
        expression += value * key
    return expression

# Applies substitutions to all expressions in a dictionary valued in symabelgroup expressions
def subs(expr):
    dictionary_of_symbolic_expressions, evaluation = expr
    out = dict()
    for key, value in dictionary_of_symbolic_expressions.items():
        try:
            out[key] = value.subs(evaluation)
        except:
            out[key] = value
    return out