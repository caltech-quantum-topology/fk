"""
Symbol class for linear algebra in FK computation.

This module provides a Symbol class that represents linear expressions
with integer coefficients, used for symbolic manipulation of constraints.
"""

from __future__ import annotations

from string import ascii_lowercase, ascii_uppercase
from typing import Dict, List, Union, Optional

import numpy as np


def symbols(n: int) -> List['Symbol']:
    """Create n symbols indexed from 1 to n."""
    return [Symbol(j) for j in range(1, n + 1)]


def _from_var(var: np.ndarray) -> 'Symbol':
    """Construct a Symbol directly from a coefficient array (no validation)."""
    sym = Symbol.__new__(Symbol)
    sym.var = var
    return sym


class Symbol:
    """
    A symbolic linear expression with integer coefficients.

    Represents expressions like: a + 2*b - 3*c + 5
    where a, b, c are variables indexed by position in the internal array.

    The internal representation is a numpy array where:
    - var[0] is the constant term
    - var[i] for i > 0 is the coefficient of variable i
    """

    def __init__(self, index: int = 0):
        """Create a symbol. index=0 creates a constant 1, index>0 creates variable x_index."""
        self.var = np.zeros(index + 1)
        self.var[-1] = 1

    def __add__(self, b: Union['Symbol', int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            var = self.var.copy()
            var[0] += b
            return _from_var(var)
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            new_var = b.var.copy()
            new_var[:s1] += self.var
        elif s2 < s1:
            new_var = self.var.copy()
            new_var[:s2] += b.var
        else:
            new_var = self.var + b.var
        return _from_var(new_var)

    def __radd__(self, b: Union[int, float]) -> 'Symbol':
        var = self.var.copy()
        var[0] += b
        return _from_var(var)

    def __sub__(self, b: Union['Symbol', int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            var = self.var.copy()
            var[0] -= b
            return _from_var(var)
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            new_var = -b.var
            new_var[:s1] += self.var
        elif s2 < s1:
            new_var = self.var.copy()
            new_var[:s2] -= b.var
        else:
            new_var = self.var - b.var
        return _from_var(new_var)

    def __rsub__(self, b: Union[int, float]) -> 'Symbol':
        var = -self.var
        var[0] += b
        return _from_var(var)

    def __mul__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            return _from_var(self.var * b)
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")

    def __rmul__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            return _from_var(self.var * b)
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")

    def __truediv__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            return _from_var(self.var / b)
        raise TypeError("Division of a Symbol object by anything other than int or float is not supported!")

    def __gt__(self, b: Union['Symbol', int]) -> bool:
        if self.is_constant():
            a = self.constant()
            if isinstance(b, int) or isinstance(b, Symbol):
                if isinstance(b, Symbol):
                    if b.is_constant():
                        b = b.constant()
                    else:
                        raise Exception("Can only compare constant Symbols!")
                return a > b
            else:
                raise Exception('Can only compare Symbol objects to other Symbol objects or integers!')
        else:
            raise Exception("Can only compare constant Symbols!")

    def __lt__(self, b: Union['Symbol', int]) -> bool:
        if self.is_constant():
            a = self.constant()
            if isinstance(b, int) or isinstance(b, Symbol):
                if isinstance(b, Symbol):
                    if b.is_constant():
                        b = b.constant()
                    else:
                        raise Exception("Can only compare constant Symbols!")
                return a < b
            else:
                raise Exception('Can only compare Symbol objects to other Symbol objects or integers!')
        else:
            raise Exception("Can only compare constant Symbols!")

    def __ge__(self, b: Union['Symbol', int]) -> bool:
        if self.is_constant():
            a = self.constant()
            if isinstance(b, int) or isinstance(b, Symbol):
                if isinstance(b, Symbol):
                    if b.is_constant():
                        b = b.constant()
                    else:
                        raise Exception("Can only compare constant Symbols!")
                return a >= b
            else:
                raise Exception('Can only compare Symbol objects to other Symbol objects or integers!')
        else:
            raise Exception("Can only compare constant Symbols!")

    def __le__(self, b: Union['Symbol', int]) -> bool:
        if self.is_constant():
            a = self.constant()
            if isinstance(b, int) or isinstance(b, Symbol):
                if isinstance(b, Symbol):
                    if b.is_constant():
                        b = b.constant()
                    else:
                        raise Exception("Can only compare constant Symbols!")
                return a <= b
            else:
                raise Exception('Can only compare Symbol objects to other Symbol objects or integers!')
        else:
            raise Exception("Can only compare constant Symbols!")

    def subs(self, dictionary: Dict['Symbol', Union['Symbol', int]]) -> 'Symbol':
        """Substitute values for symbols."""
        new_var = self.var.copy()
        for (key, value) in dictionary.items():
            if self.index(key) < len(self.var):
                if not isinstance(key, int):
                    if isinstance(key, Symbol):
                        key = key.index()
                else:
                    raise TypeError('index is neither an int nor Symbol!')
                if key == one:
                    raise RuntimeError("Attempted to substitute some value for a constant!")
                if isinstance(value, Symbol) and value.var[key] != 0:
                    raise ValueError('Tautological Substitution!')
                if isinstance(value, int):
                    new_var[0] = new_var[0] + new_var[key] * value
                else:
                    s1 = new_var.size
                    s2 = value.var.size
                    if s1 < s2:
                        var = np.concatenate((new_var, np.zeros(s2 - s1)))
                        new_var = var + var[key] * value.var
                    elif s2 < s1:
                        var = np.concatenate((value.var, np.zeros(s1 - s2)))
                        new_var = new_var + new_var[key] * var
                    else:
                        new_var = new_var + new_var[key] * value.var
                new_var[key] = 0
        return _from_var(new_var)

    def as_coefficients_dict(self) -> Dict['Symbol', float]:
        """Return dictionary mapping symbols to their coefficients."""
        return {Symbol(key): value for key, value in enumerate(self.var) if value != 0}

    def constant(self) -> float:
        """Return the constant term."""
        return self.var[0]

    def is_constant(self) -> bool:
        """Check if this symbol is a constant (no variables)."""
        for index in range(1, len(self.var)):
            if self.var[index] != 0:
                return False
        return True

    def free_symbols(self) -> List['Symbol']:
        """Return list of symbols with non-zero coefficients."""
        return [Symbol(index) for index in range(1, len(self.var)) if self.var[index] != 0]

    def __getitem__(self, index: int) -> float:
        return self.var[index]

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Symbol):
            return False
        a, b = self.var, value.var
        n = min(a.size, b.size)
        return (
            np.array_equal(a[:n], b[:n])
            and not a[n:].any()
            and not b[n:].any()
        )

    def __hash__(self):
        v = self.var.tolist()
        n = len(v)
        while n and v[n - 1] == 0:
            n -= 1
        return hash(tuple(v[:n]))

    def __repr__(self):
        syms = [''] + list(ascii_lowercase.replace("q", "").replace("x", "").replace("y", "").replace("w", "").replace("z", "")) + list(ascii_uppercase)
        string = ''
        for index in range(len(self.var)):
            coeff = self.var[index]
            if coeff == 0:
                continue
            if coeff == int(coeff):
                coeff = int(coeff)
            if string != '':
                string += ' + ' if coeff > 0 else ' - '
            elif coeff < 0:
                string = '-'
            if index == 0:
                string += str(abs(coeff))
            elif abs(coeff) != 1:
                string += str(abs(coeff)) + syms[index]
            else:
                string += syms[index]
        if string == '':
            string = '0'
        return string

    def __str__(self):
        return self.__repr__()

    def index(self, a: Optional['Symbol'] = None) -> int:
        """
        Get the index of this symbol (if elementary) or of symbol a.

        An elementary symbol has exactly one non-zero coefficient (which is 1).
        """
        if a is None:
            bool_found = False
            for index_ in range(len(self.var)):
                val = self.var[index_]
                if not bool_found and val == 1:
                    bool_found = True
                    ind = index_
                elif bool_found and val == 1:
                    raise Exception('Tried to index non-elementary symbol!')
            return ind
        else:
            if isinstance(a, int):
                return a
            bool_found = False
            for index_ in range(len(a.var)):
                val = a.var[index_]
                if val != 1 and val != 0:
                    raise Exception('Tried to index non-elementary symbol!')
                elif not bool_found and val == 1:
                    bool_found = True
                    ind = index_
                elif bool_found and val == 1:
                    raise Exception('Tried to index non-elementary symbol!')
            return ind


# Special constant symbols
one = Symbol()
zero = one - 1
neg_one = zero - 1


def solve(symbol: Symbol, index: Union[int, Symbol]) -> List[Symbol]:
    """
    Solve for a variable in a symbolic expression.

    Given symbol = 0, solve for the variable at index.
    Returns a list containing the solution.
    """
    if not isinstance(index, int):
        if isinstance(index, Symbol):
            index = index.index()
        else:
            raise TypeError('index is neither an int nor Symbol!')
    new_symbol = _from_var(symbol.var.copy())
    if new_symbol.var[index] != 0:
        new_symbol.var /= (-new_symbol.var[index])
    else:
        raise ZeroDivisionError()
    new_symbol.var[index] = 0
    return [new_symbol]
