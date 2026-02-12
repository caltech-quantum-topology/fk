"""
Symbol class for linear algebra in FK computation.

This module provides a Symbol class that represents linear expressions
with integer coefficients, used for symbolic manipulation of constraints.
"""

from __future__ import annotations

import copy
from string import ascii_lowercase, ascii_uppercase
from typing import Dict, List, Union, Optional

import numpy as np


def symbols(n: int) -> List['Symbol']:
    """Create n symbols indexed from 1 to n."""
    return [Symbol(j) for j in range(1, n + 1)]


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
            sym = copy.deepcopy(self)
            sym.var[0] += b
            return sym
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            var = np.concatenate((self.var, np.zeros(s2 - s1)))
            new_var = var + b.var
        elif s2 < s1:
            var = np.concatenate((b.var, np.zeros(s1 - s2)))
            new_var = self.var + var
        else:
            new_var = self.var + b.var
        sym = Symbol()
        sym.var = new_var
        return sym

    def __radd__(self, b: Union[int, float]) -> 'Symbol':
        sym = copy.deepcopy(self)
        sym.var[0] += b
        return sym

    def __sub__(self, b: Union['Symbol', int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            sym = copy.deepcopy(self)
            sym.var[0] -= b
            return sym
        s1 = self.var.size
        s2 = b.var.size
        if s1 < s2:
            var = np.concatenate((self.var, np.zeros(s2 - s1)))
            new_var = var - b.var
        elif s2 < s1:
            var = np.concatenate((b.var, np.zeros(s1 - s2)))
            new_var = self.var - var
        else:
            new_var = self.var - b.var
        sym = Symbol()
        sym.var = new_var
        return sym

    def __rsub__(self, b: Union[int, float]) -> 'Symbol':
        sym = copy.deepcopy(self)
        sym.var *= -1
        sym.var[0] += b
        return sym

    def __mul__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            sym = copy.deepcopy(self)
            sym.var *= b
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")

    def __rmul__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            sym = copy.deepcopy(self)
            sym.var *= b
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")

    def __truediv__(self, b: Union[int, float]) -> 'Symbol':
        if isinstance(b, (int, float)):
            sym = copy.deepcopy(self)
            sym.var /= b
            return sym
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
        new_var = copy.deepcopy(self.var)
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
        sym = Symbol()
        sym.var = new_var
        return sym

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
        if isinstance(value, Symbol):
            s1 = self.var.size
            s2 = value.var.size
            if s1 < s2:
                var = np.concatenate((self.var, np.zeros(s2 - s1)))
                for index in range(s2):
                    if var[index] != value.var[index]:
                        return False
            elif s2 < s1:
                var = np.concatenate((value.var, np.zeros(s1 - s2)))
                for index in range(s1):
                    if self.var[index] != var[index]:
                        return False
            else:
                for index in range(s1):
                    if self.var[index] != value.var[index]:
                        return False
            return True
        else:
            return False

    def __hash__(self):
        to_hash = list(self.var)
        for index in reversed(range(len(to_hash))):
            if to_hash[index] == 0:
                to_hash.pop(index)
            else:
                break
        return hash(tuple(to_hash))

    def __repr__(self):
        syms = [''] + list(ascii_lowercase.replace("q", "").replace("x", "").replace("y", "").replace("w", "").replace("z", "")) + list(ascii_uppercase)
        string = ''
        for index in range(len(self.var)):
            if self.var[index] != 0:
                try:
                    to_print = int(self.var[index])
                except:
                    pass
                sign = 2 * (to_print > 0) - 1
                if string != '':
                    if sign == 1:
                        string += ' + '
                    elif sign == -1:
                        string += ' - '
                    else:
                        raise Exception('Sign should only be 1 or -1!')
                if index == 0:
                    string += str(to_print)
                else:
                    if sign < 0 and string == '':
                        string = '-'
                    if abs(to_print) != 1:
                        string += str(abs(to_print)) + syms[index]
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
    new_symbol = copy.deepcopy(symbol)
    if new_symbol.var[index] != 0:
        new_symbol.var /= (-new_symbol.var[index])
    else:
        raise ZeroDivisionError()
    new_symbol.var[index] = 0
    return [new_symbol]
