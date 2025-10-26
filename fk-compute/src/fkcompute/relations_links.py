from fkcompute.utils import sort_any
from fkcompute.states import StateLiteral, ZERO_STATE, NUNITY_STATE

from string import ascii_lowercase, ascii_uppercase
import copy
import numpy as np
import itertools

import gurobipy as gp
from gurobipy import GRB

class Leq:
    def __init__(self, first, second):
        self.first = first
        self.second = second

    def variables(self):
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} <= {self.second}'

    def __eq__(self, other):
        if isinstance(other, Leq):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))
    
class Less:
    def __init__(self, first, second):
        self.first = first
        self.second = second

    def variables(self):
        return [self.first, self.second]

    def __repr__(self):
        return f'Inequality {self.first} < {self.second}'

    def __eq__(self, other):
        if isinstance(other, Less):
            return self.first == other.first and self.second == other.second
        return False

    def __hash__(self):
        return hash((self.first, self.second))

class Zero:
    def __init__(self, state):
        self.state = state

    def variables(self):
        return [self.state]

    def __repr__(self):
        return f'Zero {self.state} = [0]'

    def __eq__(self, other):
        if isinstance(other, Zero):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)
    
class Nunity:
    def __init__(self, state):
        self.state = state

    def variables(self):
        return [self.state]

    def __repr__(self):
        return f'Zero {self.state} = [-1]'

    def __eq__(self, other):
        if isinstance(other, Nunity):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)

class Alias:
    def __init__(self, state, alias):
        state, alias = sort_any([state, alias])
        self.state = state
        self.alias = alias

    def variables(self):
        return [self.state, self.alias]

    def __repr__(self):
        return f'Alias {self.alias} := {self.state}'

    def __eq__(self, other):
        if isinstance(other, Alias):
            return self.state == other.state and self.alias == other.alias
        return False

    def __hash__(self):
        return hash((self.state, self.alias))

class Conservation:
    # sum of inputs = sum of outputs
    def __init__(self, inputs, outputs):
        inputs = sort_any(inputs)
        outputs = sort_any(outputs)
        inputs, outputs = sort_any([inputs, outputs])
        self.inputs = inputs
        self.outputs = outputs

    def variables(self):
        return self.inputs + self.outputs

    def try_sum_alias(self):
        
        # case a + b = c + d (no variable is bound)
        if len(self.inputs) != 1 and len(self.outputs) != 1:
            return None
        
        # case a = b + c, return a
        elif len(self.inputs) == 1 and type(self.inputs[0]) != StateLiteral:
            return self.inputs[0], self.outputs

        # case a + b = c, return c
        elif len(self.outputs) == 1 and type(self.outputs[0]) != StateLiteral:
            return self.outputs[0], self.inputs

    def __repr__(self):
        input_sum = " + ".join(str(x) for x in self.inputs)
        output_sum = " + ".join(str(x) for x in self.outputs)
        return f'Conservation {input_sum} = {output_sum}'

    def __eq__(self, other):
        if isinstance(other, Conservation):
            return self.inputs == other.inputs and self.outputs == other.outputs
        return False

    def __hash__(self):
        return hash((tuple(self.inputs), tuple(self.outputs)))

def symbols(n:int):
    return [Symbol(j) for j in range(1, n + 1)]

class Symbol():
    def __init__(self, index=0):
        self.var = np.zeros(index + 1)
        self.var[-1] = 1

    def __add__(self, b):
        if isinstance(b, int) or isinstance(b, float):
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
    def __radd__(self, b):
        sym = copy.deepcopy(self)
        sym.var[0] += b
        return sym
    def __sub__(self, b):
        if isinstance(b, int) or isinstance(b, float):
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
    def __rsub__(self, b):
        sym = copy.deepcopy(self)
        sym.var *= -1
        sym.var[0] += b
        return sym
    def __mul__(self, b):
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var *= b
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")
    def __rmul__(self, b):
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var *= b
            return sym
        raise TypeError("Multiplication of a Symbol object by anything other than int or float is not supported!")
    def __truediv__(self, b):
        if isinstance(b, int) or isinstance(b, float):
            sym = copy.deepcopy(self)
            sym.var /= b
            return sym
        raise TypeError("Division of a Symbol object by anything other than int or float is not supported!")
    def __gt__(self, b):
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
    def __lt__(self, b):
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
    def __ge__(self, b):
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
    def __le__(self, b):
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
    def subs(self, dictionary):
        new_var = copy.deepcopy(self.var)
        for (key, value) in dictionary.items():
            if self.index(key) < len(self.var):
                if not isinstance(key, int):
                    if isinstance(key, Symbol):
                        key = key.index()
                else:
                    raise TypeError('index is neither an int of Symbol!')
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
    def as_coefficients_dict(self):
        return {Symbol(key):value for key, value in enumerate(self.var) if value != 0}
    def constant(self):
        return self.var[0]
    def is_constant(self):
        for index in range(1, len(self.var)):
            if self.var[index] != 0:
                return False
        return True
    def free_symbols(self):
        return [Symbol(index) for index in range(1, len(self.var)) if self.var[index] != 0]
    def __get_item__(self, index):
        return self.var[index]
    def __eq__(self, value) -> bool:
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
        syms = [''] + list(ascii_lowercase.replace("q","").replace("x","").replace("y","").replace("w","").replace("z","")) + list(ascii_uppercase)
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
        syms = [''] + list(ascii_lowercase.replace("q","").replace("x","").replace("y","").replace("w","").replace("z","")) + list(ascii_uppercase)
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
    def index(self, a=None): # converts symbol into an index
        if a == None:
            bool = False
            for index_ in range(len(self.var)):
                val = self.var[index_]
                if not bool and val == 1:
                    bool = True
                    ind = index_
                elif bool and val == 1:
                    raise Exception('Tried to index non-elementary symbol!')
            return ind
        else:
            if isinstance(a, int):
                return a
            bool = False
            for index_ in range(len(a.var)):
                val = a.var[index_]
                if val != 1 and val != 0:
                    raise Exception('Tried to index non-elementary symbol!')
                elif not bool and val == 1:
                    bool = True
                    ind = index_
                elif bool and val == 1:
                    raise Exception('Tried to index non-elementary symbol!')
            return ind

one = Symbol()
zero = one - 1
nunity = zero - 1

def solve(symbol, index):
    if not isinstance(index, int):
        if isinstance(index, Symbol):
            index = index.index()
        else:
            raise TypeError('index is neither an int of Symbol!')
    new_symbol = copy.deepcopy(symbol)
    if new_symbol.var[index] != 0:
        new_symbol.var /= (-new_symbol.var[index])
    else:
        raise ZeroDivisionError()
    new_symbol.var[index] = 0
    return [new_symbol]

def print_conservations(relations, assignment=dict()):
    printed = []
    for r in relations:
        if isinstance(r, Conservation):
            sum_alias = r.try_sum_alias()
            if sum_alias is not None:
                alias, sum_vars = sum_alias
                if all([v in list(assignment.keys()) for v in valid_list]):
                    expression = assignment[alias] - assignment[sum_vars[0]] - assignment[sum_vars[1]]
            else:
                valid_list = r.inputs + r.outputs
                if all([v in list(assignment.keys()) for v in valid_list]):
                    expression = assignment[r.inputs[0]] + assignment[r.inputs[1]] - assignment[r.outputs[0]] - assignment[r.outputs[1]]
            if expression not in printed:
                print("0", "=", expression)
                printed.append(expression)


def print_inequalities(relations, assignment=dict()):
    printed = []
    for r in relations:
        if isinstance(r, Less):
            if r.first == ZERO_STATE:
                first = zero
            elif r.first == NUNITY_STATE:
                first = nunity
            else:
                if r.first in assignment.keys():
                    first = assignment[r.first]
            if r.second == ZERO_STATE:
                second = zero
            elif r.second == NUNITY_STATE:
                second = nunity
            else:
                if r.second in assignment.keys():
                    second = assignment[r.second]
            if second - first not in printed:
                print(0, "<", second - first)
                printed.append(second - first)
        elif isinstance(r, Leq):
            if r.first == ZERO_STATE:
                first = zero
            elif r.first == NUNITY_STATE:
                first = nunity
            else:
                if r.first in assignment.keys():
                    first = assignment[r.first]
            if r.second == ZERO_STATE:
                second = zero
            elif r.second == NUNITY_STATE:
                second = nunity
            else:
                if r.second in assignment.keys():
                    second = assignment[r.second]
            if second - first not in printed:
                print(0, "<=", second - first)
                printed.append(second - first)

def print_angle_inequalities(criteria, relations, mapping, angle_assignment, n_crossings, assignment):

    acc = np.zeros(1 + n_crossings)
    dict_ = criteria[0].as_coefficients_dict()
    for key in dict_:
        acc += angle_assignment[tuple(mapping[key])] * dict_[key]
    # print("\nAngle Criterion:")  # Commented out to avoid unwanted output
    acc[1:] *= -1
    # print(acc)  # Commented out to avoid unwanted output
    # print("\n")  # Commented out to avoid unwanted output
        
    # quit()

    printed = []
    for r in relations:
        if isinstance(r, Less) or isinstance(r, Leq):
            if r.second == NUNITY_STATE:
                second = nunity
            elif r.second == ZERO_STATE:
                second = zero
            else:
                second = assignment[r.second]
            if r.first == NUNITY_STATE:
                first = nunity
            elif r.first == ZERO_STATE:
                first = zero
            else:
                first = assignment[r.first]
            dict_ = (second - first).as_coefficients_dict()
            acc = np.zeros(1 + n_crossings)
            for key in dict_:
                acc += angle_assignment[tuple(mapping[key])] * dict_[key]
            # print(expression_minimum(criteria[0], singlesigns))
            if list(acc) not in printed:
                # print(acc)  # Commented out to avoid unwanted output
                printed.append(list(acc[:]))

def print_inequality_in_angles (inequality, angle_assignment, n_crossings, mapping):
        dict_ = inequality.as_coefficients_dict()
        acc = np.zeros(1 + n_crossings)
        for key in dict_:
            acc += angle_assignment[tuple(mapping[key])] * dict_[key]
        # print(expression_minimum(criteria[0], singlesigns))
        # print(acc)  # Commented out to avoid unwanted output

def get_zeros(relations):
    return [r.state for r in relations if type(r) == Zero]

def get_nunities(relations):
    return [r.state for r in relations if type(r) == Nunity]

def get_aliases(relations):
    return {r.alias:r.state for r in relations if type(r) == Alias}

def get_sum_aliases(relations):
    return [x[0] for x in [r.try_sum_alias() for r in relations if type(r) == Conservation] if x is not None]

def propagate_zero_aliases(relations, verbose=False):
    zeros = get_zeros(relations)
    aliases = get_aliases(relations)

    new_relations = []
    for r in relations:
        if type(r) == Alias:
            if r.state in zeros:
                if verbose:
                    print('reduction:', r,Zero(r.state),'===>',Zero(r.alias))
                if r.alias not in zeros:
                    new_relations.append(Zero(r.alias))
            elif r.alias in zeros:
                if verbose:
                    print('reduction:',r,Zero(r.alias), '===>',Zero(r.state))
                if r.state not in zeros:
                    new_relations.append(Zero(r.state))
            else:
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations

def propagate_nunity_aliases(relations, verbose=False):
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []
    for r in relations:
        if type(r) == Alias:
            if r.state in nunities:
                if r.alias not in nunities:
                    new_relations.append(Nunity(r.alias))
            elif r.alias in nunities:
                if verbose:
                    print('reduction:',r,Nunity(r.alias), '===>',Nunity(r.state))
                if r.state not in nunities:
                    new_relations.append(Nunity(r.state))
            else:
                new_relations.append(r)
        else:
            new_relations.append(r)

    return new_relations

def de_alias_inequalities(relations, verbose=False):
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)

    new_relations = []
    
    for r in relations:
        if type(r) == Leq:
            f = r.first
            s = r.second
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in nunities:
                f = NUNITY_STATE
            elif r.first in aliases:
                f = aliases[r.first]
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in nunities:
                s = NUNITY_STATE
            elif r.second in aliases:
                s = aliases[r.second]
            r = Leq(f,s)
        elif type(r) == Less:
            f = r.first
            s = r.second
            if r.first in zeros:
                f = ZERO_STATE
            elif r.first in aliases:
                f = aliases[r.first]
            elif r.first in nunities:
                f = NUNITY_STATE
            if r.second in zeros:
                s = ZERO_STATE
            elif r.second in aliases:
                s = aliases[r.second]
            elif r.second in nunities:
                s = NUNITY_STATE
            r = Less(f,s)

        new_relations.append(r)
        
    return new_relations

def symmetric_inequality(relations, verbose=False):
    new_relations = [r for r in relations]
    inequalities = [r for r in relations if type(r) == Leq]

    for r1 in inequalities:
        for r2 in inequalities:
            if r1.second == r2.first and r2.second == r1.first:
                if r1.first == ZERO_STATE:
                    if verbose:
                        print(r1,r2,'=====>',Zero(r1.second))
                    new_relations.append(Zero(r1.second))
                elif r1.first == NUNITY_STATE:
                    new_relations.append(Nunity(r1.second))
                elif r1.second == ZERO_STATE:
                    if verbose:
                        print(r1,r2,'=====>',Zero(r1.first))
                    new_relations.append(Zero(r1.first))
                elif r1.second == NUNITY_STATE:
                    new_relations.append(Nunity(r1.first))
                else:
                    if verbose:
                        print(r1,r2,'=====>',Alias(r1.first,r1.second))
                    new_relations.append(Alias(r1.first,r1.second))
    return new_relations

def conservation_zeros(relations, verbose=False):
    zeros = get_zeros(relations)
    new_relations = []
    for r in relations:
        if type(r) == Conservation:
            inputs = [v for v in r.inputs if v != ZERO_STATE and v not in zeros]
            outputs = [v for v in r.outputs if v!= ZERO_STATE and v not in zeros]
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('remove zeros:', r,'=====>',Conservation(inputs,outputs))
            new_relations.append(Conservation(inputs,outputs))
        else:
            new_relations.append(r)
    return new_relations
            
def conservation_alias(relations, verbose=False):
    aliases = get_aliases(relations)
    new_relations = []
    for r in relations:
        if type(r) == Conservation:
            inputs = [aliases.get(v) or v for v in r.inputs]
            outputs = [aliases.get(v) or v for v in r.outputs]
            if verbose and (inputs != r.inputs or outputs != r.outputs):
                print('conservation alias:', r,'=====>',Conservation(inputs,outputs))
            new_relations.append(Conservation(inputs,outputs))
        else:
            new_relations.append(r)
    return new_relations

def unary_conservation_is_alias(relations, verbose=False):
    new_relations = []
    for r in relations:
        if type(r) == Conservation and len(r.inputs) == 1 and len(r.outputs) == 1:
            alias = Alias(r.inputs[0], r.outputs[0])
            new_relations.append(alias)
            if verbose:
                print('unary conservation is alias', r,'=====>',alias)
        else:
            new_relations.append(r)

    return new_relations
                
def is_vacuous(relation):
    if type(relation) == Leq:
        if relation.first == relation.second:
            return True
    if type(relation) == Alias:
        if relation.alias == relation.state:
            return True
    if type(relation) == Zero and relation.state == ZERO_STATE:
        return True
    if type(relation) == Nunity and relation.state == NUNITY_STATE:
        return True

    return False

def remove_vacuous(relations,verbose=False):
    new_relations = []
    for r in relations:
        if not is_vacuous(r):
            new_relations.append(r)
        elif verbose:
            print('removing vacuous:', r)
    return new_relations

def reduce_relations(relations, verbose=False):
    reduction_rules = [
        de_alias_inequalities,
        symmetric_inequality,
        propagate_zero_aliases,
        propagate_nunity_aliases,
        conservation_alias,
        conservation_zeros,
        unary_conservation_is_alias,
        remove_vacuous
    ]
    for rule in reduction_rules:
        relations = rule(relations, verbose=verbose)
        
    relations = list(set(relations))
    return sort_any(relations)

def full_reduce(relations, verbose=False, max_depth=20):
    if max_depth == 0:
        raise Exception('max reduction recursion depth exceeded')
        
    new_relations = reduce_relations(relations,verbose=verbose)
    if set(relations) == set(new_relations):
        return new_relations
    else:
        return full_reduce(new_relations, verbose=verbose, max_depth = max_depth-1)

def all_variables(relations):
    all_variables = []
    for r in relations:
        for v in r.variables():
            if v != ZERO_STATE and v != NUNITY_STATE:
                all_variables.append(v)
    all_variables = list(set(all_variables))
    return all_variables

def free_variables(relations):
    all_vars = all_variables(relations)
    zeros = get_zeros(relations)
    nunities = get_nunities(relations)
    aliases = get_aliases(relations)
    sum_aliases = get_sum_aliases(relations)
    res = []
    for v in all_vars:
        if v not in zeros and v not in nunities and v not in aliases and v not in sum_aliases:
            res.append(v)
    return sort_any(res)

def equivalence_assignment (assignment, braid_states):
    update = {}
    for (key, value) in assignment.items():
        for state in braid_states.state_equivalence_classes[braid_states.get_state(key)]:
            update[state] = value
    assignment.update(update)
    return assignment

def find_expressions(relations, assignment, braid_states):
    expressions = []
    assignment = equivalence_assignment(assignment, braid_states)
    for relation in relations:
        if isinstance(relation, Conservation):
            sum_alias = relation.try_sum_alias()
            expression = None
            if sum_alias is not None:
                alias, sum_vars = sum_alias
                valid_list = [alias] + sum_vars
                if all([v in list(assignment.keys()) for v in valid_list]):
                    if len(sum_vars) > 0: # modifying: bad coding practice
                        expression = assignment[alias] - assignment[sum_vars[0]] - assignment[sum_vars[1]]
            else:
                if len(relation.inputs) > 0:# modifying: bad coding practice
                    valid_list = relation.inputs + relation.outputs
                    if all([v in list(assignment.keys()) for v in valid_list]):
                        expression = assignment[relation.inputs[0]] + assignment[relation.inputs[1]] - assignment[relation.outputs[0]] - assignment[relation.outputs[1]]
            if expression != None and expression != 0:
                expressions.append(expression)
    return expressions
        
def minimal_free(expressions, new={}, verbose=False):
    if not expressions:
        return new
    considering = expressions.pop()
    syms = list(considering.free_symbols()) 
    if syms:
        update = {syms[0]: solve(considering, syms[0])[0]}
        for j in range(len(expressions)):
            expressions[j] = expressions[j].subs(update)
        keys, values = list(new.keys()), list(new.values())
        for (key, value) in zip(keys, values):
            new[key] = value.subs(update)
        new.update(update)

    return minimal_free(expressions, new)

def extend_variable_assignment(reduced_relations, partial_assignment, braid_states, verbose=False):
    assigned = partial_assignment.keys()
    unassigned = [x for x in all_variables(reduced_relations) if x not in assigned]
    if verbose:
        print("EXTENDING VARIABLE ASSIGNMENT")
        print(f"\tpartial assignment: {partial_assignment}")
        print(f"\tunassigned variables: {unassigned}")
    match unassigned:
        case []:
            return partial_assignment
        case _:
            next_assignment = None
            for relation in reduced_relations:
                relation_type = type(relation)
                if relation_type == Zero:
                    state = relation.state
                    if verbose:
                        print (f"\tconsidering zero {state} := 0")
                    if state in assigned:
                        continue
                    else:
                        next_assignment = [state, 0]
                        break
                elif relation_type == Nunity:
                    state = relation.state
                    if verbose:
                        print (f"\tconsidering nunity {state} := -1")
                    if state in assigned:
                        continue
                    else:
                        next_assignment = [state, -1]
                        break
                elif relation_type == Alias:
                    alias = relation.alias
                    state = relation.state
                    if verbose:
                        print(f"\tconsidering alias {alias} := {state}")
                    if state in assigned and alias in unassigned:
                        next_assignment = [alias, partial_assignment[state]]
                        break
                    else:
                        continue      
                elif relation_type == Conservation:
                    sum_alias = relation.try_sum_alias()
                    if sum_alias is not None:
                        alias, sum_vars = sum_alias
                        if verbose:
                            print(f"\tconsidering sum alias {alias} := {sum_vars}")
                        if alias in unassigned and all([x in assigned for x in sum_vars]):
                            if verbose:
                                print(f"restoring sumalias {alias} {sum_vars}")
                            next_assignment = [alias, sum([partial_assignment[x] for x in sum_vars])]
                            break
                        else:
                            continue
                else:
                    continue
            if next_assignment is not None:
                partial_assignment[next_assignment[0]] = next_assignment[1]
                return extend_variable_assignment(reduced_relations, partial_assignment, braid_states, verbose)
            else:
                raise Exception(f"could not find next assignment, even though not all locations have been assigned.\n\tunassigned variables: {unassigned}\n\tassignments: {partial_assignment}")

def violates_relation(assignment, relation, verbose=False):
    relation_type = type(relation)
    if relation_type == Leq:
        first = relation.first
        second = relation.second
        first_type = type(first)
        second_type = type(second)
        if first_type == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if second_type == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
        # if (first_value > second_value):
        #     print("Leq")
        return (first_value > second_value)
    elif relation_type == Less:
        first = relation.first
        second = relation.second
        first_type = type(first)
        second_type = type(second)
        if first_type == StateLiteral:
            first_value = first.state
        else:
            first_value = assignment[first]
        if second_type == StateLiteral:
            second_value = second.state
        else:
            second_value = assignment[second]
        # if (first_value >= second_value):
        #     print('Less')
        return (first_value >= second_value)
    elif relation_type == Alias:
        alias = relation.alias
        state = relation.state
        alias_value = assignment[alias]
        state_value = assignment[state]
        if (alias_value != state_value):
            print('Alias')
        return (alias_value != state_value)
    elif relation_type == Conservation:
        sum_alias = relation.try_sum_alias()
        if sum_alias is not None:
            alias, sum_vars = sum_alias
            alias_value = assignment[alias]
            sum_value = sum([assignment[x] for x in sum_vars])
            if (alias_value != sum_value):
                print("Sum Alias")
            return (alias_value != sum_value)
        else: 
            inputs = relation.inputs
            outputs = relation.outputs
            if verbose:
                pass
            input_sum = sum([assignment[x] for x in inputs])
            output_sum = sum([assignment[x] for x in outputs])
            if verbose:
                pass
            if (input_sum != output_sum):
                print("Conservation")
            return (input_sum != output_sum)
    elif relation_type == Zero:
        state = relation.state
        state_value = assignment[state]
        if (state_value != 0):
            print("Zero")
        return (state_value != 0)
    elif relation_type == Nunity:
        state = relation.state
        state_value = assignment[state]
        if (state_value != -1):
            print("Nunity")
        return (state_value != -1)

def violates_any_relation(assignment, relations, verbose=False):
    return any([violates_relation(assignment, r, verbose) for r in relations])

def symbolic_variable_assignment(relations, braid_states):

    vars = free_variables(relations)
    assignment = dict(zip(vars, symbols(len(vars))))

    assignment = extend_variable_assignment(relations, assignment, braid_states)
    expressions = list(set(find_expressions(relations, assignment, braid_states)))
    minimalizer = minimal_free(expressions, {})
    for key, value in zip(list(assignment.keys()), list(assignment.values())):
        if value != 0 and value != -1:
            assignment[key] = value.subs(minimalizer)

    return assignment

def unreduced_variable_assignment(braid_states):
    return dict(zip(braid_states.state_locations, symbols(len(braid_states.state_locations))))

def subs(expr):
    dictionary_of_symbolic_expressions, evaluation = expr
    out = dict()
    for key, value in dictionary_of_symbolic_expressions.items():
        try:
            out[key] = value.subs(evaluation)
        except:
            out[key] = value
    return out

def R1R2_first (in2, out1):
    return (out1 + in2 + 1) / 4

def R1R2_second (in1, out2, inverted):
    if inverted:
        return (3 * in1 - out2 + 1) / 4
    return (3 * out2 - in1 + 1) / 4

def R3R4_first (in2, out1, inverted):
    if inverted:
        return (in2 - 3 * out1 + 1) / 4
    return (out1 - 3 * in2 + 1) / 4

def R3R4_second(in1, out2):
    return (out2 + in1 + 1) / 4

def generate_criteria(assignment, braid_states, exp_signs):
    criteria = {}
        
    conditions = {val:zero for val in range(braid_states.n_components)}
    for index in range(0, braid_states.n_strands): # starting at 0 seems to be correct for links with >= 2 components
       conditions[braid_states.closed_strand_components[index]] -= 1/2

    for index in range(braid_states.n_crossings):
        crossing_type = braid_states.r_matrices[index]
        in1 = braid_states.top_input_state_locations[index]
        in2 = (in1[0] + 1, in1[1])
        out1 = (in1[0], in1[1] + 1)
        out2 = (out1[0] + 1, out1[1])
        in1 = braid_states.get_state(in1)
        in2 = braid_states.get_state(in2)
        out1 = braid_states.get_state(out1)
        out2 = braid_states.get_state(out2)
        if crossing_type == "R1" or crossing_type == "R2":
            conditions[braid_states.top_crossing_components[index]] += R1R2_first(assignment[in2],
                                                                                  assignment[out1])
            conditions[braid_states.bottom_crossing_components[index]] +=  R1R2_second(assignment[in1],
                                                                                       assignment[out2],
                                                                                       exp_signs[braid_states.bottom_crossing_components[index]])
        elif crossing_type == "R3" or crossing_type == "R4":
            conditions[braid_states.top_crossing_components[index]] -= R3R4_first(assignment[in2],
                                                                                  assignment[out1],
                                                                                  exp_signs[braid_states.top_crossing_components[index]]) 
            conditions[braid_states.bottom_crossing_components[index]] -= R3R4_second(assignment[in1],
                                                                                      assignment[out2])
        else:
            raise Exception("Crossing type is not one of the four acceptable values: 'R1', 'R2', 'R3', or 'R4'.")

    return criteria


def minimum_degree_symbolic(assignment, braid_states, verbose=False):
    conditions = {val:zero for val in range(braid_states.n_components)}
    for index in range(0, braid_states.n_strands): # starting at 0 seems to be correct for links with >= 2 components
        if verbose:
            print(braid_states.closed_strand_components)
        conditions[braid_states.closed_strand_components[index]] -= 1/2
    for index in range(braid_states.n_crossings):
        crossing_type = braid_states.r_matrices[index]
        in1 = braid_states.top_input_state_locations[index]
        in2 = (in1[0] + 1, in1[1])
        out1 = (in1[0], in1[1] + 1)
        out2 = (out1[0] + 1, out1[1])
        in1 = braid_states.get_state(in1)
        in2 = braid_states.get_state(in2)
        out1 = braid_states.get_state(out1)
        out2 = braid_states.get_state(out2)
        if crossing_type == "R1" or crossing_type == "R2":
            conditions[braid_states.top_crossing_components[index]] += (assignment[out1] + assignment[in2] + 1) / 4
            conditions[braid_states.bottom_crossing_components[index]] += (3 * assignment[out2] - assignment[in1] + 1) / 4
        elif crossing_type == "R3" or crossing_type == "R4":
            conditions[braid_states.top_crossing_components[index]] -= (3 * assignment[in2] - assignment[out1] + 1) / 4
            conditions[braid_states.bottom_crossing_components[index]] -= (assignment[in1] + assignment[out2] + 1) / 4
        else:
            raise Exception("Crossing type is not one of the four acceptable values: 'R1', 'R2', 'R3', or 'R4'.")
    if verbose:
        for value in conditions.values():
            print(value.var)
    # quit()
    return conditions

def inequality_manager(relations, assignment, braid_states):

    singles = []
    multiples = []

    for inequality in relations:
        if isinstance(inequality, Leq) or isinstance(inequality, Less):
            if inequality.first == ZERO_STATE:
                a = 0
            elif inequality.first == NUNITY_STATE:
                a = -1
            elif isinstance(inequality.first, tuple):
                a = assignment[braid_states.get_state(inequality.first)]
            if inequality.second == ZERO_STATE:
                b = 0
            elif inequality.second == NUNITY_STATE:
                b = -1
            elif isinstance(inequality.second, tuple):
                b = assignment[braid_states.get_state(inequality.second)]
            if isinstance(inequality, Leq):
                c = "<="
            elif isinstance(inequality, Less):
                c = "<"
            if isinstance(a, Symbol):
                a_dict = a.as_coefficients_dict()
            elif a == 0:
                a_dict = {}
            elif a == -1:
                a_dict = {one: -1}
            else:
                raise Exception(f'Expected variable "a" to be a Symbol, 0, or -1 (the latter two being integers), but "a" was {a}!')
            if isinstance(b, Symbol):
                b_dict = b.as_coefficients_dict()
            elif b == 0:
                b_dict = {}
            elif b == -1:
                b_dict = {one: -1}
            else:
                raise Exception(f'Expected variable "b" to be a Symbol, 0, or -1 (the latter two being integers), but "b" was {b}!')
            c_dict = {}
            for key in b_dict.keys():
                if key in a_dict.keys():
                    c_dict[key] = b_dict[key] - a_dict[key]
                else:
                    c_dict[key] = b_dict[key]
            for key in a_dict.keys():
                if key not in b_dict.keys():
                    c_dict[key] = -a_dict[key]
            bad_keys = []
            for (key, value) in zip(c_dict.keys(), c_dict.values()):
                if value == 0:
                    bad_keys.append(key)
            for key in bad_keys:
                c_dict.pop(key)
            if c_dict:
                if len(set(c_dict.keys()) - set([one])) == 1:
                    expression = expr_from_dict(c_dict)
                    singles.append(expression)
                else:
                    expression = expr_from_dict(c_dict)
                    multiples.append(expression)

    return list(set(singles)), list(set(multiples))

def is_pathological(symbols, signs, expression):
    dict_ = expression.as_coefficients_dict()
    for symbol in symbols:
        if symbol not in dict_.keys():
            return True
        elif signs[symbol] * dict_[symbol] < 0:
            return True
    return False

def expr_from_dict (dict_):
    expression = 0
    for (key, value) in zip(dict_.keys(), dict_.values()):
        expression += value * key
    return expression

def breadth_first_search(criteria, multiples, singlesigns, n_components, return_branching=False):
    n_inequalities = len(multiples)
    queue = [[copy.deepcopy(criteria), {integer:[0 for _ in range(n_inequalities)] for integer in range(n_components)}, 0]]
    depth = 0
    paths = []
    while True:
        if queue:
            criteria = queue[0][0]
            if queue[0][2] > depth:
                paths = []
            depth = queue[0][2]
        else:
            return None
        status, data = find_pathos(criteria, multiples, singlesigns)

        if status == 200:
            # print(criteria[0])  # Commented out to avoid unwanted output
            if return_branching:
                return_list = []
                for index in range(len(queue[0][1][0])):
                    if queue[0][1][0][index] > 0:
                        return_list.append(multiples[index])
                return return_list
            return criteria, data

        branching = data
        path = queue[0][1]
        for key, value in branching.items():
            for index in range(n_inequalities):
                if value[index] == True:
                    new_history = copy.deepcopy(path)
                    new_history[key][index] += 1
                    if new_history not in paths:
                        new_criteria = copy.deepcopy(criteria)
                        new_criteria[key] -= multiples[index] / 2
                        queue.append([copy.deepcopy(new_criteria), copy.deepcopy(new_history), depth + 1])
                        paths.append(new_history)

        queue.pop(0)

def find_pathos(joint, multiples, singlesigns):
    pathology_syms = set()
    bounded_syms = set()
    bounding_expressions = set()
    specific_pathologies_components = dict()
    all_vars = set(singlesigns.keys())
    for (key_, value_) in list(joint.items()):
        pathology_syms_ = set()
        pathology_found = False
        dict_ = value_.as_coefficients_dict()
        for key in list(singlesigns.keys()):
            if key in dict_.keys() and singlesigns[key] * dict_[key] < 0:
                pathology_syms_.add(key)
                pathology_found = True
        specific_pathologies_components[key_] = pathology_syms_
        if not pathology_found:
            bounded_syms.update(set(value_.free_symbols()))
            bounding_expressions.add(key_)
        else:
            pathology_syms.update(pathology_syms_)
    
    vanishing = all_vars - bounded_syms
    # print(joint[0], vanishing)  # Commented out to avoid unwanted output
    bounding_inequalities = []
    index = 0
    multiples_ = copy.deepcopy(multiples)
    n_inequalities = len(multiples)
    while index < n_inequalities:
        reset = False
        for v in vanishing:
            syms = multiples_[index].free_symbols()
            if v in syms and not len(set(syms) - bounded_syms - set([v])) and singlesigns[v] * multiples_[index].as_coefficients_dict()[v] < 0:
                vanishing.remove(v)
                bounded_syms.add(v)
                bounding_inequalities.append((v, multiples_[index]))
                reset = True
                multiples_.pop(index)
                n_inequalities -= 1
                break
        if reset:
            index = 0
        else:
            index += 1
        
    if all_vars - bounded_syms == set():
        return [200, [bounding_expressions, bounding_inequalities]]
    
    specific_pathologies_inequalities = []
    for multiple in multiples:
        specific_pathologies_ = set()
        for sym in set(multiple.free_symbols()).intersection(pathology_syms): 
            if multiple.as_coefficients_dict()[sym] * singlesigns[sym] < 0:
                specific_pathologies_.add(sym)
        specific_pathologies_inequalities.append(specific_pathologies_)
    branching = {key:[len(value.intersection(z)) > 0 for z in specific_pathologies_inequalities] for key, value in specific_pathologies_components.items()}
    return [404, branching]

# use GB for this, too, if you want to go back to using the symbolic relation reduction
def breadth_first_search_reduction(base, reduction, singlesigns, max_depth=5):
    n_inequalities = len(reduction)
    queue = [[copy.deepcopy(base), [0 for _ in range(n_inequalities)], 0]]
    depth = 0
    paths = []
    while True:
        if queue:
            criteria = queue[0][0]
            if queue[0][2] > depth:
                paths = []
            depth = queue[0][2]
        else:
            return False
        status, data = find_pathos_reduction(criteria, reduction, singlesigns)

        if status == 200:
            return True

        if max_depth == -1 or depth < max_depth:
            branching = data
            path = queue[0][1]
            for index in range(n_inequalities):
                if branching[index] == True:
                    new_history = copy.deepcopy(path)
                    new_history[index] += 1
                    if new_history not in paths:
                        new_criteria = copy.deepcopy(criteria)
                        new_criteria -= reduction[index]
                        queue.append([copy.deepcopy(new_criteria), copy.deepcopy(new_history), depth + 1])
                        paths.append(new_history)
        queue.pop(0)

def find_pathos_reduction(joint, multiples, singlesigns):
    pathology_syms = set()
    bounded_syms = set()
    specific_pathologies_components = set()
    
    pathology_syms_ = set()
    pathology_found = False
    shift = 0
    dict_ = joint.as_coefficients_dict()
    for key in list(singlesigns.keys()):
        if key in dict_.keys():
            if singlesigns[key] * dict_[key] < 0:
                pathology_syms_.add(key)
                pathology_found = True
            else:
                if singlesigns[key] < 0:
                    shift -= dict_[key]
    specific_pathologies_components = pathology_syms_
    if not pathology_found:
        bounded_syms.update(set(joint.free_symbols()))
    else:
        pathology_syms.update(pathology_syms_)
        
    if not pathology_found:
        if joint.var[0] + shift >= 0:
            return [200, None]
        else:
            return [404, [False for _ in range(len(specific_pathologies_inequalities))]]
    
    specific_pathologies_inequalities = []
    for multiple in multiples:
        specific_pathologies_ = set()
        for sym in set(multiple.free_symbols()).intersection(pathology_syms): 
            if multiple.as_coefficients_dict()[sym] * singlesigns[sym] < 0:
                specific_pathologies_.add(sym)
        specific_pathologies_inequalities.append(specific_pathologies_)
    branching = [len(set(specific_pathologies_components).intersection(z)) > 0 for z in specific_pathologies_inequalities]
    return [404, branching]

def reduce_inequalities_symbolic(multiples, singles, singlesigns, depth=-1):
    for index in reversed(range(len(multiples))):
        multiple = multiples[index]
        condition = breadth_first_search_reduction(base=multiple, reduction=list(set(multiples) - set([multiple])) + singles, singlesigns=singlesigns, max_depth=depth)
        if condition:
            multiples.pop(index)
    return multiples

def process_assignment(assignment, braid_states, relations, exp_signs):
    criteria = generate_criteria(assignment, braid_states, exp_signs)
    singles, multiples = inequality_manager(relations, assignment, braid_states)
    singlesigns = {}
    for entry in singles:
        dict_ = entry.as_coefficients_dict()
        singlesigns[list(set(dict_.keys()) - set([one]))[0]] = list(dict_.values())[0]
    multiples = list(set(multiples))
    # multiples = reduce_inequalities_symbolic(multiples, singles, singlesigns)
    return criteria, multiples, singlesigns

def condition_satisfier(expression_as_coefficients_dict, bound, singlesigns, current=[]):
    out = []
    key = list(expression_as_coefficients_dict.keys())[0]
    sign = singlesigns[key]
    if sign == 1:
        values = range(int(bound // expression_as_coefficients_dict[key] + 1))
    elif sign == -1:
        values = range(-1, int(bound // expression_as_coefficients_dict[key] - 1), -1)
    else:
        raise ValueError("Symbol signs should only be '1' or '-1'!")
    for x in values:
        temp_dict = copy.deepcopy(expression_as_coefficients_dict)
        coeff = temp_dict.pop(key)
        temp_current = copy.deepcopy(current)
        temp_current.append((key, x))
        if len(temp_dict.keys()) == 0:  
            out.append(dict(temp_current))
        else:
            out.append(condition_satisfier(temp_dict, bound - x * coeff, singlesigns, temp_current))
    return out

def expression_minimum(bounding_expression, singlesigns):
    summand = 0
    for key, value in bounding_expression.as_coefficients_dict().items():
        if key == one:
            summand += value
        else:
            if singlesigns[key] == -1:
                summand -= value
    return summand

def bounding(expr, bounding_expression_indices, bounding_inequalities, degree, singlesigns):
    assignments = [{}]
    for index in bounding_expression_indices:
        new_assignments = []
        for index_ in reversed(range(len(assignments))):
            criterion = expr[index].subs(assignments[index_])
            if expression_minimum(criterion, singlesigns) >= degree:
                assignments.pop(index_)
            else:
                criterion_dict = criterion.as_coefficients_dict()
                degree_ = copy.deepcopy(degree)
                if one in criterion_dict.keys():
                    degree_ -= criterion_dict.pop(one)
                if len(criterion_dict.keys()) == 0:
                    new_assignments.append(copy.deepcopy(assignments[index_]))
                else:
                    new_sub_assignments = condition_satisfier(criterion_dict, degree_, singlesigns, [])
                    for _ in range(len(criterion_dict.keys()) - 1):
                        new_sub_assignments = list(itertools.chain(*new_sub_assignments))
                    mass_ = copy.deepcopy(assignments[index_])
                    for mass in new_sub_assignments:
                        mass_.update(mass)
                        new_assignments.append(copy.deepcopy(mass_))
        assignments = copy.deepcopy(new_assignments)
    pvars = [x[0] for x in bounding_inequalities]
    while pvars != []:
        var = pvars[0]
        temp = []
        for inequality in bounding_inequalities:
            if inequality[0] == var:
                ineq = inequality[1]
        if singlesigns[var] == 1:
            for dict_ in assignments:
                upper_bound = ineq.subs(dict_).subs({var: 0}) / abs(ineq.as_coefficients_dict()[var])
                if upper_bound.is_constant():
                    upper_bound = upper_bound.constant()
                else:
                    raise Exception('Upper bound is not a constant symbol!')
                values = range(int(upper_bound + 1))
                for value in values:
                    dict_[var] = value
                    temp.append(copy.deepcopy(dict_))
        elif singlesigns[var] == -1:
            for dict_ in assignments:
                lower_bound = -1 * (ineq.subs(dict_).subs({var: 0})) / abs(ineq.as_coefficients_dict()[var])
                if lower_bound.is_constant():
                    lower_bound = lower_bound.constant()
                else:
                    raise Exception('Lower bound is not a constant symbol!')
                values = range(int(lower_bound - 1), 0)
                for value in values:
                    dict_[var] = value
                    temp.append(copy.deepcopy(dict_))
        else:
            raise ValueError("Symbol signs should only be '1' or '-1'!")
        assignments = temp
        pvars.pop(0)
        for index in reversed(range(len(bounding_inequalities))):
            if bounding_inequalities[0] == var:
                bounding_inequalities.pop(index)
    return assignments

def check_sign_assignment(degree, relations, braid_states, fk_case = {"case": 1}):
    assignment = symbolic_variable_assignment(relations, braid_states)

    if fk_case["case"] == 1:
        exp_signs = {0: False, 1: False}
    elif fk_case["case"] == 2:
        exp_signs = {0: True, 1: False}
    elif fk_case["case"] == 3:
        exp_signs = {0: True, 1: False}
    elif fk_case["case"] == 4:
        exp_signs = {0: False, 1: True}
    elif fk_case["case"] == 5:
        exp_signs = {0: False, 1: True}

    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations, exp_signs)
    for value in criteria.values():
        multiples.append(degree - value)
    if not integral_bounded(multiples, singlesigns):
        return None
    return {
            "criteria" : criteria,          # inequalities from degree bounding
            "multiples" : multiples,        
            "single_signs" : singlesigns, 
            "assignment" : assignment,       # symbolic assignment
        }

def czech_sign_assignment(degree, relations, braid_states, verbose=False, fk_case = {"case": 1}):
    assignment = symbolic_variable_assignment(relations, braid_states)

    if fk_case["case"] == 1:
        exp_signs = {0: False, 1: False}
    elif fk_case["case"] == 2:
        exp_signs = {0: True, 1: False}
    elif fk_case["case"] == 3:
        exp_signs = {0: True, 1: False}
    elif fk_case["case"] == 4:
        exp_signs = {0: False, 1: True}
    elif fk_case["case"] == 5:
        exp_signs = {0: False, 1: True}
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations, exp_signs)

    if verbose:
        for value in criteria.values():
            print(expression_minimum(value, singlesigns))
    
    for (key, value) in criteria.items():
        invert = exp_signs[key]
        criteria[key] = degree - (1-2*int(invert))*value

    #if not integral_bounded(multiples + list(criteria.values()), singlesigns):
    #    return None
    return {
            "criteria" : criteria,          # inequalities from degree bounding
            "multiples" : multiples,        
            "single_signs" : singlesigns, 
            "assignment" : assignment       # symbolic assignment
        }

def unreduced_information(relations, braid_states):
    assignment = unreduced_variable_assignment(braid_states)
    print("\nBraid:\n")
    print(braid_states.braid)
    print("\nInversion Data:\n")
    print(braid_states.sign_assignment)
    print("\nSymbolic Assignment:\n")
    print(assignment)
    print("\nInequalities:\n")
    print_inequalities(relations, assignment)
    print("\nConservation Relations:\n")
    print_conservations(relations, assignment)
    criteria, multiples, singlesigns = process_assignment(assignment, braid_states, relations)
    print("\nMinimal x-Degree:\n")
    print(criteria[0].as_coefficients_dict())
    print("\nSigns of Each Variable:\n")
    print(singlesigns)

def total_weight(assignment, braid_states):
    acc = 0
    for index in range(braid_states.n_strands):
        acc += assignment[(index, 0)]
    return acc

def condition_assignments(degree, relations, braid_states):
    check = check_sign_assignment(degree, relations, braid_states)
    if check is None:
        return None
    
    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']
    
    solution = breadth_first_search(criteria, multiples, singlesigns, braid_states.n_components, return_branching=False)
    print(solution)

    expr, data = solution
    bounding_expression_indices, bounding_inequalities = data
    minimal_assignments = bounding(expr, bounding_expression_indices, bounding_inequalities, degree, singlesigns)
    minimal_assignments = [(assignment, dict(t)) for t in {tuple(d.items()) for d in minimal_assignments}]
    assignments = [subs(mass) for mass in minimal_assignments]
    filtered = [assignment for assignment in assignments if not violates_any_relation(assignment, relations)]
    return filtered


def print_symbolic_relations(degree, relations, braid_states, write_to=None, verbose=False):
    """Print the reduced relations in human-readable SymPy format"""

    check = czech_sign_assignment(degree, relations, braid_states)
    if check is None:
        print("No valid sign assignment exists for this degree!")
        return None

    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']

    print(f"=== SYMBOLIC VARIABLES RELATIONS AT DEGREE {degree} ===\n")

    # Print braid information
    print(f"Braid: {braid_states.braid}")
    print(f"Number of components: {braid_states.n_components}")
    print(f"Writhe: {braid_states.writhe}")
    if hasattr(braid_states, 'strand_signs') and braid_states.strand_signs:
        print(f"Inversion data: {braid_states.strand_signs}")
    print()

    # Print variable assignments
    print("=== VARIABLE ASSIGNMENTS ===")
    for state, expr in sorted(assignment.items()):
        if isinstance(expr, int):
            print(f"{state} = {expr}")
        else:
            print(f"{state} = {expr}")
    print()

    # Print variable signs
    print("=== VARIABLE SIGNS ===")
    for var, sign in sorted(singlesigns.items(), key=lambda x: x[0].index() if hasattr(x[0], 'index') else 0):
        sign_str = "+" if sign == 1 else "-"
        print(f"{var}: {sign_str}")
    print()

    # Print degree criteria (constraints for each component)
    print("=== DEGREE CONSTRAINTS ===")
    for component, constraint in sorted(criteria.items()):
        if constraint.is_constant():
            print(f"Component {component}: 0 ≤ {constraint.constant()}")
        else:
            print(f"Component {component}: 0 ≤ {constraint}")
    print()

    # Print inequality constraints from relations
    print("=== RELATION INEQUALITIES ===")
    if multiples:
        for i, ineq in enumerate(multiples):
            if ineq.is_constant():
                print(f"Inequality {i+1}: 0 ≤ {ineq.constant()}")
            else:
                print(f"Inequality {i+1}: 0 ≤ {ineq}")
    else:
        print("No inequality constraints from relations")
    print()

    # Print original relations for context
    print("=== ORIGINAL RELATIONS ===")
    for i, relation in enumerate(relations):
        print(f"{i+1}. {relation}")
    print()

    # Print free variables
    free_vars = [var for var, expr in assignment.items()
                 if isinstance(expr, Symbol) and not expr.is_constant()]
    print("=== FREE VARIABLES ===")
    if free_vars:
        print(f"Free variables: {free_vars}")
    else:
        print("All variables are determined by relations")
    print()

    if write_to:
        # Save to file as well
        import sys
        from io import StringIO

        # Capture the output
        old_stdout = sys.stdout
        sys.stdout = captured_output = StringIO()

        # Re-run the printing to capture it
        print_symbolic_relations(degree, relations, braid_states, None, verbose)

        # Restore stdout and save to file
        sys.stdout = old_stdout
        content = captured_output.getvalue()

        with open(write_to, 'w') as f:
            f.write(content)
        print(f"Relations saved to {write_to}")

    return {
        "criteria": criteria,
        "multiples": multiples,
        "single_signs": singlesigns,
        "assignment": assignment
    }

def ilp(degree, relations, braid_states, fk_case = {"case": 1}, write_to=None, verbose=False):

    check = czech_sign_assignment(degree, relations, braid_states, fk_case)
    if check is None:
        return None

    criteria = check['criteria']
    multiples = check['multiples']
    singlesigns = check['single_signs']
    assignment = check['assignment']

    if verbose:
        for value in criteria.values():
            print(value.var)
    # quit()

    criteria = list(criteria.values())

    n_criteria = len(criteria)
    n_multiples = len(multiples)

    criteria_sizes = [criterion.var.size for criterion in criteria]
    inequality_sizes = [multiple.var.size for multiple in multiples]
    segment_sizes = [segment.var.size for segment in assignment.values() if not isinstance(segment, int)]
    for_fill = max(criteria_sizes + inequality_sizes + segment_sizes)
    criteria_tableau = []
    for index in range(n_criteria):
        if criteria[index].is_constant():
            if criteria[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {criteria[index].constant()}")
        else:
            criteria_tableau.append(np.concat((criteria[index].var, np.zeros(for_fill - criteria_sizes[index]))))
    inequality_tableau = []
    for index in range(n_multiples):
        if multiples[index].is_constant():
            if multiples[index].constant() < 0:
                raise Exception(f"Impossible Inequality 0 <= {multiples[index].constant()}")
        else:
            inequality_tableau.append(np.concat((multiples[index].var, np.zeros(for_fill - inequality_sizes[index]))))
    assignment_tableau = []
    keys = sorted(assignment.keys())
    

    for key in keys:
        value = assignment[key]
        if isinstance(value, int):
            assignment_tableau.append(np.concat(([value], np.zeros(for_fill - 1))))
        else:
            assignment_tableau.append(np.concat((value.var, np.zeros(for_fill - len(value.var)))))

    if verbose:
        print(keys)
        print(assignment[keys[-1]])
        print(np.concat((assignment[keys[-1]].var, np.zeros(for_fill - assignment[keys[-1]].var.size))))
        print(assignment_tableau)
        for (key, value) in sorted(assignment.items()):
            print(key, ":", value)
            # if not isinstance(value, int):
            #     print(key, ":", expression_minimum(value, singlesigns)) # modifying 
            # else:
            #     print(key, ":", value) # modifying 
        print(singlesigns)

        print(len(singlesigns.items()))

        for (key, value) in sorted(assignment.items()):
            if not isinstance(value, int):
                print(key)
                print(value)
                print(np.concat((value.var, np.zeros(for_fill - len(value.var)))))
                print("",end="\n\n")
            else:
                print(key)
                print(value)
                print(np.concat(([value], np.zeros(for_fill - 1))))
                print("",end="\n\n")

    criteria_tableau = np.array(criteria_tableau)
    inequality_tableau = np.array(inequality_tableau)
    assignment_tableau = np.array(assignment_tableau)
    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                if criteria_tableau.shape[0] > 0:
                    criteria_tableau[:, index] *= -1
                    criteria_tableau[:, 0] += criteria_tableau[:, index]
                if inequality_tableau.shape[0] > 0:
                    inequality_tableau[:, index] *= -1
                    inequality_tableau[:, 0] += inequality_tableau[:, index]
                if assignment_tableau.shape[0] > 0:
                    assignment_tableau[:, index] *= -1
                    assignment_tableau[:, 0] += assignment_tableau[:, index]

    vars = sort_any([0] + [x.index() for x in list(singlesigns.keys())])
    if verbose:
        print(len(vars))
    # quit()

    if criteria_tableau.shape[0] > 0:
        criteria_tableau = criteria_tableau[:, vars]
    if inequality_tableau.shape[0] > 0:
        inequality_tableau = inequality_tableau[:, vars]
    if assignment_tableau.shape[0] > 0:
        assignment_tableau = assignment_tableau[:, vars]

    lines = []

    # header
    lines.append(f"{degree},")
    lines.append(f"{braid_states.n_components},")
    lines.append(f"{braid_states.writhe},")

    # braid data
    braid_line = []
    for c in range(braid_states.n_crossings):
        braid_line.append(str(abs(braid_states.braid[c])))
        braid_line.append(str(braid_states.r_matrices[c])[1])
    lines.append(",".join(braid_line)+",")

    # closed strands
    closed_line = [str(braid_states.closed_strand_components[c]) 
                   for c in range(1, braid_states.n_strands)]
    lines.append(",".join(closed_line)+",")

    # top & bottom crossing components
    cross_line = []
    for c in range(braid_states.n_crossings):
        cross_line.append(str(braid_states.top_crossing_components[c]))
        cross_line.append(str(braid_states.bottom_crossing_components[c]))
    lines.append(",".join(cross_line)+",")

    # criteria tableau
    for row in criteria_tableau:
        lines.append(",".join(str(e) for e in row)+",")
    lines.append("/")  # separator

    # inequality tableau
    for row in inequality_tableau:
        lines.append(",".join(str(int(e)) for e in row)+",")
    lines.append("/")  # separator

    # assignment tableau
    for row in assignment_tableau:
        lines.append(",".join(str(int(e)) for e in row)+",")

    # Join everything with newlines
    output = "\n".join(lines)

    # Write once
    if write_to:
        with open(write_to, "w") as f:
            f.write(output)
    return output

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()
# can streamline the readout of singlesigns and multiples by doing so directly from their creation.
def integral_bounded(multiples, singlesigns): # note that rational feasibility does not guarantee integer feasibility, though the former is a prerequisite for the latter
    n_multiples = len(multiples)
    sizes = [multiple.var.size for multiple in multiples]
    for_fill = max(sizes)
    tableau = []
    for index in range(n_multiples):
        if multiples[index].is_constant() and multiples[index].constant() == -1:
            return False
        tableau.append(np.concat((multiples[index].var, np.zeros(for_fill - sizes[index]))))
    tableau = np.array(tableau)
    for index in range(1, for_fill):
        if Symbol(index) in singlesigns.keys():
            if singlesigns[Symbol(index)] == -1:
                tableau[:, index] *= -1
                tableau[:, 0] += tableau[:, index]
    model = gp.Model(env=env)
    model.setParam(gp.GRB.Param.PoolSearchMode, 1)
    x = model.addVars(singlesigns.keys(), vtype=GRB.INTEGER) # you can merge below and above, bt-dubs
    for index in range(n_multiples):
        value = tableau[index][0]
        for index_ in range(1, for_fill):
            if tableau[index][index_] != 0:
                value += tableau[index][index_] * x[Symbol(index_)]
        try:
            model.addConstr(0 <= value)
        except:
            pass
    for key in singlesigns.keys():
        model.setObjective(x[key], sense=GRB.MAXIMIZE)
        model.optimize()
        if model.Status != 2:
            return False
    if model.SolCount == 0:
        return False
    return True

def stratified_condition_assignments(weight, relations, braid_states):
    assignment = symbolic_variable_assignment(relations, braid_states)
    symbolic_weight = {0:total_weight(assignment, braid_states)}
    _, multiples, singlesigns = process_assignment(assignment, braid_states, relations)
    multiples.append(weight - symbolic_weight[0])
    if not integral_bounded(multiples, singlesigns):
        raise Exception("It's impossible to use the total weight to bound the braid segment labels!")
    solution = breadth_first_search(symbolic_weight, multiples, singlesigns, braid_states.n_components)
    expr, data = solution
    bounding_expression_indices, bounding_inequalities = data
    minimal_assignments = bounding(expr, bounding_expression_indices, bounding_inequalities, weight, singlesigns)
    minimal_assignments = [(assignment, dict(t)) for t in {tuple(d.items()) for d in minimal_assignments}]
    assignments = [subs(mass) for mass in minimal_assignments]
    filtered = [assignment for assignment in assignments if not violates_any_relation(assignment, relations)]
    return filtered

# there is a redundancy in the above solver: you are checking all criteria when only those with only nonpathological variables can start a complete bounding

