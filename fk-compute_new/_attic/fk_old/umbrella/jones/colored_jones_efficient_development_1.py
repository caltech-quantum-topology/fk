import numpy as np
from sympy import *
import itertools
from bisect import bisect_left
from operator import add, mul
import symengine

def countInversions(arr, n):
    v = []
    for i in range(n):
        v.append(i)
    ans = 0
    for i in range(n):
        itr = bisect_left(v, arr[i])
        ans += itr
        v = v[:itr] + v[itr + 1 :]
    return ans    

def generator_burau_matrix(generator, index):
    if generator[1] == 1:
        burau22 = Matrix([[
            Symbol(f'a+{index}', commutative=False), 
            Symbol(f'b+{index}', commutative=False)], [
            Symbol(f'c+{index}', commutative=False), 0]])
    elif generator[1] == -1:
        burau22 = Matrix([[0, 
            Symbol(f'c-{index}', commutative=False)], [
            Symbol(f'b-{index}', commutative=False), 
            Symbol(f'a-{index}', commutative=False)]])
    return burau22

def reduced_burau_matrix(braid_sequence, n):
    out = Identity(n + 1)
    for l in range(len(braid_sequence)):
        generator = braid_sequence[l]
        index = int(generator[0])
        factor = BlockDiagMatrix(
            Identity(index - 1), 
            generator_burau_matrix(generator, l + 1), 
            Identity(n - index))
        factor = factor.as_explicit() 
        out = (out * factor)
    return out[1:, 1:]

def walks(q, right_quantum_matrix, n):
    summation = 0
    integers = range(n)
    for m in range(1, n + 1):
        family = itertools.combinations(integers, m)
        shuffles = [(pi, countInversions(pi, m)) for pi in itertools.permutations(range(m))]
        summation_ = 0
        for subset in family:
            subset = list(subset)
            right_quantum_submatrix = right_quantum_matrix[subset, subset]
            for (pi, inversions) in shuffles:
                summand = (-q)**(inversions)
                for index in range(m):
                    summand *= right_quantum_submatrix[pi[index], index]
                summation_ += summand
        summation += (-1)**(m - 1) * summation_
    return summation

def braid_monomial_exponent_counter(q, monomial, NCrs): # works for braid sequences of arbitrary length
    l1 = []; l2 = []; l3 = []; l4 = []; l5 = []
    c1 = [0 for _ in range(NCrs)] 
    c2 = [0 for _ in range(NCrs)]  
    c3 = [0 for _ in range(NCrs)] 
    c4 = [0 for _ in range(NCrs)]
    c5 = [0 for _ in range(NCrs)]
    d = [0 for _ in range(NCrs)]  # because each crossing is either negative or positive, but not both, you can merge the unprimed and primed lists (not yet implemented)
    s = [0 for _ in range(NCrs)] 
    r = [0 for _ in range(NCrs)]
    prefactor = monomial.args[0]
    if prefactor.args:
        coefficient = 1
        if isinstance(prefactor, Pow):
            prefactor = prefactor.args[1]
        else:
            prefactor = 1
    else:
        if prefactor == q:
            coefficient = 1
            prefactor = 1
        else:
            coefficient = prefactor
            prefactor = monomial.args[1]
            if isinstance(prefactor, Pow):
                    prefactor = prefactor.args[1]
            else:
                prefactor = 1
    monomial = str(monomial)
    for j in range(len(monomial)):
        if monomial[j] in ['a', 'b', 'c']:
            index = j - 1
            break
    count = monomial[index + 1 :].count('*')
    for _ in range(count):
        monomial = monomial[index + 1:]
        index = monomial.find('*')
        integer = int(monomial[2:index])
        if monomial[0] == 'a':
            if monomial[1] == '+':
                c1[integer - 1] += 1
                d[integer - 1] += 1
            else:
                c3[integer - 1] += 1
                c4[integer - 1] += 1
                d[integer - 1] += 1
        elif monomial[0] == 'b':
            if monomial[1] == '+':
                s[integer - 1] += 1
                l2.append(c2[integer - 1])
            else:
                s[integer - 1] += 1
                l3.append(c3[integer - 1])
                l5.append(c5[integer - 1])
        else:
            if monomial[1] == '+':
                c2[integer - 1] += 1
                r[integer - 1] += 1
                l1.append(c1[integer - 1])
            else:
                c5[integer - 1] += 1
                r[integer - 1] += 1
                l4.append(c4[integer - 1])
    monomial = monomial[index + 1:]
    index = len(monomial)
    integer = int(monomial[2:index])
    if monomial[0] == 'a':
        if monomial[1] == '+':
            c1[integer - 1] += 1
            d[integer - 1] += 1
        else:
            c3[integer - 1] += 1
            c4[integer - 1] += 1
            d[integer - 1] += 1
    elif monomial[0] == 'b':
        if monomial[1] == '+':
            s[integer - 1] += 1
            l2.append(c2[integer - 1])
        else:
            s[integer - 1] += 1
            l3.append(c3[integer - 1])
            l5.append(c5[integer - 1])
    else:
        if monomial[1] == '+':
            c2[integer - 1] += 1
            r[integer - 1] += 1
            l1.append(c1[integer - 1])
        else:
            c5[integer - 1] += 1
            r[integer - 1] += 1
            l4.append(c4[integer - 1])
    return [sum(l1) - 2*sum(l2) + 2*sum(l3) - sum(l4) + 2*sum(l5), d, s, r, coefficient, prefactor]

def SWC(q, N, NCrs, signs, alt3):
    summation = 0
    SWC.walks = SWC.walks.expand()
    if isinstance(SWC.walks, Add):
        SWC.walks = list(SWC.walks.args)
    else:
        SWC.walks = [SWC.walks]
    for index in range(len(SWC.walks) -1, -1, -1):
        z, d, s, r, coefficient, prefactor = braid_monomial_exponent_counter(q, SWC.walks[index], NCrs)                        # don't get rid of z computation here
        if (d + np.maximum(s, r) >= 2).any():
            SWC.walks.pop(index)
        else:
            summand = q**(prefactor + z - sum(map(mul, map(mul, d, r), alt3)) + sum(map(mul, r, signs)) * (N - 1))                   # the sums over + and - r's go to dot product between signs and r under merge
            for i in range(NCrs):
                for h in range(d[i]):
                    summand *= (1 - q**(signs[i] * (N - 1 - r[i] - h)))
            summation += coefficient * summand
            SWC.walks[index] = [coefficient, z + prefactor, d, s, r]
    return summation

def SWC_mirror(q, N, NCrs, signs, alt3):
    summation = 0
    SWC_mirror.walks = SWC_mirror.walks.expand()
    if isinstance(SWC_mirror.walks, Add):
        SWC_mirror.walks = list(SWC_mirror.walks.args)
    else:
        SWC_mirror.walks = [SWC_mirror.walks]
    for index in range(len(SWC_mirror.walks) -1, -1, -1):
        z, d, s, r, coefficient, prefactor = braid_monomial_exponent_counter(q, SWC_mirror.walks[index], NCrs)                        # don't get rid of z computation here
        if (d + np.maximum(s, r) >= 2).any():
            SWC_mirror.walks.pop(index)
        else:
            summand = q**(prefactor + z - sum(map(mul, map(mul, d, r), alt3)) + sum(map(mul, r, signs)) * (N - 1))                   # the sums over + and - r's go to dot product between signs and r under merge
            for i in range(NCrs):
                for h in range(d[i]):
                    summand *= (1 - q**(signs[i] * (N - 1 - r[i] - h)))
            summation += coefficient * summand
            SWC_mirror.walks[index] = [coefficient, z + prefactor, d, s, r]
    return summation

# merge the above two functions into one function by using different attribute names and "if mirror_optimization ..."

def product(q, N, NCrs, signs, alt1, alt2, alt3):
    summation = 0
    for index in range(len(product.walks) -1, -1, -1):
        d1 = product.walks[index][0][2]; d2 = product.walks[index][1][2]
        s1 = product.walks[index][0][3]; s2 = product.walks[index][1][3]
        r1 = product.walks[index][0][4]; r2 = product.walks[index][1][4]
        d, s, r = np.add([d1, s1, r1], [d2, s2, r2])
        if (d + np.maximum(s, r) >= N).any():
            product.walks.pop(index)
        else:
            z = sum(map(mul, map(mul, d1, r2), signs)) - sum(map(mul, map(mul, r1, s2), alt1)) - sum(map(mul, map(mul, d1, s2), alt2))
            coefficient = product.walks[index][0][0] * product.walks[index][1][0]
            prefactor = product.walks[index][0][1] + product.walks[index][1][1]
            summand = q**(prefactor + z - sum(map(mul, map(mul, d, r), alt3)) + sum(map(mul, r, signs)) * (N - 1))
            for i in range(NCrs):
                for h in range(d[i]):
                    summand *= (1 - q**(signs[i] * (N - 1 - r[i] - h)))
            summation += symengine.sympify(coefficient * summand)
            product.walks[index] = [coefficient, z + prefactor, d, s, r]
    print('oof')
    return summation

def Jones(braid_word, N = 2, mirror_optimization = True):
    q = Symbol('q', commutative=True)
    NCrs = len(braid_word)
    braid_sequence = np.array([np.abs(braid_word), np.sign(braid_word)]).T
    m = max(braid_sequence[:, 0])
    signs = braid_sequence[:, 1]
    alt3 = list(map(lambda x: (x + 1) / 2, signs))
    writhe = sum(signs)
    matrix = reduced_burau_matrix(braid_sequence, m) # reducing to simple walks here ("after each multiplication") is supposed to speed up computation time (not yet implemented). Where are other places we can insert SWC?
    SWC.walks = walks(q, q * matrix, m)
    non_mirror = SWC(q, N, NCrs, signs, alt3)
    mirroring = False
    if mirror_optimization:
        signs *= -1
        alt3 = list(map(lambda x: 1 - x, alt3))
        matrix = reduced_burau_matrix(braid_sequence, m)
        SWC_mirror.walks = walks(q, q * matrix, m)
        mirror = SWC_mirror(q, N, NCrs, signs, alt3)
        if len(SWC_mirror.walks) < len(SWC.walks):
            SWC.walks = SWC_mirror.walks
            writhe *= -1
            mirroring = True
        else:
            signs *= -1
            alt3 = list(map(lambda x: 1 - x, alt3))
    polynomial = 1
    if mirroring:
        polynomial += mirror
    else:
        polynomial += non_mirror
    if not SWC.walks:
        return q**((N - 1) * (writhe - m) / 2) * polynomial.expand()
    product.walks = SWC.walks
    alt1 = list(map(lambda x: 2 * x, signs)); alt2 = list(map(lambda x: x - 1, signs)); 
    while True:
        product.walks = list(itertools.product(SWC.walks, product.walks))
        polynomial += product(q, N, NCrs, signs, alt1, alt2, alt3)
        if not product.walks:
            polynomial = (q**((N - 1) * (writhe - m) / 2) * simplify(polynomial)).expand()
            if mirror_optimization:
                if mirroring:
                    polynomial = sympify(polynomial).subs({q:q**(-1), q**(-1):q}, simultaneous=True)
            return polynomial

def format(polynomial):
    data = polynomial.as_coefficients_dict()
    keys = list(data.keys())
    for index in range(len(keys)):
        if isinstance(keys[index], Pow):
            keys[index] = int(keys[index].args[1])
        else:
            if str(keys[index]) == 'q':
                keys[index] = 1
            else:
                keys[index] = 0
    values = list(data.values())
    index = np.argsort(keys)
    return list(zip(
        [keys[i] for i in index], 
        [values[j] for j in index]
    ))

if __name__ == '__main__':

    braid_word = [1, 1, 1, 2, -1, 2] # -5_2
    # braid_word = [1, 1, 2, -1, 2, 3, -2, -4, 3, -4] # 8_1
    # braid_word = [1, 1, -2, 1, 3, 2, 2, 2, 3] # 8_15
    # braid_word = [1, 1, 1, -2, 1, 1, -2, -2] # 8_10
    #braid_word = [1, -2, 1, -2, 3, -2, 3] # 7_7
    # braid_word = [1, 1, 1, 2, -1, 2, 2, 2, 3, -2, 3] # 9_18
    #braid_word = [1, 1,1]
    braid_word = [1, 1, 1, 1, 1, 2, -1, 2, -3, 2, -3]

    import time
    start_time = time.time()

    print(format(Jones(braid_word, N=4, mirror_optimization=False)))
    
    end_time = time.time()
    elapsed_time = (end_time - start_time)
    print(f"Time Consumption: {elapsed_time} seconds")

    # In this file, I store algebraic expression using combined lists for the + and - Burau variables. This actually seems to slow down the algorithm.
    # For some reason, the polynomials returned by this algorithm are sometimes not simplified (e.g. terms in 1/q and q**(-1) would not be added),
    # so I had to add a sympy.simplify call before the final expression is returned. The algorithm is still slower without the simplify call, though.

    # does turning q to symengine once remove the performance nerfs for small knots?