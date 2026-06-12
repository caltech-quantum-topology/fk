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

def braid_monomial_exponent_counter(q, monomial, NCrs): 
    l1 = []; l2 = []; l3 = []; l4 = []; l5 = []
    c1 = [0 for _ in range(NCrs)] 
    c2 = [0 for _ in range(NCrs)]  
    c3 = [0 for _ in range(NCrs)] 
    c4 = [0 for _ in range(NCrs)]
    c5 = [0 for _ in range(NCrs)]
    d = [0 for _ in range(NCrs)]  
    d_ = [0 for _ in range(NCrs)] 
    s = [0 for _ in range(NCrs)] 
    s_ = [0 for _ in range(NCrs)]
    r = [0 for _ in range(NCrs)]
    r_ = [0 for _ in range(NCrs)]
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
                d_[integer - 1] += 1
        elif monomial[0] == 'b':
            if monomial[1] == '+':
                s[integer - 1] += 1
                l2.append(c2[integer - 1])
            else:
                s_[integer - 1] += 1
                l3.append(c3[integer - 1])
                l5.append(c5[integer - 1])
        else:
            if monomial[1] == '+':
                c2[integer - 1] += 1
                r[integer - 1] += 1
                l1.append(c1[integer - 1])
            else:
                c5[integer - 1] += 1
                r_[integer - 1] += 1
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
            d_[integer - 1] += 1
    elif monomial[0] == 'b':
        if monomial[1] == '+':
            s[integer - 1] += 1
            l2.append(c2[integer - 1])
        else:
            s_[integer - 1] += 1
            l3.append(c3[integer - 1])
            l5.append(c5[integer - 1])
    else:
        if monomial[1] == '+':
            c2[integer - 1] += 1
            r[integer - 1] += 1
            l1.append(c1[integer - 1])
        else:
            c5[integer - 1] += 1
            r_[integer - 1] += 1
            l4.append(c4[integer - 1])
    return [sum(l1) - 2*sum(l2) + 2*sum(l3) - sum(l4) + 2*sum(l5), d, d_, s, s_, r, r_, coefficient, prefactor]

def highest_nonzero(array1, array2, array3, array4, array5, array6):
    done1 = False; done2 = False
    in1 = 0; in2 = 0                                                                          
    for index in range(len(array1) - 1, -1, -1):
        if not done1: 
            if array1[index] != 0:
                in1 = index
                done1 = True
            if array2[index] != 0:
                in1 = index
                done1 = True
            if array3[index] != 0:
                in1 = index
                done1 = True
        if not done2:
            if array4[index] != 0:
                in2 = index
                done2 = True
            if array5[index] != 0:
                in2 = index
                done2 = True
            if array6[index] != 0:
                in2 = index
                done2 = True
    return [in1 + 1, in2 + 1]

def SWC(q, N, NCrs):
    summation = 0
    SWC.walks = SWC.walks.expand()
    if isinstance(SWC.walks, Add):
        SWC.walks = list(SWC.walks.args)
    else:
        SWC.walks = [SWC.walks]
    for index in range(len(SWC.walks) -1, -1, -1):
        z, d, d_, s, s_, r, r_, coefficient, prefactor = braid_monomial_exponent_counter(q, SWC.walks[index], NCrs)                    
        if (d + np.maximum(s, r) >= 2).any() or (d_ + np.maximum(s_, r_) >= 2).any():
            SWC.walks.pop(index)
        else:
            u, v = highest_nonzero(d, s, r, d_, s_, r_)
            summand = q**(prefactor + z - sum(map(mul, d, r)) + (sum(r) - sum(r_)) * (N - 1))
            for i in range(u):
                for h in range(0, d[i]):
                    summand *= (1 - q**(N - 1 - r[i] - h))
            for j in range(v):
                for l in range(0, d_[j]):
                    summand *= (1 - q**(r_[j] + l + 1 - N))
            summation += coefficient * summand
            SWC.walks[index] = [coefficient, z + prefactor, d, d_, s, s_, r, r_]
    return summation

def SWC_mirror(q, N, NCrs):
    summation = 0
    SWC_mirror.walks = SWC_mirror.walks.expand()
    if isinstance(SWC_mirror.walks, Add):
        SWC_mirror.walks = list(SWC_mirror.walks.args)
    else:
        SWC_mirror.walks = [SWC_mirror.walks]
    for index in range(len(SWC_mirror.walks) -1, -1, -1):
        z, d, d_, s, s_, r, r_, coefficient, prefactor = braid_monomial_exponent_counter(q, SWC_mirror.walks[index], NCrs)                    
        if (d + np.maximum(s, r) >= 2).any() or (d_ + np.maximum(s_, r_) >= 2).any():
            SWC_mirror.walks.pop(index)
        else:
            u, v = highest_nonzero(d, s, r, d_, s_, r_)
            summand = q**(prefactor + z - sum(map(mul, d, r)) + (sum(r) - sum(r_)) * (N - 1))
            for i in range(u):
                for h in range(0, d[i]):
                    summand *= (1 - q**(N - 1 - r[i] - h))
            for j in range(v):
                for l in range(0, d_[j]):
                    summand *= (1 - q**(r_[j] + l + 1 - N))
            summation += coefficient * summand
            SWC_mirror.walks[index] = [coefficient, z + prefactor, d, d_, s, s_, r, r_]
    return summation

'''def product(q, N, NCrs):
    summation = 0
    for index in range(len(product.walks) -1, -1, -1):
        d, d_, s, s_, r, r_ = [
            list(map(add, product.walks[index][0][2], product.walks[index][1][2])), 
            list(map(add, product.walks[index][0][3], product.walks[index][1][3])), 
            list(map(add, product.walks[index][0][4], product.walks[index][1][4])), 
            list(map(add, product.walks[index][0][5], product.walks[index][1][5])), 
            list(map(add, product.walks[index][0][6], product.walks[index][1][6])), 
            list(map(add, product.walks[index][0][7], product.walks[index][1][7]))
        ]
        if (d + np.maximum(s, r) >= N).any() or (d_ + np.maximum(s_, r_) >= N).any():
            product.walks.pop(index)
        else:
            z = sum(map(mul, product.walks[index][0][2], product.walks[index][1][6])) - 2 * sum(map(mul, product.walks[index][0][6], product.walks[index][1][4])) + 2 * sum(map(mul, map(add, product.walks[index][0][3], product.walks[index][0][7]), product.walks[index][1][5])) - sum(map(mul, product.walks[index][0][3], product.walks[index][1][7]))
            coefficient = product.walks[index][0][0] * product.walks[index][1][0]
            prefactor = product.walks[index][0][1] + product.walks[index][1][1]
            u, v = highest_nonzero(d, s, r, d_, s_, r_)
            summand = q**(prefactor + z - np.dot(d, r) + (sum(r) - sum(r_)) * (N - 1))
            for i in range(u):
                for h in range(0, d[i]):
                    summand *= (1 - q**(N - 1 - r[i] - h))
            for j in range(v):
                for l in range(0, d_[j]):
                    summand *= (1 - q**(r_[j] + l + 1 - N))
            summation += coefficient * summand
            product.walks[index] = [coefficient, z + prefactor, d, d_, s, s_, r, r_]
    return summation'''

def f(walks, q, N, NCrs):
    a = []; b = []; c = []; e = []; f = []; g = []
    z = 0
    summand = 0
    mult = 1
    execute1 = False
    execute2 = False
    for j in range(NCrs):
        d = walks[0][2][j] + walks[1][2][j]
        d_ = walks[0][3][j] + walks[1][3][j]
        s = walks[0][4][j] + walks[1][4][j]
        s_ = walks[0][5][j] + walks[1][5][j]
        r = walks[0][6][j] + walks[1][6][j]
        r_ = walks[0][7][j] + walks[1][7][j]
        a.append(d); b.append(d_); c.append(s); e.append(s_); f.append(r); g.append(r_)
        if not execute1:
            if d != 0 or s != 0 or r != 0:
                execute1 = True
        if execute1:
            if d + max(s, r) >= N:
                return None
            for h in range(0, d):
                mult *= (1 - q**(N - 1 - r - h))
        if not execute2:
            if d_ != 0 or s_ != 0 or r_ != 0:
                execute2 = True
        if execute2:
            if d_ + max(s_, r_) >= N:
                return None
            for l in range(0, d_):
                mult *= (1 - q**(r_ + l + 1 - N))
        z += walks[0][2][j] * walks[1][6][j] - 2 * walks[0][6][j] * walks[1][4][j] + 2 * (walks[0][3][j] + walks[0][7][j]) * walks[1][5][j] - walks[0][3][j] * walks[1][7][j]
        summand += (r - r_) * (N - 1) - d * r
    coefficient = walks[0][0] * walks[1][0]
    prefactor = walks[0][1] + walks[1][1]
    summand += prefactor + z
    summand = q**(summand)
    summand *= mult
    new_term = [coefficient, z + prefactor, a, b, c, e, f, g]
    return [symengine.sympify(coefficient * summand), new_term]

import multiprocessing
def product(q, N, NCrs):
    print('pre-summation!')
    results = []
    pool = multiprocessing.Pool()
    for index in range(len(product.walks) -1, -1, -1):
        # result = pool.apply_async(f, args=(product.walks[index], q, N, NCrs))
        result = f(product.walks[index], q, N, NCrs)
        results.append(result)
    pool.close()
    pool.join()
    print('made it to summation!')
    summation = 0
    print(len(product.walks))
    for index in range(len(product.walks) -1, -1, -1):
        # result = results[index].get()
        result = results[index]
        if result == None:
                product.walks.pop(index)
        else:
            summation += result[0].expand()
            product.walks[index] = result[1]
    return summation

def Jones(braid_word, N = 2, mirror_optimization = True):
    q = Symbol('q', commutative=True)
    NCrs = len(braid_word)
    braid_sequence = np.array([np.abs(braid_word), np.sign(braid_word)]).T
    m = max(braid_sequence[:, 0])
    writhe = sum(braid_sequence[:, 1])
    matrix = reduced_burau_matrix(braid_sequence, m)
    SWC.walks = walks(q, q * matrix, m)
    non_mirror = SWC(q, N, NCrs)
    mirroring = False
    if mirror_optimization:
        braid_sequence[:, 1] *= -1
        matrix = reduced_burau_matrix(braid_sequence, m)
        SWC_mirror.walks = walks(q, q * matrix, m)
        mirror = SWC_mirror(q, N, NCrs)
        if len(SWC_mirror.walks) < len(SWC.walks):
            SWC.walks = SWC_mirror.walks
            writhe *= -1
            mirroring = True
    polynomial = 1
    if mirroring:
        polynomial += mirror
    else:
        polynomial += non_mirror
    q = symengine.sympify(q)
    if not SWC.walks:
        return q**((N - 1) * (writhe - m) / 2) * polynomial.expand()
    product.walks = SWC.walks
    while True:
        product.walks = list(itertools.product(SWC.walks, product.walks))
        new = time.time(); current = new 
        polynomial += product(q, N, NCrs) # product function is main time consumer
        new = time.time(); print("addition of product to polynomial:", new - current); current = new 
        if not product.walks:
            polynomial = (q**((N - 1) * (writhe - m) / 2) * polynomial).expand()
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

    braid_word = [1, 1, 1, 1, 1, 2, -1, 2, -3, 2, -3]
    braid_word = [1, 1, 1, 2, -1, 2]

    import time
    start_time = time.time() 
    print((Jones(braid_word, N=4, mirror_optimization=True)))
    
    end_time = time.time()
    elapsed_time = (end_time - start_time)
    print(f"Time Consumption: {elapsed_time} seconds")