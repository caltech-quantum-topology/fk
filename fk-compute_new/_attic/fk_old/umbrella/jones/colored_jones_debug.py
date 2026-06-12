import numpy as np
from sympy import *
from itertools import permutations, combinations
from bisect import bisect_left

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

def walks(right_quantum_matrix, n):
    summation = 0
    integers = range(n)
    q = Symbol('q', commutative=True)
    for m in range(1, n + 1):
        family = combinations(integers, m)
        shuffles = [(pi, countInversions(pi, m)) for pi in permutations(range(m))]
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

def braid_monomial_exponent_counter(monomial, NCrs): 
    monomial = str(monomial)
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
    return [sum(l1) - 2*sum(l2) + 2*sum(l3) - sum(l4) + 2*sum(l5), d, d_, s, s_, r, r_]

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

def monomial_evaluation(monomial, N, NCrs):
    q = Symbol('q', commutative=True)
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
    z, d, d_, s, s_, r, r_ = braid_monomial_exponent_counter(monomial, NCrs)
    if (d + np.maximum(s, r) >= N).any() or (d_ + np.maximum(s_, r_) >= N).any(): 
       return [0, True]
    u, v = highest_nonzero(d, s, r, d_, s_, r_)
    out = q**(prefactor + z - np.dot(d, r) + (sum(r) - sum(r_)) * (N - 1))
    for index in range(u):
        for h in range(0, d[index]):
            out *= (1 - q**(N - 1 - r[index] - h))
    for index in range(v):
        for l in range(0, d_[index]):
            out *= (1 - q**(r_[index] + l + 1 - N))
    return [coefficient * out, False]

def evaluation(N, NCrs):
    summation = 0
    for index in range(len(evaluation.terms) - 1, -1, -1):
        score, delete = monomial_evaluation(evaluation.terms[index], N, NCrs)
        if delete:
            evaluation.terms.pop(index)
        summation += score
    return summation

def SWC(walks, NCrs):
    walks = walks.expand()
    if isinstance(walks, Add):
        walks = list(walks.args)
    else:
        walks = [walks]
    for index in range(len(walks) -1, -1, -1):
        z, d, d_, s, s_, r, r_ = braid_monomial_exponent_counter(walks[index], NCrs)                     
        if (d + np.maximum(s, r) >= 2).any() or (d_ + np.maximum(s_, r_) >= 2).any():
            walks.pop(index)
    return sum(walks)

def Jones(braid_word, N = 2, mirror_optimization = True):
    NCrs = len(braid_word)
    braid_sequence = np.array([np.abs(braid_word), np.sign(braid_word)]).T
    m = max(braid_sequence[:, 0])
    writhe = sum(braid_sequence[:, 1])
    matrix = reduced_burau_matrix(braid_sequence, m) 
    q = Symbol('q', commutative=True)
    C = walks(q * matrix, m)
    C = SWC(C, NCrs)
    if mirror_optimization:
        mirroring = False
        braid_sequence[:, 1] *= -1
        matrix = reduced_burau_matrix(braid_sequence, m)
        C_mirror = walks(q * matrix, m)
        C_mirror = SWC(C_mirror, NCrs)
        if len(C_mirror.expand().args) < len(C.expand().args):
            C = C_mirror
            writhe *= -1
            mirroring = True
    polynomial = 1
    poly = expand(C)


    if isinstance(poly, Add):
        evaluation.terms = list(poly.args)
    else:
        evaluation.terms = [poly]
    polynomial += evaluation(N, NCrs)
    C = sum(evaluation.terms)
    poly = sum(evaluation.terms)


    if C == 0:
        return q**((N - 1) * (writhe - m) / 2) * polynomial.expand()

    depth = 1
    while True:
        poly = (poly * C).expand() 
        if isinstance(poly, Add):
            evaluation.terms = list(poly.args)
        else:
            evaluation.terms = [poly]
        polynomial += evaluation(N, NCrs)
        poly = sum(evaluation.terms)
        if poly == 0:
            polynomial = (q**((N - 1) * (writhe - m) / 2) * polynomial).expand()
            if mirror_optimization:
                if mirroring:
                    polynomial = polynomial.subs({q:q**(-1), q**(-1):q}, simultaneous=True)
            return polynomial
        depth += 1

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
    braid_word = [1, 1, 2, -1, 2, 3, -2, -4, 3, -4] # 8_1
    braid_word = [1, 1, -2, 1, 3, 2, 2, 2, 3] # 8_15
    braid_word = [1, 1, 1, -2, 1, 1, -2, -2] # 8_10
    #braid_word = [1, -2, 1, -2, 3, -2, 3] # 7_7
    braid_word = [-1, -1, -1]
    #braid_word = list(map(lambda x: -x, braid_word))

    import time
    start_time = time.time()

    # mirror optimization giving different results for some knots; check this (tried on 8_10)
    print(format(Jones(braid_word, N=3, mirror_optimization=True)))
    
    end_time = time.time()
    elapsed_time = (end_time - start_time)
    print(f"Time Consumption: {elapsed_time} seconds")

    '''
    Job 1. I still have to find out if the formula presented works for links.
    For the Hopf Link I get -J_2(Hopf). This is not the same as merely changing orientation,
    because the sign change is an overall coefficient and not in the exponent.
    Also, Jones([1, -1]) is equal to q^(-1/2) + q^(1/2), and not -q^2 - q^(-2), as would
    be expected by computing the Kaufmann bracket of two disjoint unknots.
    If not, what is the generalization of the walk-computational method to links,
    if any such (efficient) method could be found?

    Job 2. The run time of this algorithm greatly decreases as the length of the 
    input braid word decreases. Consider producing a database of minimal braid 
    sequences.

    Job 3. It would be great to have the outputs of the "walks" function pre-computed
    for knots we are interested in. Then, the "Jones" function could be modified
    to take these as input, and the polynomials for any colors for that knot could
    be quickly produced. If this is complete, the "mirroring" variable from above
    could be stored with the walks and we only need the shortest of the two walks
    between that of the knot and its mirror.
    '''