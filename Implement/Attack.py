# -*- coding: utf-8 -*-

from sage.all import shuffle, randint, ceil, next_prime, log, cputime, mean, variance, set_random_seed, sqrt
from copy import copy
from sage.all import GF, ZZ, RR, RealField
from sage.all import random_matrix, random_vector, vector, matrix, identity_matrix, transpose
from sage.structure.element import parent
from sage.symbolic.all import pi
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian
import sys

def generate_table(S, q):    

    T = {} # empty dictionary

    for v in S:
        sgnvec = power(sgn(v, q))
        if sgnvec in T:
            T[sgnvec].append(v)
        else:
            T[sgnvec] = []
            T[sgnvec].append(v)

    return T

def noisy_search(q, T, bound, query):

    sgn_, index = sgnvar(query, bound, q)

    print 'Current number of x = ', len(index)

    v = check_collision(query, sgn_, T, index, 0, bound)

    sys.stdout.write("\033[F")
    sys.stdout.write("\033[K")

    return v


def sgn(v, q):
    length = len(v)
    bin = []
    for i in range(length):
        if v[i] >= q//2:
            bin.append(1)
        else:
            bin.append(0)
    return bin

def sgnvar(v, bound, q):
    
    length = len(v)
    bin = []
    index = []
    
    for i in range(length):
        tmp = abs(ZZ(v[i]) - q//2)
        if bound < tmp < q//2 - bound:
            if ZZ(v[i]) >= q//2:
                bin.append(1)
            else:
                bin.append(0)
        else:
            bin.append(-1) # dummy sign
            index.append(i)

    return bin, index

def check_collision(query, sgn_, T, index, index_num, bound):
    '''

    :Param index:       A set of positions marked by 1 ( x in paper )
    :Param index_num:   Current position where check_collision is

    '''

    if index_num == len(index): 
        if power(sgn_) in T:
            for vec in T[power(sgn_)]:
                vec = vector(ZZ, vec)
                query = vector(ZZ, query)
                if max_norm(query - vec) <= bound:
                    return vec

    else:
        for i in [0, 1]:
            sgn_[index[index_num]] = i
            res = check_collision(query, sgn_, T, index, index_num + 1, bound)
            if res is not None:
                return res
            

'''def hybrid_mitm(A, c, beta, h, tau, k, ell, m=None, scale=1, float_type="double"):
    
    # Lattice reduction stage
    E, f = dim_error_tradeoff(A, c, beta, h, m, tau, scale, float_type, k)

    # MITM stage
    S = data_gen(E, ell)

    T = generate_table(S, q) # T : binary table sgn(S)

    for query in S:
        res = noisy_search(q, T, bound, f - query)
        if res != None:
            return True

    return False'''

def data_gen(A, ell):

    k = A.ncols()
    m = A.nrows()    
    A = transpose(A)
    S = set()
    Count = vector(ZZ, [1] * k)
    A = A.augment(Count)
    S.add(tuple([0] * (m+1)))

    for i in range(k):
        TMP = set()
        TMP.add(tuple(A[i]))   
        for v in S:
            if v[-1] < ell:
                tmp = vector(ZZ, v)
                TMP.add(tuple(tmp + A[i]))
        S = S.union(TMP)

    S_ = set()
    for v in S:
        tmp = vector(ZZ, v)
        tmp = tmp[:m]
        S_.add(tuple(tmp))
    return S_

'''def data_recursive(A, S, s, position, ell, num_one):
    # FIXME
    if position < len(s) and num_one <= ell :
        s[position] = 0
        data_recursive(A, S, s, position + 1, ell)

        s[position] = 1
        num_one += 1
        S.add(A[position])
        for v in S:
            S.add(v + A[position])
        data_recursive(A, S, s, position + 1, ell + 1)   '''     

def max_norm(v):
    length = len(v)
    max = abs(v[0])
    for i in range(1, length):
        if max < abs(v[i]):
            max = abs(v[i])
    return max

def power(v):
    length = len(v)
    pow = v[0]
    for i in range(1, length):
        pow += 2**i * v[i]
    return pow

def fhe_experiment(A, c, bound, h, q):

    if h % 2 == 0:
        halfh = h // 2
    else:
        halfh = h // 2 + 1

    S = data_gen(A, halfh) # S = {As : Hw(s) <= ell}
    T = generate_table(S, q)
    is_LWE = False
    count = 0
    for v in S:
        count += 1
        print 'Working...'

        v = vector(ZZ, v)
        c = vector(ZZ, c)
        query = c - v
        
        for i in range(len(query)):
            if query[i] < 0: query[i] += q
        res = noisy_search(q, T, bound, query)
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K")

        if res is not None:
            is_LWE = True
            #print 'query            :', query
            #print 'collision in S   :', res
            #print
            s = solveLA(A, c, query, res, q)
            s = vector(ZZ, s)
            if max_norm(s) <= 1:
                break

    
    if is_LWE == True:
        print 'Input is LWE samples with secret s:'
        print s
    else:
        print 'Input is NOT LWE samples'

def solveLA(A, c, query, res, q):
    n = A.rank()
    error = query - res
    for i in range(len(error)):
        if error[i] < 0: error[i] += q

    c -= error
    for i in range(len(c)):
        if c[i] < 0: c[i] += q
    c = vector(GF(q, proof=False), c)

    A_ = A.matrix_from_rows([i for i in range(n)])
    if n != A_.rank():
        print 'A_ is not full-rank'
    else:
        c = c[:n]
        return A_.inverse() * c

def gen_instance(n, q, h, alpha=None, m=None, seed=None):
    """
    Generate FHE-style LWE instance

    :param n:     dimension
    :param q:     modulus
    :param alpha: noise rate (default: 8/q)
    :param h:     hamming weight of the secret (default: 2/3n)
    :param m:     number of samples (default: n)

    """
    if seed is not None:
        set_random_seed(seed)

    q = next_prime(ceil(q)-1, proof=False)
    if alpha is None:
        alpha = ZZ(8)/q
        stddev = 3.2

    #RR = parent(alpha*q)
    #stddev = alpha*q/RR(sqrt(2*pi))

    if m is None:
        m = 2*n
    K = GF(q, proof=False)

    while 1:
        A = random_matrix(K, m, n)
        if A.rank() == n:
            break


    if h is None:
        s = random_vector(ZZ, n, x=-1, y=1)
    else:
        s = [1 for i in range(h)]
        s += [0 for _ in range(n-h)]
        shuffle(s)
        s = vector(ZZ, s)
    c = A*s

    D = DiscreteGaussian(stddev)

    for i in range(m):
        c[i] += D()

    u = random_vector(K, m)
    u = vector(ZZ, u)

    return A, c, u, s
