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


def dual_instance0(A):
    """
    Generate dual attack basis.

    :param A: LWE matrix A

    """
    q = A.base_ring().order()
    B0 = A.left_kernel().basis_matrix().change_ring(ZZ)
    m = B0.ncols()
    n = B0.nrows()
    r = m-n
    B1 = matrix(ZZ, r, n).augment(q*identity_matrix(ZZ, r))
    B = B0.stack(B1)
    return B


def dual_instance1(A, scale=1):
    """
    FIXME : m != n 인 경우에 해결하기

    Generate dual attack basis for LWE normal form.

    :param A: LWE matrix A
    """
    q = A.base_ring().order()
    n = A.ncols()
    m = A.nrows()
    B = A.matrix_from_rows(range(m-n, m)).inverse().change_ring(ZZ)
    L = identity_matrix(ZZ, n).augment(matrix(ZZ,n,m-n))
    L = L.augment(B)

    L = L.stack(matrix(ZZ, m-n, n).augment(A.left_kernel().basis_matrix().change_ring(ZZ)))
    

    L = L.stack(matrix(ZZ, n, m).augment(q*identity_matrix(ZZ, n)))

    for i in range(0, m + n):
        for j in range(n, m + n):
            L[i, j] = scale*L[i, j]

    return L


def balanced_lift(e):
    """
    Lift e mod q to integer such that result is between -q/2 and q/2

    :param e: a value or vector mod q

    """
    from sage.rings.finite_rings.integer_mod import is_IntegerMod

    q = e.base_ring().order()
    if is_IntegerMod(e):
        e = ZZ(e)
        if e > q//2:
            e -= q
        return e
    else:
        return vector(balanced_lift(ee) for ee in e)


def apply_short1(y, A, c, scale=1):
    """
    Compute `y*A`, `y*c` where y is a vector in the integer row span of
    ``dual_instance(A)``

    :param y: (short) vector in scaled dual lattice
    :param A: LWE matrix
    :param c: LWE vector
    """
    m = A.nrows()
    y = vector(ZZ, 1/ZZ(scale) * y[-m:])
    a = balanced_lift(y*A)
    e = balanced_lift(y*c)
    return a, e


def log_mean(X):
    return log(mean([abs(x) for x in X]), 2)


def log_var(X):
    return log(variance(X).sqrt(), 2)


def dim_error_tradeoff(A, c, beta, h, k, s, num_sample = None, float_type="double"):
    """

    :param A:    LWE matrix
    :param c:    LWE vector
    :param beta: BKW block size
    :param m:    number of given LWE samples
    :param tau:  number of new samples to generate
    :param scale: scale rhs of lattice by this factor
    :param k:    LWE dim after tradeoff

    """
    from fpylll import BKZ, IntegerMatrix, LLL, GSO
    from fpylll.algorithms.bkz2 import BKZReduction as BKZ2

    n = A.ncols()
    q = A.base_ring().order()
    K = GF(q, proof=False)

    scale = round(8 * sqrt(n) / sqrt(2*pi*h))
    scale = ZZ(scale)

    A1 = A.matrix_from_columns([i for i in range(n - k)])
    
    L = dual_instance1(A1, scale=scale)
    L = IntegerMatrix.from_matrix(L)
    L = LLL.reduction(L, flags=LLL.VERBOSE)
    M = GSO.Mat(L, float_type=float_type)
    bkz = BKZ2(M)
    param = BKZ.Param(block_size=beta,
                      strategies=BKZ.DEFAULT_STRATEGY,
                      auto_abort=True,
                      max_loops=16,
                      flags=BKZ.VERBOSE|BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
    bkz(param)
    
    H = copy(L)
    A = A.matrix_from_columns([n - k + i for i in range(k)])
	
    E = matrix(ZZ, 1, k)
    f = []
    V = set()
    y_i = vector(ZZ, tuple(L[0]))
    E[0], ft = apply_short1(y_i, A, c, scale=scale)
    f.append(ft)
    	
    if num_sample is None:
    	num_sample = 15

    while len(V) < num_sample:
        M = GSO.Mat(L, float_type=float_type)
        bkz = BKZ2(M)
        bkz.randomize_block(0, L.nrows, density=5)
        LLL.reduction(L)
        y_i = vector(ZZ, tuple(L[0]))
        l_n = L[0].norm()
        if L[0].norm() > H[0].norm():
            L = copy(H)

        count = len(V)
        V.add(y_i.norm())
        if count != len(V):
        	print '|V| : ', len(V)
    		Etmp, ftmp = apply_short1(y_i, A, c, scale=scale)
        	E = E.stack(Etmp)
        	f.append(ftmp)

    E = E.change_ring(K)
    f = vector(K, f)
    e = vector(f - E * s, ZZ)
    e = balanced_lift(e)
    B = max_norm(e)
    print '|E| = ', E.nrows(), 'Bound = ', B

    return E, f, B

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
            

def hybrid_mitm(A, c, beta, h, k, s, num_sample = None, ell = None, m=None, float_type="double"):
    # Lattice reduction stage
    E, f, bound = dim_error_tradeoff(A, c, beta, h, k, s, num_sample = num_sample)

    if ell is None:
    	ell = h

	q = A.base_ring().order()
	
    # MITM stage
    Mitm_on_LWE(E, f, bound, ell, q)
    

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

def Mitm_on_LWE(A, c, bound, h, q):

    if h % 2 == 0:
        halfh = h // 2
    else:
        halfh = h // 2 + 1

    S = data_gen(A, halfh) # S = {As : Hw(s) <= ell}
    T = generate_table(S, q)
    is_LWE = False

    for v in S:
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
            #s = solveLA(A, c, query, res, q)
            #s = vector(ZZ, s)
            #if max_norm(s) <= 1:
            break

    
    if is_LWE == True:
    	print 'Input is from LWE samples'
        #print 'Input is from LWE samples with secret s:'
        #print s
    else:
        print 'Input is NOT from LWE'

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
        m = n
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

    return A, c, u, s
    

    