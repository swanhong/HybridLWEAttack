# -*- coding: utf-8 -*-

from sage.all import shuffle, randint, ceil, next_prime, log, cputime, mean, variance, set_random_seed, sqrt
from copy import copy
from sage.all import GF, ZZ
from sage.all import random_matrix, random_vector, vector, matrix, identity_matrix
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian
from estimator.estimator import preprocess_params, stddevf


def gen_fhe_instance(n, q, alpha=None, h=None, m=None, seed=None):
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

    n, alpha, q = preprocess_params(n, alpha, q)

    stddev = stddevf(alpha*q)

    if m is None:
        m = n
    K = GF(q, proof=False)
    A = random_matrix(K, m, n)

    if h is None:
        s = random_vector(ZZ, n, x=-1, y=1)
    else:
        S = [-1, 1]
        s = [S[randint(0, 1)] for i in range(h)]
        s += [0 for _ in range(n-h)]
        shuffle(s)
        s = vector(ZZ, s)
    c = A*s

    D = DiscreteGaussian(stddev)

    for i in range(m):
        c[i] += D()

    return A, c


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
    Generate dual attack basis for LWE normal form.

    :param A: LWE matrix A

    """
    q = A.base_ring().order()
    n = A.ncols()
    B = A.matrix_from_rows(range(0, n)).inverse().change_ring(ZZ)
    L = identity_matrix(ZZ, n).augment(B)
    L = L.stack(matrix(ZZ, n, n).augment(q*identity_matrix(ZZ, n)))

    for i in range(0, 2*n):
        for j in range(n, 2*n):
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


def dim_error_tradeoff(A, c, beta, h, m=None, tau, scale=1, float_type="double", k):
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

    if m is None:
        m = A.nrows()

    A1 = concatenate(A, n-k)

    L = dual_instance1(A, scale=scale)
    L = IntegerMatrix.from_matrix(L)
    L = LLL.reduction(L, flags=LLL.VERBOSE)
    M = GSO.Mat(L, float_type=float_type)
    bkz = BKZ2(M)
    t = 0.0
    param = BKZ.Param(block_size=beta,
                      strategies=BKZ.DEFAULT_STRATEGY,
                      auto_abort=True,
                      max_loops=16,
                      flags=BKZ.VERBOSE|BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
    bkz(param)
    t += bkz.stats.total_time

    H = copy(L)

    # import pickle
    # pickle.dump(L, open("L-%d-%d.sobj"%(L.nrows, beta), "wb"))

    # E = []
    E = matrix(K, tau, k)
    f = vector(K, tau)
    Y = set()
    # V = set()
    y_i = vector(ZZ, tuple(L[0]))
    Y.add(tuple(y_i))

    temp = apply_short1(y_i, A, c, scale=scale)
    E[0] = temp[0]
    f[0] = temp[1]

    v = L[0].norm()
    v_ = v/sqrt(L.ncols)
    v_r = 3.2*sqrt(L.ncols - A.ncols())*v_/scale
    v_l = sqrt(h)*v_

    # fmt = u"{\"t\": %5.1fs, \"log(sigma)\": %5.1f, \"log(|y|)\": %5.1f, \"log(E[sigma]):\" %5.1f}"

    print
    print fmt%(t,
               log(abs(E[-1]), 2),
               log(L[0].norm(), 2),
               log(sqrt(v_r**2 + v_l**2), 2))
    print
    for i in range(1, tau):
        t = cputime()
        M = GSO.Mat(L, float_type=float_type)
        bkz = BKZ2(M)
        t = cputime()
        bkz.randomize_block(0, L.nrows, stats=None, density=3)
        LLL.reduction(L)
        y_i = vector(ZZ, tuple(L[0]))
        l_n = L[0].norm()
        if L[0].norm() > H[0].norm():
            L = copy(H)
        t = cputime(t)

        Y.add(tuple(y_i))
        # V.add(y_i.norm())
        temp = apply_short1(y_i, A, c, scale=scale)
        E[i] = temp[0]
        f[i] = temp[1]
        #if len(V) >= 2:
        #    fmt =  u"{\"i\": %4d, \"t\": %5.1fs, \"log(|e_i|)\": %5.1f, \"log(|y_i|)\": %5.1f,"
        #    fmt += u"\"log(sigma)\": (%5.1f,%5.1f), \"log(|y|)\": (%5.1f,%5.1f), |Y|: %5d}"
        #    print fmt%(i+2, t, log(abs(E[-1]), 2), log(l_n, 2), log_mean(E), log_var(E), log_mean(V), log_var(V), len(Y))

    return E, f

def generate_table(S, tau):    

    T = {} # empty dictionary

    for vec in S:
        sgnvec = Power(sgn(vec)) # dictionary의 key 부분이 무지막지하게 (2^tau 수준) 커져도 괜찮나?
        if sgnvec in T:
            T[sgnvec].append(vec)
        else:
            T[sgnvec] = []
            T[sgnvec].append(vec)

    return T

def noisy_search(q, T, bound, query):
        
    tau = len(query)
    sgn_ = vector(ZZ, tau)
    index = []
    for i in range(tau):
        if bound < abs(query[i]) < q//2 - bound:
            sgn_[i] = query[i] // (q//2) # 0 if query[i] < 0, 1 otherwise
        else:
            sgn_[i] = 1
            index.append(i)

    return check_collision(sgn_, T, index, 0)

def check_collision(v, T, index, index_num):
    '''

    :Param index:       A set of positions marked by 1 ( x in paper )
    :Param index_num:   Current position where check_collision is

    '''

    if index_num == len(index): 
        if v in T:
            for vec in T[v]:
                if max_norm(query - vec) <= bound:
                    return vec
        return None

    else:
        for i in [0, 1]:
            v[index[index_num]] = i
            vec = check_collision(v, T, index, index_num + 1)
            if vec is not None:
                return vec
    return None
            

def hybrid_mitm(A, c, beta, h, m=None, tau, scale=1, float_type="double", k, ell):
    
    # Lattice reduction stage
    E, f = dim_error_tradeoff(A, c, beta, h, m, tau, scale, float_type, k)

    # MITM stage
    S = data_gen(E, ell)

    T = generate_table(S, tau) # T : binary table sgn(S)

    for query in S:
        res = noisy_search(q, T, bound, f - query)
        if res != None:
            return True

    return False

def data_gen(A, ell):
    k = A.ncols()
    s = vector(ZZ, k)
    s[0] = 1
    S = set()
    S.add = A[0]
    S_recursive(A, S, s, 1, ell, 1)
    return S

def data_recursive(A, S, s, position, ell, num_one):
    # FIXME
    if position < len(s) and num_one <= ell :
        s[position] = 0
        S_recursive(A, S, s, position + 1, ell)

        s[position] = 1
        num_one += 1
        S.add(A[position])
        for vec in S:
            S.add(vec + A[position])
        S_recursive(A, S, s, position + 1, ell + 1)        

def concatenate(A, k):
    m = A.nrows()
    A1 = matrix(ZZ, m, k)
    for i in range(k):
        A1[i] = A[i]
    return A1

def max_norm(vec):
    max = vec[0]
    for i in range(len(1, vec)):
        if max < vec[i]:
            max = vec[i]
    return max