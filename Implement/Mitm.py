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


def gen_instance(n, q, h, alpha=None, m=None, seed=None, s=None):
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
        stddev = 3.2
    else:
        stddev = alpha * q / sqrt(2*pi)

    #RR = parent(alpha*q)
    #stddev = alpha*q/RR(sqrt(2*pi))

    if m is None:
        m = n
    K = GF(q, proof=False)

    while 1:
        A = random_matrix(K, m, n)
        if A.rank() == n:
            break

    if s is not None:
        c = A*s

    else:
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

    u = random_vector(K, m)

    print '(A, c) is n-dim LWE samples (with secret s) / (A, u) is uniform samples'

    return A, c, u, s


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


def apply_short1(y, A, c, u = None, scale=1):
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

    if u is not None:
    	e2 = balanced_lift(y*u)
    	return a, e, e2

    return a, e


def log_mean(X):
    return log(mean([abs(x) for x in X]), 2)


def log_stddev(X):
    return log(variance(X).sqrt(), 2)


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
        pow += 2 ** i * v[i]
    return pow


def sgn(v, q):

	# v : ZZ vector represented in (0, q)

    length = len(v)
    bin = []
    for i in range(length):
        if v[i] >= 0:
            bin.append(1)
        else:
            bin.append(0)
    return bin

def sgnvar(v, bound, q):

	# v : ZZ vector represented in (-q/2, q/2)

    length = len(v)
    bin = []
    index = []
    
    for i in range(length):
        tmp = abs(v[i])
        if bound < tmp < q//2 - bound:
            if ZZ(v[i]) >= 0:
                bin.append(1)
            else:
                bin.append(0)
        else:
            bin.append(-1) # dummy sign
            index.append(i)

    return bin, index


def dim_error_tradeoff(A, c, u, beta, h, k, alpha = None, tau = None, float_type="mpfr", use_lll = True):
    """

    :param A:    LWE matrix
    :param c:    LWE vector
    :param u:	 Uniform vector
    :param beta: BKW block size
    :param h: 	 Hamming weight of secret
    :param k:    LWE dim after tradeoff
    :param tau:  number of new samples to generate
    :param use_lll: 	If True, run BKZ only once and then run LLL
    					If False, run BKZ iteratively

    * secret vector s is used to see the error term of new LWE(-like) samples.

    """

    from fpylll import BKZ, IntegerMatrix, LLL, GSO
    from fpylll.algorithms.bkz2 import BKZReduction as BKZ2

    n = A.ncols()
    q = A.base_ring().order()
    K = GF(q, proof=False)

    if alpha is None:
    	alpha = 8/q

    if tau is None:
    	tau = 30

    m = A.nrows() / n

    scale = round(alpha * q * sqrt(m) / sqrt(2 * pi * h))
    scale = ZZ(scale)

    count = 0

    A_k = matrix(ZZ, 1, k)
    c_k = []
    u_k = []
    length = 0

    while count < tau:

        r = count * m
        T = A.matrix_from_rows([i + r for i in range(m)])
        ct = c[r : r + m]
        ut = u[r : r + m]

        T1 = T.matrix_from_columns([i for i in range(n - k)])
    
        L = dual_instance1(T1, scale=scale)
        L = IntegerMatrix.from_matrix(L)
        L = LLL.reduction(L)
        M = GSO.Mat(L, float_type=float_type)
        bkz = BKZ2(M)
        param = BKZ.Param(block_size=beta,
                      	strategies=BKZ.DEFAULT_STRATEGY,
                      	auto_abort=True,
                      	max_loops=16,
                      	flags=BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
        bkz(param)

        H = copy(L)
    
        y = vector(ZZ, tuple(L[0]))
        length += y.norm()

        T2 = T.matrix_from_columns([n - k + i for i in range(k)])
        
        A_kt, c_kt, u_kt = apply_short1(y, T2, ct, ut, scale=scale)
        if r == 0:
        	A_k[0] = A_kt
        else:
        	A_k = A_k.stack(A_kt)
        c_k.append(c_kt)
        u_k.append(u_kt)

        count += 1

    length = float(length / tau)
    A_k = A_k.change_ring(K)    
    c_k = vector(K, c_k)
    u_k = vector(K, u_k)

    B = float(2 + 1/sqrt(2*pi)) * (alpha * q)
    B = B * B * m / (m + n)
    B = sqrt(B) * length / scale

    print '(A_k, c_k) is k-dim LWE samples (with secret s[-k:]) / (A_k, u_k) is uniform samples. '

    return A_k, c_k, u_k, B


def generate_table(S, q):   

	# Generate (binary) hash table T induced from locality hash function 'sgn'

    #print 'Generating table ... '

    T = {}
    
    for v in S:
        v = vector(ZZ, v)
        for i in range(len(v)):
            if v[i] > q//2: v[i] -= q
        sgnvec = power(sgn(v, q))
        if sgnvec in T:
            T[sgnvec].append(v)
        else:
            T[sgnvec] = []
            T[sgnvec].append(v)

    #print ' - Done'
    #print

    return T


def noisy_search(q, T, bound, query):
    
    # Find v in S within distance 'bound' from query
    # Input T is a hash table from S
    # If no such v, return None

    sgn_, index = sgnvar(query, bound, q)

    if len(index) == len(sgn_):
        return [0] * len(sgn_)

    #print ' - Current number of x in sgn\'(query) = ', len(index)

    v = check_collision(query, sgn_, T, index, 0, bound)

    #sys.stdout.write("\033[F")
    #sys.stdout.write("\033[K")

    return v


def check_collision(query, sgn_, T, index, index_num, bound):
    '''

	:param sgn_:	A ternary (-1,0,1) vector / See 'sgnvar' function
	:param index:	The index set where sgn_(i) = -1
	:param index_num:	A flag for recursive call which lies in (0, len(index))


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


def data_gen(A, ell):

	# Generate a set S = {As : Hw(s) <= ell}

    k = A.ncols()
    m = A.nrows()    
    A = transpose(A)
    S_ = set()
    S = set()
    Count = vector(ZZ, [0] * k)
    A = A.augment(Count)
    S_.add(tuple([0] * (m+1)))
    S.add(tuple([0] * m))
    #print 'Making a dataset S = { YA_2 * s : HW(s) < ell }...'

    for _ in range(ell):
        #print ' - Processing with hamming weight ', _, '... '

        TMP = set()
        TMP2 = set()
        for v in S_:
            for i in range(v[-1], k):
                for j in [0,1]:
                    tmp = vector(ZZ, v) + (-1)**j * A[i]
                    tmp[-1] = i + 1
                    TMP.add(tuple(tmp))
                    tmp = tmp[:m]
                    TMP2.add(tuple(tmp))
        S_ = S_.union(TMP)
        S = S.union(TMP2)

        #sys.stdout.write("\033[F")
        #sys.stdout.write("\033[K")

    #print ' - Done'
    #print

    return S


def solveLA(A, c, query, res, q):

    query = vector(ZZ, query)
    res = vector(ZZ, res)
    n = A.rank()
    error = query - res

    c -= error
    c = vector(GF(q, proof=False), c)

    A_ = A.matrix_from_rows([i for i in range(n)])
    if n != A_.rank():
        print 'A_ is not full-rank'
    else:
        c = c[:n]
        return balanced_lift(A_.inverse() * c)


def Mitm_on_LWE(A, c, u, bound, ell, check_unif = True):

	# Solve LWE on input (A, c) with error bound 'bound' by Mitm strategy
	# check_unif flag : whether try with input (A,u)

    q = A.base_ring().order()
    K = GF(q, proof = False)

    if ell % 2 == 0:
        half = ell // 2
    else:
        half = ell // 2 + 1

    S = data_gen(A, half) 
    T = generate_table(S, q)
    is_LWE = False

    print 'Table size = %d' % len(S)
    print

    print '** Mitm on (A_k, c_k) ** '
    print 

    count = 0
    print count
    for v in S:

        v = vector(K, v)
        query = c - v
        query = balanced_lift(query)
        
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K")

        count += 1
        print 'Number of noisy searches = %d' % count

        res = noisy_search(q, T, bound, query)        

        if res is not None:
            is_LWE = True
            if A.nrows() >= A.ncols():
                #print 'query            :', query
                #print 'collision in S   :', res
                #print
                s = solveLA(A, c, query, res, q)
                s = vector(ZZ, s)
            break
    
    if is_LWE == True:
        if A.nrows() < A.ncols():
            print ' - Input is LWE'
        else:
            print ' - Input is LWE with secret'
            print s

    else:
        print ' - Input is uniform'


    if check_unif is True:
        print
        print '** Mitm on (A_k, u_k) **'
        print

        is_LWE = False
        count = 0
        print count
        for v in S:

            v = vector(K, v)
            query = u - v
            query = balanced_lift(query)

            sys.stdout.write("\033[F")
            sys.stdout.write("\033[K")

            count += 1
            print 'Number of noisy searches = %d' % count

            res = noisy_search(q, T, bound, query)

            if res is not None:
                is_LWE = True
                break

        if is_LWE == True:
        	print ' - Input is LWE'

        else:
            print ' - Input is uniform'

def hybrid_mitm(n, q, h, beta, k, alpha = None, tau = None, ell = None, float_type="mpfr", check_unif = True):

    if tau is None:
        tau = 30

    K = GF(q, proof = False)

    # Lattice reduction stage
    A, c, u, s = gen_instance(n, q, h, m = tau * n)
    # Generate n-dim samples (A, c) (LWE) and (A, u) (Uniform).


    A_k, c_k, u_k, B = dim_error_tradeoff(A, c, u, beta, h, k, float_type=float_type)
    # Now have new k-dim samples (A_k, c_k) (LWE) and (A_k, u_k) (Uniform)
    # The error of (A_k, c_k) is bounded by B.

    # 'New k-dim samples have error bound = %5.2f, where q = %d' %(sqrt(2*pi) * 2 * bound, A.base_ring().order())

    if ell is None:
    	ell = h

    # MITM stage
    result = Mitm_on_LWE(A_k, c_k, u_k, B, ell, check_unif = check_unif) # Bound FIX

    