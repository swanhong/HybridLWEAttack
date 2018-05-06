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
        pow += 2**i * v[i]
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


def dim_error_tradeoff(A, c, u, beta, h, s, k, num_sample = None, float_type="double", use_lll = True):
    """

    :param A:    LWE matrix
    :param c:    LWE vector
    :param u:	 Uniform vector
    :param beta: BKW block size
    :param h: 	 Hamming weight of secret
    :param k:    LWE dim after tradeoff
    :param num_sample:  number of new samples to generate
    :param use_lll: 	If True, run BKZ only once and then run LLL
    					If False, run BKZ iteratively

    * secret vector s is used to see the error term of new LWE(-like) samples.

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
    L = LLL.reduction(L)
    M = GSO.Mat(L, float_type=float_type)
    bkz = BKZ2(M)
    if num_sample > 1:
        param = BKZ.Param(block_size=beta,
                      strategies=BKZ.DEFAULT_STRATEGY,
                      auto_abort=True,
                      max_loops=16,
                      flags=BKZ.VERBOSE|BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
    else:
        param = BKZ.Param(block_size=beta,
                      strategies=BKZ.DEFAULT_STRATEGY,
                      auto_abort=True,
                      max_loops=16,
                      flags=BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
    bkz(param)

    if num_sample > 1:
        print
        print '** Generating short vectors by LLL **'
        print

    H = copy(L)
    
    A2 = A.matrix_from_columns([n - k + i for i in range(k)])
    E = matrix(ZZ, 1, k)
    f_c = []
    f_u = []
    Y = set()
    y_i = vector(ZZ, tuple(L[0]))
    E[0], ft, ft2 = apply_short1(y_i, A2, c, u, scale=scale)
    f_c.append(ft)
    f_u.append(ft2)

    y = L[0].norm()
    y_ = y/sqrt(L.ncols)
    y_r = 3.2*sqrt(L.ncols - A2.ncols())*y_/scale
    y_l = sqrt(h)*y_

    count = 0

    if num_sample > 1:
        print 'Expected log(E[error]) = %5.1f, BKZ vector log(||y_1||) = %5.1f' % (log(sqrt(y_r**2 + y_l**2), 2), log(y, 2))
        print
        print '# of LLL : %d, Current |Y| : %d' % (count, len(Y))
    
        if num_sample is None:
            num_sample = 30
    
        if use_lll is True:
            while len(Y) < num_sample:
                count += 1
                M = GSO.Mat(L, float_type=float_type)
                bkz = BKZ2(M)
                bkz.randomize_block(0, L.nrows, density=3)
                LLL.reduction(L)
                y_i = vector(ZZ, tuple(L[0]))
                if L[0].norm() > H[0].norm():
                    L = copy(H)

                sys.stdout.write("\033[F")
                sys.stdout.write("\033[K")

                tmp = len(Y)
                Y.add(y_i.norm())
                if tmp != len(Y):
                    Etmp, ftmp, ftmp2 = apply_short1(y_i, A2, c, u, scale=scale)
                    E = E.stack(Etmp)
                    f_c.append(ftmp)
                    f_u.append(ftmp2)
                print '# of LLL : %d, Current |Y| : %d' % (count, len(Y))

        else:
            while len(Y) < num_sample:
                count += 1
                M = GSO.Mat(L, float_type=float_type)
                bkz = BKZ2(M)
                bkz.randomize_block(0, L.nrows, density=3)
                bkz(param)
                y_i = vector(ZZ, tuple(L[0]))
                if L[0].norm() > H[0].norm():
                    L = copy(H)

                sys.stdout.write("\033[F")
                sys.stdout.write("\033[K")

                tmp = len(Y)
                Y.add(y_i.norm())
                if tmp != len(Y):
                    Etmp, ftmp, ftmp2 = apply_short1(y_i, A2, c, u, scale=scale)
                    E = E.stack(Etmp)
                    f_c.append(ftmp)
                    f_u.append(ftmp2)

                print '# of BKZ : %d, Current |Y| : %d' % (count, len(Y))

    E = E.change_ring(K)    
    f_c = vector(K, f_c)
    f_u = vector(K, f_u)
    s = s[-k:]
    e = f_c - E * s
    e = balanced_lift(e)
    
    # print e

    bound = sqrt(y_r**2 + y_l**2)

    if num_sample > 1:
        print
        print 'log(E[error]) : (%5.1f, %5.1f), log(E[||y_i||]) = (%5.1f, %5.1f) ' \
            % (log_mean(e), log_stddev(e), log_mean(Y), log_stddev(Y))

        print

    return E, f_c, f_u, bound


def generate_table(S, q):   

	# Generate (binary) hash table T induced from locality hash function 'sgn'

    print 'Generating table ... '

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

    print ' - Done'
    print

    return T


def noisy_search(q, T, bound, query):
    
    # Find v in S within distance 'bound' from query
    # Input T is a hash table from S
    # If no such v, return None

    sgn_, index = sgnvar(query, bound, q)

    if len(index) == len(sgn_):
        return [0] * len(sgn_)

    print ' - Current number of x in sgn\'(query) = ', len(index)

    v = check_collision(query, sgn_, T, index, 0, bound)

    sys.stdout.write("\033[F")
    sys.stdout.write("\033[K")

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
    print 'Making a dataset S = { YA_2 * s : HW(s) < ell }...'

    for _ in range(ell):
        print ' - Processing with hamming weight ', _, '... '

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

        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K")

    print ' - Done'
    print

    return S


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
        return balanced_lift(A_.inverse() * c)


def Mitm_on_LWE(A, c, u, bound, ell, check_unif = False):

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

    print '** Mitm with (YA_2, Yc) ** '

    for v in S:
        v = vector(K, v)
        query = c - v
        query = balanced_lift(query)
        
        res = noisy_search(q, T, bound, query)

        if res is not None:
            is_LWE = True
            if A.nrows() >= A.ncols():
                print 'query            :', query
                print 'collision in S   :', res
                print
                s = solveLA(A, c, query, res, q)
                s = vector(ZZ, s)
            break
    
    if is_LWE == True:
        if A.nrows() < A.ncols():
            print ' - Input is LWE'
        else:
            print 'Input is from LWE samples with secret s:'
            print s

    else:
        print ' - Input is uniform'


    if check_unif is True:
        print
        print '** Mitm with (YA_2, Yu) **'

        is_LWE = False

        for v in S:
            v = vector(K, v)
            query = u - v
            query = balanced_lift(query)

            res = noisy_search(q, T, bound, query)

            if res is not None:
                is_LWE = True
                break

        if is_LWE == True:
        	print ' - Input is LWE'

        else:
            print ' - Input is uniform'
   

def hybrid_mitm(A, c, u, beta, h, s, k, num_sample = None, ell = None, float_type="double", check_unif = False):
    
    # if flag = False, then don't check for uniform vector u

    # Lattice reduction stage
    E, f_c, f_u, bound = dim_error_tradeoff(A, c, u, beta, h, s, k, num_sample = num_sample, float_type=float_type)

    print 'Bound setting = %5.2f, q = %d' %(sqrt(2*pi) * 2 * bound, A.base_ring().order())
    print

    if ell is None:
    	ell = h

    # MITM stage
    Mitm_on_LWE(E, f_c, f_u, sqrt(2*pi) * 2 * bound, ell, check_unif = check_unif)

def hybrid_mitm_(n, q, beta, h, k, alpha = None, tau = None, ell = None, float_type="double", check_unif = False):

    if tau is None:
        tau = 30

    K = GF(q, proof = False)

    # Lattice reduction stage
    A, c, u, s = gen_instance(n, q, h)
    E, f_c, f_u, bound = dim_error_tradeoff(A, c, u, beta, h, s, k, num_sample = 1, float_type=float_type)
    f_c = f_c.list()
    f_u = f_u.list()

    print s

    for i in range(1, tau):
        print 'Generating ...' 
        A2, c2, u2, s2 = gen_instance(n, q, h, s = s)
        E2, f_c2, f_u2, bound = dim_error_tradeoff(A2, c2, u2, beta, h, s, k, num_sample = 1, float_type=float_type)
        E = E.stack(E2)
        f_c.append(f_c2[0])
        f_u.append(f_u2[0])

    f_c = vector(K, f_c)
    f_u = vector(K, f_u)

    print 'Bound setting = %5.2f, q = %d' %(sqrt(2*pi) * 2 * bound, A.base_ring().order())
    print

    if ell is None:
    	ell = h

    # MITM stage
    Mitm_on_LWE(E, f_c, f_u, sqrt(2*pi) * 2 * bound, ell, check_unif = check_unif)


    


def silke(A, c, beta, h, max_loops = 16, num_sample = None, float_type="double", use_lll = True):
    """

    :param A:    LWE matrix
    :param c:    LWE vector
    :param beta: BKW block size
    :param h: 	 Hamming weight of secret
    :param k:    LWE dim after tradeoff
    :param num_sample:  number of new samples to generate

    """

    from fpylll import BKZ, IntegerMatrix, LLL, GSO
    from fpylll.algorithms.bkz2 import BKZReduction as BKZ2

    n = A.ncols()
    q = A.base_ring().order()
    K = GF(q, proof=False)

    scale = round(8 * sqrt(n) / sqrt(2*pi*h))
    scale = ZZ(scale)
    
    L = dual_instance1(A, scale=scale)
    L = IntegerMatrix.from_matrix(L)

    L1 = copy(L)
    L1 = LLL.reduction(L1, flags=LLL.VERBOSE)
    M = GSO.Mat(L1, float_type=float_type)
    bkz = BKZ2(M)
    param = BKZ.Param(block_size=beta,
    				  strategies=BKZ.DEFAULT_STRATEGY,
    				  auto_abort=True,
    				  max_loops=max_loops,
    				  flags=BKZ.VERBOSE|BKZ.AUTO_ABORT|BKZ.MAX_LOOPS)
    bkz(param)
    
    H = copy(L1)
    
    f_c = []
    V = set()
    y_i = vector(ZZ, tuple(L[0]))
    ft = apply_short1(y_i, A, c, scale=scale)[1]
    f_c.append(ft)

    v = L1[0].norm()
    v_ = v/sqrt(L1.ncols)
    v_r = 3.2*sqrt(L1.ncols - A.ncols())*v_/scale
    v_l = sqrt(h)*v_

    print 'Expected log(E[error]) = %5.1f, BKZ vector log(||y_1||) = %5.1f' % (log(sqrt(v_r**2 + v_l**2), 2), log(v, 2))

    print 'Current |Y| : ', len(V)
    	
    if num_sample is None:
    	num_sample = 30

    count = 0
    if use_lll is True:
        while len(V) < num_sample:
            count += 1
            M = GSO.Mat(L, float_type=float_type)
            bkz = BKZ2(M)
            bkz.randomize_block(0, L.nrows, density=5)
            LLL.reduction(L)
            y_i = vector(ZZ, tuple(L[0]))
            if L[0].norm() > H[0].norm():
                L = copy(H)

            sys.stdout.write("\033[F")
            sys.stdout.write("\033[K")

            flag = len(V)
            V.add(y_i.norm())
            if flag != len(V):
                ftmp = apply_short1(y_i, A, c, scale=scale)[1]
                f_c.append(ftmp)
            print '# of LLL : %d, Current |Y| : %d' % (count, len(V))

    else:
        while len(V) < num_sample:
            count += 1
            M = GSO.Mat(L, float_type=float_type)
            bkz = BKZ2(M)
            bkz.randomize_block(0, L.nrows, density=3)
            bkz(param)
            y_i = vector(ZZ, tuple(L[0]))
            if L[0].norm() > H[0].norm():
                L = copy(H)

            sys.stdout.write("\033[F")
            sys.stdout.write("\033[K")

            flag = len(V)
            V.add(y_i.norm())
            if flag != len(V):
                ftmp = apply_short1(y_i, A, c, scale=scale)[1]
                f_c.append(ftmp)
            print '# of BKZ : %d, Current |Y| : %d' % (count, len(V))


    print	
    print 'log(E[error]) : (%5.1f, %5.1f), log(E[||y_i||]) = (%5.1f, %5.1f) ' \
    	% (log_mean(f_c), log_stddev(f_c), log_mean(V), log_stddev(V))

    print

    return f_c


    