from estimator import *
import pprint

def hyb_estimate(f, n, alpha, q, bound, secret_distribution=True, **kwds):
    '''
    * Input 'bound' means,
    the lattice reduction will finds short vectors of (log) size < 'log(q) - bound'.
       
    * Our Strategy:
    take s = s1 + s2
    If HW(s1) <= l, HW(s2) <= k - l, then HW(s) <= k
            
    1. Store S = {(s1, b - As1) : s1} first
    sizeS = 2 ** (dim / 2)
    cost_post = tau * sizeS

    2. Compute As2 and check noisy collision.
    '''

    best = None
    best_k = None
    
    # too small a step size leads to an early abort, too large a step
    # size means stepping over target
    step_size = int(n / 8)

    if not SDis.is_ternary(secret_distribution):
        raise NotImplementedError("Only ternary secrets are currently supported.")

    a, b = SDis.bounds(secret_distribution)
    h = SDis.nonzero(secret_distribution, n)
    k = int(n / 2 - step_size)
    C = binomial(n, h)
    B = (1 + 1/sqrt(2 * pi)) * alpha * q * q / (sqrt(2) * (2 ** bound))

    while True:
        # <constants>
        # - cost_lat      : cost of lattice reduction (dual lattice attack)
        # - cost_post     : cost of postprocess (want to compute)
        # - cost_mem = N  : size of MITM table T
        # - mitm_num      : Table T consists of all c_1 with HW <= mitm_num
        cost_lat = 0
        cost_post = 0
        cost_mem = 0
        mitm_num = 0
        best_h1 = None

        k1 = k / 2
        k2 = k - k1             # Divide the matrix equally.      
        mem_bound = 2 ** 80     # Bound the memory capacity for MITM by 2**80.
        h1 = 5

        while h1 < k1: 
            h2 = 5

            while h2 <= k2 and h2 <= h1 + 4:
                current = f(n - k, alpha, q, bound, k1, k2, h1, h2, secret_distribution=secret_distribution, **kwds)
                cost_lat = current["red"] 
                repeat = current["repeat"]               
                                 
                cost_tableConstruct = 0
                for i in range(1, h1 + 1):
                    cost_tableConstruct += binomial(k1, i) * (b-a) ** i * ( i * repeat + repeat)

                cost_query = 0
                for j in range(1, h2 + 1):
                    cost_query += binomial(k2, j) * (b-a)**j * ( 2 ** (repeat * 4 * B / q) )                    

                probability = 0
                for i in range(0, h1 + 1):
                    for j in range(0, h2 + 1):
                        probability += (binomial(n - k1 - k2, h - i - j) * binomial(k1, i) * binomial(k2, j)) / C
                    
                cost_post = cost_tableConstruct + cost_query

                total_cost = (cost_post + cost_lat) / probability

                current["rop"] = total_cost
                current["post"] = cost_post
                current["prob_inv"] = int(1/probability)
                current["k"] = k
                current["h1"] = h1
                current["h2"] = h2
                current = current.reorder(["rop"])

                logging.getLogger("guess").debug(current)

                if int(cost_query) > int(cost_lat):
                    best_h1 = current
                    best_h1["rop"] = oo
                    #print 'mem = %5.1f, cost_tableConstruct = %5.1f, cost_query = %5.1f, cost_lat = %5.1f' % (log(mem,2), log(cost_tableConstruct,2),log(cost_query,2), log(cost_lat,2))
                    break
                
                if best_h1 is None or current["rop"] < best_h1["rop"]:
                    best_h1 = current
                h2 += 1

            if best_k is None or best_h1["rop"] < best_k["rop"]:
                best_k = best_h1

            for i in range(0, h1 + 1):
                cost_mem += best_k["repeat"] * binomial(k1, i) * (b-a) ** i # bits       
                best_k["mem"] = cost_mem 

            if cost_mem + best_k["repeat"] * binomial(k1, h1 + 1) * (b-a) ** h1 > mem_bound:                
                break

            #print 'k = %4d, rop = %5.1f, h1 = %4d, h2 = %4d' % (best_h1["k"], log(best_h1["rop"], 2), best_h1["h1"], best_h1["h2"])
            h1 += 1
                                
        if best is None or best_k["rop"] < best["rop"]:
            best = best_k
            k += step_size
        else:
            # we go back
            k = best["k"] - step_size
            k += step_size / 2
            if k <= 0:
                k = step_size / 2
            # and half the step size
            step_size = step_size / 2

        if step_size == 0:
            break
    
    #print("Optimal k = " + str(k) + ", k1 = " + str(k1) +", k2 = " + str(k2) + ", h1 = " + str(h1_best) + ", h2 = " + str(h2_h1_best))
    return best

def dual_scale_hyb(n, alpha, q, bound, k1, k2, h1, h2, secret_distribution,
                m=oo, reduction_cost_model=reduction_default_cost, c=None, use_lll=True):
    
    T = binomial(k1, h1) * (2 ** h1)
    Q = binomial(k2, h2) * (2 ** h2)
    B = (1 + 1/sqrt(2 * pi)) * alpha * q * q / (sqrt(2) * (2 ** bound))

    # stddev of the error
    e = RR(stddevf(alpha * q))

    if SDis.is_ternary(secret_distribution):
        a, b = SDis.bounds(secret_distribution)
        h = SDis.nonzero(secret_distribution, n)
        e_ = RR(1)
        if c is None:
            c = RR(e * sqrt(2 * n - n) / sqrt(h))
    else:
        if not SDis.is_small(secret_distribution):
            m = m - n
        c = RR(1)
        e_ = e
        h = n

    best = dual(n=n, alpha=alpha, q=q, m=m, reduction_cost_model=reduction_cost_model)
    delta_0 = best["delta_0"]

    if use_lll:
        scale = 2
    else:
        scale = 1

    while True:
        m_optimal = lattice_reduction_opt_m(n, q / c, delta_0)
        m_ = ZZ(min(m_optimal, m + n))

        v = scale * delta_0 ** m_ * (q / c) ** (n / m_)

        if log(v, 2) > log(q, 2) - bound:
            delta_0 = delta_0 - RR(0.00005)

        else:
            ret = lattice_reduction_cost(reduction_cost_model, delta_0, m_, B=log(q, 2))

            repeat = log(T * Q, 2) / (1 - 4 * B / q)

            if use_lll:
                lll = BKZ.LLL(m_, log(q, 2))
            else:
                lll = None
            ret = ret.repeat(times=repeat, lll=lll)

            ret[u"m"] = m_
            ret[u"repeat"] = repeat
            ret[u"d"] = m_
            ret[u"c"] = c

            ret = ret.reorder(["rop", "m"])
            logging.getLogger("repeat").debug(ret)
            best = ret
            break

    return best

def MITM_estimator(n, alpha, q, h = 64, start_bound = 10, Max_bound = 13, step_size = 1, reduction_cost_model = reduction_default_cost):
    bound = float(start_bound)
    mitm_hyb = partial(hyb_estimate, dual_scale_hyb)
    best = None
    bound_best = bound
    print 'Chosen Parameters : '
    print '     n = %5d, log(q) = %5.1f, stddev = %5.2f, HW(s) = %4d' % (n, log(q,2), RR(stddevf(alpha * q)), h)
    print
    Level = 0
    step = float(step_size)

    print 'Start Estimation . . .'
    print 

    while Level <= 2:

        res = mitm_hyb(n, alpha, q, bound, secret_distribution=((-1, 1), h), reduction_cost_model=reduction_cost_model)
        
        print "Optimizing with beta = %4d . . ." % res["beta"]

        if best is None or res["rop"] < best["rop"]:            
            best = res
            bound_best = bound
            while bound + step > Max_bound and Level <= 2:
                step = float(step/2)
                Level += 1
            bound += step

        else:
            step = float(step/2)
            Level += 1
            if Level <=2:
                bound -= step
    print
    print '== Bit-security : %5.1f with optimal parameters' % log(best["rop"], 2)
    print '     k = %5d, h1 = %2d, h2 = %2d, beta = %4d, mem = %5.1f' % (best["k"], best["h1"], best["h2"], best["beta"], log(best["mem"],2))
    print '             (For simplicity, we set k1 = k2 = k/2)'
    print
    print "== Details"
    pprint.pprint(best)
