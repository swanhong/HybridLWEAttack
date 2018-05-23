from estimator import *

def drop_and_solve_hyb(f, n, alpha, q, bound, secret_distribution=True, **kwds):
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

        l1 = k / 2
        #l1_best = 0
        #step_size_l1 = k / 16
            
        cost_best_l1 = 0


        l2 = k - l1                
        mem_bound = 2 ** 80 # We bound the memory capacity for MITM by 2**80.
        h1 = 6
        h1_best = 0
        h2_best = 0
        h2_h1_best = 0
        mem_h1_best = 0
                
        cost_best_h1 = cost_best_l1

        sum_of_candidates = 0

        while h1 < l1:
            h1 += 1
            repeat = 0                                
            cost_best_h2 = cost_best_h1
            mem_h2_best = 0

            flag = True
            h2 = 8          

            while h2 <= l2 and h2 <= h1 + 5:
                h2 += 1
                current = f(n - k, alpha, q, bound, l1, l2, h1, h2, secret_distribution=secret_distribution, **kwds)
                cost_lat = current.values()[0] 
                repeat = current["repeat"]                     

                cost_tableConstruct = 0
                for i in range(1, h1 + 1):
                    cost_tableConstruct += binomial(l1, i) * (b-a) ** i * ( i * repeat + repeat)

                cost_query = 0
                for j in range(1, h2 + 1):
                    cost_query += binomial(l2, j) * (b-a)**j * ( 2 ** (repeat * 4 * B / q) )                    

                probability = 0
                for i in range(0, h1 + 1):
                    for j in range(0, h2 + 1):
                        probability += (binomial(n - l1 - l2, h - i - j) * binomial(l1, i) * binomial(l2, j)) / C
                        
                cost_post = cost_tableConstruct + cost_query

                total_cost = (cost_post + cost_lat) / probability

                mem = 0
                for i in range(0, h1 + 1):
                    mem += repeat * binomial(l1, i) * (b-a) ** i # bits
                    if mem > mem_bound:
                        mem -= repeat * binomial(l1, i) * (b-a) ** i
                    
                if cost_best_h2 is 0 or total_cost < cost_best_h2:
                    h2_best = h2
                    cost_best_h2 = total_cost
                    mem_h2_best = mem
                    
                if cost_query > cost_lat:
                    #print 'mem = %5.1f, cost_tableConstruct = %5.1f, cost_query = %5.1f, cost_lat = %5.1f' % (log(mem,2), log(cost_tableConstruct,2),log(cost_query,2), log(cost_lat,2))
                    h2 -= 1
                    break

                #print("l1 = " + str(l1) + ", l2 = " + str(l2) + ", h1 = " + str(h1) + ", h2 = " + str(h2))
                #print 'Total cost = %5.1f' % log(total_cost, 2)   

            if cost_best_h1 is 0 or cost_best_h2 < cost_best_h1:
                h1_best = h1
                h2_h1_best = h2_best
                cost_best_h1 = cost_best_h2
                mem_h1_best = mem_h2_best

            for i in range(0, h1 + 1):
                mem += repeat * binomial(l1, i) * (b-a) ** i # bits
                if mem > mem_bound:
                    mem -= repeat * binomial(l1, i) * (b-a) ** i     
                    flag = False                       

            if flag is False:
                break
            
        current["rop"] = cost_best_h1
        current["post"] = cost_post
        current["k"] = k
        current["postprocess"] = h1_best + h2_h1_best
        current["mem"] = mem_h1_best
        current = current.reorder(["rop"])

        logging.getLogger("guess").debug(current)

        key = list(current)[0]
        if best is None:
            best = current
            k += step_size
            continue

        if current[key] < best[key]:
            best = current
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
    
    #print("Optimal k = " + str(k) + ", l1 = " + str(l1) +", l2 = " + str(l2) + ", h1 = " + str(h1_best) + ", h2 = " + str(h2_h1_best))
    return best

def dual_scale_hyb(n, alpha, q, bound, l1, l2, h1, h2, secret_distribution,
                m=oo, reduction_cost_model=reduction_default_cost, c=None, use_lll=True):
    
    T = binomial(l1, h1) * (2 ** h1)
    Q = binomial(l2, h2) * (2 ** h2)
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

def MITM_estimator(n, alpha, q, start_bound, step_size, Max_bound):
    bound = float(start_bound)
    duald = partial(drop_and_solve_hyb, dual_scale_hyb)
    best = None
    print("Chosen Constants : ")
    print("     n = " + str(n) + ", q = " + str(q))
    print("")
    cost_best = 0
    Level = 0

    print("== Start Estimation ==")
    print("")
    
    while Level <= 2:
        print("Current bound = " + str(bound))
        res = duald(n, alpha, q, bound, secret_distribution=((-1, 1), 64))
        
        if best is None:
            best = res

        rop_key = list(res)[0]
        rop = res[rop_key]

        if cost_best is 0 or rop < cost_best:
            cost_best = rop
            bound += step_size
            best = res
            while bound >= Max_bound and step_size >= 1/4:
                step_size = float(step_size/2)
                bound -= step_size
        else:
            step_size = float(step_size/2)
            Level += 1
            if Level <=2:
                bound -= step_size

    print("")
    print("Optimal bound = " + str(bound))
    pprint(best)
