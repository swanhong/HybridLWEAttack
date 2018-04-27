from estimator import *


def drop_and_solve_hyb(f, n, alpha, q, bound, secret_distribution=True, success_probability=0.99,
                       postprocess=True, decision=True, **kwds):
    '''
    <Additional Inputs>
    - bound = Maximum memory. memory <= 2 ** bound.
    '''
    n, alpha, q, success_probability = Param.preprocess(n, alpha, q, success_probability)

    RR = parent(alpha)

    best = None

    if not decision and postprocess:
        raise ValueError("Postprocessing is only defined for the dual attack which solves the decision version.")

    # too small a step size leads to an early abort, too large a step
    # size means stepping over target
    step_size = int(n / 16)

    if not SDis.is_ternary(secret_distribution):
        raise NotImplementedError("Only ternary secrets are currently supported.")

    a, b = SDis.bounds(secret_distribution)
    h = SDis.nonzero(secret_distribution, n)
    k = 0
    C = binomial(n, h)

    while True:
        probability = RR(success_probability_drop(n, h, k))
        # increase precision until the probability is meaningful
        while success_probability ** probability == 1:
            success_probability = RealField(64 + success_probability.prec())(success_probability)

        current = f(n - k, alpha, q, bound,
                    success_probability=success_probability ** probability if decision else success_probability,
                    secret_distribution=secret_distribution, **kwds)

        # <constants>
        # - cost_lat      : cost of lattice reduction (dual lattice attack)
        # - cost_post     : cost of postprocess (want to compute)
        # - cost_mem = N  : size of MITM table T
        # - mitm_num      : Table T consists of all c_1 with HW <= mitm_num
        cost_lat = current.values()[0]
        cost_post = 0
        cost_mem = 0
        mitm_num = 0

        if postprocess:
            repeat = current["repeat"]
            dim = current["d"]
            '''
            repeat <- tau
            A2 : m by k matrix
            dim <- m
            length of y <= q / (2 ** bound)
            '''
            ### <STEP1> Constructing the table of T of size (cost_mem)

            # - compute YA_2
            # # of inner product = repeat * k
            # # of rop in each inner product = (dim) mult + (dim - 1) sum ~= 2 * dim
            cost_post = (repeat * k) * (2 * dim)

            # determine HW = mitm_num to store
            for i in range(1, k):
                # number of vectors of HW = i
                cost_mem_i = repeat * binomial(k / 2, i) * (b-a)**i

                # if current memory >= bound, then break
                # < y, e > <= B
                B = 4 * alpha * q * q / (2 ** bound)
                if cost_mem + cost_mem_i >= 2 ** (repeat / 2 * (1 - B / q)):
                    mitm_num = i - 1
                    postprocess = 2 * mitm_num
                    break
                cost_mem += cost_mem_i

                # compute YA_{2,l}s_l and YA_{2,r}s_r for HW = i
                # column vectors of YA_2 are already stored
                # To compute YA_{2,1}s_1, pick i columns in YA_{2,1} and do i vector sum for all s_1's.
                # Do same for YA_{2,2}s_2.
                cost_post += binomial(k / 2, i) * (b-a)**i * ( i * repeat ) * 2

            # probability = Prob[ nonzero elements of s(k slots) splits to 2 blocks,
            #                     so that each block contains <=l nonzero elements ]
            if mitm_num != 0:
                probability = 0
                for i in range(0, mitm_num+ 1):
                    for j in range(0, mitm_num + 1):
                        if k > i + j:
                            probability += succ_prob_fcn(n, h, k, i, j, C)

            ### <STEP2> Sorting the table T
            # we can use Heap Sort algorithm, it needs only 1 memory
            if cost_mem > 0:
                cost_post += cost_mem * log(cost_mem, 2)

            ### <STEP3> Searching B-noisy collision
            # overall binary search time
            if cost_mem != 0:
                cost_post += cost_mem * log(cost_mem, 2) * (2 ** (repeat * 4 * B / q))

        current["rop"] = (cost_lat + cost_post) / probability
        current["post"] = cost_post
        # current = current.repeat(1 / probability, select={"m": False})
        current["k"] = k
        current["postprocess"] = postprocess
        current["mem"] = cost_mem
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

    print("Exact k", str(k))
    return best

def drop_and_solve_hybOpt(f, n, alpha, q, bound, secret_distribution=True, success_probability=0.99,
                    postprocess=True, decision=True, **kwds):
    '''
    <Additional Inputs>
    - bound = Maximum memory. memory <= 2 ** bound.
    '''
    n, alpha, q, success_probability = Param.preprocess(n, alpha, q, success_probability)

    RR = parent(alpha)

    best = None

    if not decision and postprocess:
        raise ValueError("Postprocessing is only defined for the dual attack which solves the decision version.")

    # too small a step size leads to an early abort, too large a step
    # size means stepping over target
    step_size = int(n / 32)

    if not SDis.is_ternary(secret_distribution):
        raise NotImplementedError("Only ternary secrets are currently supported.")

    a, b = SDis.bounds(secret_distribution)
    h = SDis.nonzero(secret_distribution, n)
    k = step_size
    C = binomial(n, h)

    while True:
        # print("k = " + str(k))
        # print("k = " + str(k))
        probability = RR(success_probability_drop(n, h, k))
        # increase precision until the probability is meaningful
        while success_probability ** probability == 1:
            success_probability = RealField(64 + success_probability.prec())(success_probability)

        current = f(n - k, alpha, q, bound,
                    success_probability=success_probability ** probability if decision else success_probability,
                    secret_distribution=secret_distribution, **kwds)

        # <constants>
        # - cost_lat      : cost of lattice reduction (dual lattice attack)
        # - cost_post     : cost of postprocess (want to compute)
        # - cost_mem = N  : size of MITM table T
        # - mitm_num      : Table T consists of all c_1 with HW <= mitm_num
        cost_lat = current.values()[0]
        cost_post = 0
        cost_mem = 0
        mitm_num = 0
        total_post_cost = 0

        if postprocess:

            '''
            <Strategy>
            take s = s1 + s2
            If HW(s1) <= l, HW(s2) <= k - l, then HW(s) <= k
            
            1. Store S = {(s1, b - As1) : s1} first
            sizeS = 2 ** (dim / 2)
            cost_post = tau * sizeS

            2. Compute As2 and check noisy collision.
            '''

            repeat = current["repeat"]
            dim = current["d"]
            
            '''
            < STEP1 > Constucting hash table H.
            repeat <- tau
            dim <- lattice dimension

            store all vector with Hamming weight <= mitm_num
            '''

            l1 = k / 2
            l1_best = 0
            step_size_l1 = k / 16
            
            cost_best_l1 = 0

            while 1: # do for l1
                # print("l1 = " + str(l1))
                l2 = k - l1
                
                mem_bound = 2 ** 100
                B = 2 * alpha * q * q / (2 ** bound)

                h1 = 0
                h1_best = 0
                h2_best = 0
                cost_mem = 0

                cost_tableConstruct = 0
                
                cost_best_h1 = cost_best_l1
                
                while h1 < l1: # do for h1
                    h1 += 1
                    # number of coord = repeat
                    # candidates of s with HW i = (k choose i) * (b-a)^i 
                    candidates_i = binomial(l1, h1) * (b-a) ** h1
                    cost_mem_i = repeat * candidates_i # bits
                    
                    if cost_mem + cost_mem_i > mem_bound:
                        h1 -= 1
                        break

                    cost_mem += cost_mem_i
                    
                    # Computation of As1 for HW = i
                    # Computing As1 : for each s1, pick i columns of A and do i number of vector sum.
                    # c - As1 : tau times
                    cost_tableConstruct += binomial(l1, h1) * (b-a) ** h1 * ( h1 * repeat + repeat)
                    h2 = 0

                    cost_query = 0

                    cost_best_h2 = cost_best_h1

                    while cost_query <= cost_tableConstruct and h2 <= l2:
                        h2 += 1
                        cost_query += binomial(l2, h2) * (b-a)**h2 * ( h2 * repeat )
                        cost_query += binomial(l2, h2) * (b-a)**h2 * ( 2 ** (repeat * 4 * B / q) )
                    
                    # print("(h1, h2) = (" + str(h1) + ", " + str(h2) + ")")

                        probability = 0
                        for i in range(h1 + 1):
                            for j in range(h2 + 1):
                                if k > i + j:
                                    probability += (binomial(n - l1 - l2, h - i - j) * binomial(l1, i) * binomial(l2, j)) / C
                        
                        total_post_cost = cost_tableConstruct + cost_query
                        total_cost = (total_post_cost + cost_lat) / probability

                        if cost_best_h2 is 0 or total_cost < cost_best_h2:
                            h2_best = h2
                            cost_best_h2 = total_cost

                    # print("cost_mem", int(log(cost_mem, 2)))
                    # print("mem_bound", int(log(mem_bound, 2)))
                    # print("h1, h2", h1, h2)    
                    # print("total cost", int(log(total_cost, 2)))
                    # print("table construct", int(log(cost_tableConstruct, 2)))
                    # print("query", int(log(cost_query, 2)))
                    # print("prob", int(log(probability, 2)))
                    if cost_best_h1 is 0 or cost_best_h2 < cost_best_h1:
                        h1_best = h1
                        cost_best_h1 = cost_best_h2

                h1 = h1_best
                h2 = h2_best
                # print("best h1, h2", h1, h2)
                # print(" ======= ======= ======= ======= ======= ======= ")
                # print("")

                if cost_best_l1 is 0 or cost_best_h1 < cost_best_l1:
                    l1_best = l1
                    cost_best_l1 = cost_best_h1
                    l1 += step_size_l1
                    # print("l1", l1, "jumped", step_size_l1)
                    if l1 >= k:
                        while l1 >= k and step_size_l1 > 1:
                            step_size_l1 /= 2     
                            l1 -= step_size_l1

                else:
                    step_size_l1 /= 2
                    l1 -= step_size_l1
                    # print("l1", l1, "back", step_size_l1)

                if step_size_l1 <= 1:
                    break
            
            l1 = l1_best
            l2 = k - l1
            
        cost_post = total_post_cost

        current["rop"] = (cost_lat + cost_post) / probability
        current["post"] = cost_post
        # current = current.repeat(1 / probability, select={"m": False})
        # current = current.repeat(1, select={"m": False})
        current["k"] = k
        current["postprocess"] = h1 + h2
        current["mem"] = cost_mem
        # current["l2"] = query_HW
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
    
    print("Optimal k = " + str(k) + ", l1 = " + str(l1) +", l2 = " + str(l2) + ", h1 = " + str(h1) + ", h2 = " + str(h2))
    return best

def dual_scale_hyb(n, alpha, q, bound, secret_distribution,
                m=oo, success_probability=0.99,
                reduction_cost_model=reduction_default_cost,
                c=None, use_lll=True):
    n, alpha, q, success_probability = Param.preprocess(n, alpha, q, success_probability)

    RR = parent(alpha)

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

    best = dual(n=n, alpha=alpha, q=q, m=m,
                reduction_cost_model=reduction_cost_model)
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

            repeat = 160
            # B = 4 * alpha * q * (q / 2 ** bound)
            # B = 4 * alpha * q * (q / 2 ** bound)
            # C = ret["red"] # BKZ coeff
            # repeat = 1 / (1 + 4 * B / q) * 2 * (log(C, 2) - log((1 - 4 * B / q) / 2, 2))

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
    duald = partial(drop_and_solve_hybOpt, dual_scale_hyb)
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
