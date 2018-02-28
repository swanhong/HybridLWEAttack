# HybridLWEAttack

Our estimator uses functions in "estimator.py" made by Albrecht [Alb17].
Most of our functions are simply modifed from "estimator.py", just fit our algorithm.

### How to Run

We recommand to use our **MITM_estimator** function in "estimator_Hybrid.py". This function estimates the cost of ring operations and memory.

EXAMPLE :\
    * sage: from estimator_hybrid import *\
    * sage: n = 8192\
    * sage: q = next_prime(2^125)\
    * sage: alpha = 8 / q\
    * sage: MITM_estimator(n, alpha, q, 11, 1, 13)\

        Chosen Constants :
            n = 8192, q = 42535295865117307932921825928971026459

        == Start Estimation ==

        Current bound 11.0
        Optimal k = 5368, l1 = 3020, l2 = 2350, h1 = 6, h2 = 7
        Current bound 12.0
        Optimal k = 5343, l1 = 2671, l2 = 2672, h1 = 6, h2 = 7
        Current bound 11.5
        Optimal k = 5372, l1 = 2687, l2 = 2687, h1 = 6, h2 = 7
        Current bound 11.25
        Optimal k = 5311, l1 = 2986, l2 = 2325, h1 = 6, h2 = 7
        
        Optimal bound = 11.25
        rop:  2^132.6
        m:   2^12.6
        red:   2^87.1
        delta_0: 1.006505
        beta:      189
        repeat:      160
        d:   2^12.6
        c:   21.197
        post:  2^132.6
        k:   2^12.4
        postprocess:        1
        mem:   2^73.2
