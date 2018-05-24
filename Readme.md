# HybridLWEAttack

Our "estimator_hybrid.py" estimates the bit-security of the hybrid MITM attack, given parameter algorithm.
We uses functions in "estimator.py" made by Albrecht et al [[APS15]](https://eprint.iacr.org/2015/046)
and the high-level structure of our estimator "estimator_hybrid.py" is from "estimator.py".

### How to Run

We recommend to use our **MITM_estimator** function in "estimator_hybrid.py". This function estimates the cost of ring operations, given a bounded memory capacity 2^80.

EXAMPLE :\
    * sage: from estimator_hybrid import *\
    * sage: n = 8192\
    * sage: q = next_prime(2^125)\
    * sage: alpha = 8 / q\
    * sage: MITM_estimator(n, alpha, q)
    
    
    

    Chosen Parameters :
             n =  8192, log(q) = 125.0, stddev =  3.19, HW(s) =   64
     
    Start Estimation . . .

    Optimizing with beta =  240 . . .
    Optimizing with beta =  247 . . .
    Optimizing with beta =  243 . . .
    Optimizing with beta =  244 . . .

    == Bit-security : 129.9 with optimal parameters
         k =  4861, h1 =  7, h2 =  9, beta =  240, mem =  80.8
                 (For simplicity, we set k1 = k2 = k/2)

    == Details
    rop:  2^129.9
    m:   2^12.8
    red:  2^102.3
    delta_0: 1.005603
    beta:      240
    repeat:  170.419
    d:   2^12.8
    c:   23.025
    post:   2^97.0
    prob_inv:   2^27.6
    k:   2^12.2
    h1:        7
    h2:        9
    mem:   2^80.8