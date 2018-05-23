# HybridLWEAttack

Our "estimator_hybrid.py" estimates the bit-security, given parameter algorithm.
We uses functions in "estimator.py" made by Albrecht [Alb17]
and the high-level structure of our estimator "estimator_hybrid.py" is from "estimator.py".

### How to Run

We recommend to use our **MITM_estimator** function in "estimator_hybrid.py". This function estimates the cost of ring operations, given a bounded memory capacity 2^80.

EXAMPLE :\
    * sage: from estimator_hybrid import *\
    * sage: n = 8192\
    * sage: q = next_prime(2^125)\
    * sage: alpha = 8 / q\
    * sage: MITM_estimator(n, alpha, q, 11, 13)