## Attack Implementation

"mitm.py" is an attack code for Sage.
This is just for a proof-of-concept of the hybrid attack, and hence it is definitely not optimized.

### EXAMPLE :

    > sage: from Mitm import *
    > sage: n = 45; q = next_prime(2^14); h = 6; beta = 20; k = 30; tau = 30;
    > sage: A, c, u, s = gen_instance(n, q, h, m = tau * n)
    
    (A, c) is n-dim LWE samples (with secret s) / (A, u) is uniform samples
    
    > sage: A_k, c_k, u_k, B = dim_error_tradeoff(A, c, u, beta, h, k)
    
    (A_k, c_k) is k-dim LWE samples (with secret s[-k:]) / (A_k, u_k) is uniform samples. 
    
    > sage: Mitm_on_LWE(A_k, c_k, u_k, B, h)
    
    Table size = 34281

    ** Mitm on (A_k, c_k) ** 

    Number of noisy searches = 193
     - Input is LWE with secret
    (0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)

    ** Mitm on (A_k, u_k) **

    Number of noisy searches = 34281
     - Input is uniform
     
    > sage: s[-k:]
    
    (0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)


