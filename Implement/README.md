## Attack Implementation

"mitm.py" is an attack code for Sage.
This is just for a proof-of-concept of the hybrid attack, and hence it is definitely not optimized.

### How to run

The most high-level function is **hybrid_mitm**.
It takes the following input

    n, q, alpha (=8/q default) 		: LWE parameters
    h 					: the hamming weight of secret vector
    beta 					: BKZ blocksize
    k					: LWE dim after Dim-Error tradeoff
    tau (= 30 default)			: LWE sample number after Dim-Error tradeoff
    ell (= h default)			: MITM parameter (hamming weight bound)

and generate LWE sample (A, **c**) and Uniform vector **u**,
and perform the hybrid MITM attack on (A, **c**), and then on (A, **u**).
  
EXAMPLE :\
    * sage: from Mitm import *\
    * sage: n = 50; q = next_prime(2^18); h = 8; beta = 20; k = 30;\
    * sage: hybrid_mitm(n, q, h, beta, k)
    
    Performing Dimension-error trade-off . . .

    New k-dim samples have error bound = 405.96, where q = 16411

    Start MITM on the k-dim samples . . .

    Table size = 451
    
    ** Mitm with (YA_2, Yc) ** 
    Number of noisy searches = 48
     - Input is LWE

    ** Mitm with (YA_2, Yu) **
    Number of noisy searches = 451
     - Input is uniform    


