### Task

- Dimension을 높이는 경우, 너무 낮은 block size를 사용하면 bkz에서 오류 발생.
Infinite loop in babai 라는데, precision 문제로 생각됨.

- One BKZ & Iterative LLL 전략은 정말 잘 작동하는걸까?
적당히 큰 Param에 대해 실험하고 싶은데 자꾸 Infinite loop 오류가 나서 확인할수가 없고, 일단 Albrecht가 주장하는 n = 100, q = next_prime(2^23), h = 20, beta = 50에 대해 LLL = 256회를 수행하면 서로 다른 벡터가 별로 안 나옴

- 적당히 작은 Param에 대한 실험은 성공함!

### How to use

## Requirements

Run in Python

- Sage
- LWE-estimator (by Alb)
- fpylll

## Example

A, c, u, s = gen_instance(n, q, h, alpha, m)
''' 
Input Params
n :		LWE dim
q :		LWE mod
h :	 	HW(s)
alpha : 	8/q (Default)
m :		n (Default)

Output Params
(A, c) : 	LWE instance
u :		Uniform vector in GF(Z_q)
s :		Secret
'''

hybrid_mitm(A, c, u, beta, h, s, k, num_sample, ell, check_unif)
'''
Solves LWE with input (A,c) by Hybrid-Mitm strategy 
(This is void, and print the result by print)

Input Params
beta :		BKZ block size
k :			LWE dim after trade-off
num_sample :	Desired LWE samples after trade-off
ell :			Presumed HW(s_2)
check_unif :	Flag standing for applying attack also with (A,u)
'''
