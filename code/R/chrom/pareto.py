import numpy as np
import pandas as pd
# from scipy.special import binom
import math
# from mpmath import *
from mpmath import mp
# import timeit
import time


# Sieve of Eratosthenes
def get_first_primes(n):
    is_prime = [True]*(n+1)
    p = 2
    while p**2 <= n:
        if is_prime[p]:
            for i in range(p**2, n+1, p):
                is_prime[i] = False
        p += 1
    return [p for p in range(2, n) if is_prime[p]]


# Compute the least common multiple of 1,2,..,n
def lcm_first(n):
    primes = get_first_primes(n)
    lcm = 1
    for p in primes:
        lcm *= p**int(math.log(n)//math.log(p))
    return lcm


def pareto_E_Z1Z2_python(k, n, dps = 100):
    start = time.time()
    mp.dps = dps  # set significant digits
    n_factorial_vec = [math.factorial(i) for i in range(n-1)]  # prepare in advance
    max_denom_lcm = lcm_first(n)**(3*k)

    e_k_n_int = int(0)
    ctr = 0
    for a in range(n-1):
        if a % 10 == 0:
            print("Run a: " + str(a) + " out of " + str(n-2))
        for b in range(n-1-a):
            for c in range(n-1-a-b):
                d = n-2-a-b-c
                numerator = (n_factorial_vec[n-2] // (n_factorial_vec[a]*n_factorial_vec[b]*n_factorial_vec[c]*n_factorial_vec[d])) * (-1)**(a+b)  * ((a+b+2*c+2)**k - 2*(b+c+1)**k)
                denominator = ((a+c+1)*(b+c+1)*(a+b+c+2))**k
                e_k_n_int += numerator * (max_denom_lcm // denominator)
                ctr += 1
#    return 0.1

    print("k=" + str(k) + ", n=" + str(n) + ", Run time: " + str(time.time()-start))
    return mpf(e_k_n_int) / mpf(max_denom_lcm)
