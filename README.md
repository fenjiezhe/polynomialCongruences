# Implementation of Coppersmith method

## Algorithm 1(plain)
**Input**:

RSA-type composite N 

h>=2

integer-coefficient polynomial f(x) with degree delta and coefficients in {0,1,...,N-1}, where 2<delta+1<(logN)/2

**Output**: all integer roots of polynomial congruence f(x)=0 mod N, the absolute value of the roots is bounded by power(N,1/delta)
