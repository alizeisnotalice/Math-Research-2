"""A zero-radius smooth quartic jet raises the core ratio across one."""
from fractions import Fraction as Q
from math import factorial

m,a,r,c,eps,sigma = Q(5,4),Q(3,2),Q(1,8),Q(16),Q(1,2**40),Q(1,2**20)
target = m+eps
loss = a*a*(1-1/target)/2
assert eps == sigma**2 and c*r*r == Q(1,4)
assert 5*c*r**4 < m < a*a  # a² > M log(M)/(M-1).

def exp_partial(x,n):
    return sum(x**k/factorial(k) for k in range(n+1))

# Corrected: source pair 2*exp(-a²/2-u²/2)*cosh(a*u).
assert Q(44,100)**2*2*Q(22,7) < target
assert loss+r*r/2 < Q(7,30)
assert exp_partial(-Q(7,30),3) > Q(79,100)
assert Q(999,1000)**2*(1+5*c*(3*sigma)**4/eps) < 1
assert 3*sigma < r/2
near = 2*Q(999,1000)*Q(1,8)*sum(exp_partial(-Q(j,8)**2/2,31) for j in range(1,25))
# On t=u/r in [1/2,1], H(32*t^2) >= 8/9 > 7/8.
assert Q(8,9)>Q(7,8)
# Thus (9/8)r²t⁴ <= v <= (5/4)r²t⁴. Substitute z=1/t in the integral.
assert Q(89,100)**2*(Q(5,4)+eps/(r*r*Q(1,2)**4)) < 1
step = Q(1,64)
broad = 2*Q(89,100)*step*sum(exp_partial(-Q(4,9)*(1+j*step)**2,23) for j in range(1,65))
lower = Q(44,100)*Q(79,100)*(near+broad)
assert lower > Q(21,20)

# Baseline: cosh(a*u) <= exp(a²r²/2) <= 1/(1-a²r²/2).
assert Q(45,100)**2*2*Q(314,100) > target
assert loss > Q(9,40)
assert exp_partial(-Q(9,40),2) < Q(801,1000)
pair_upper = 1/(1-a*a*r*r/2)
u1,u2 = Q(1,2**13),Q(1,2**9)
variance_factor = Q(257,256)
# Region 1: enlarge to the full Gaussian integral.
assert eps+c*u1**4 == variance_factor*eps
near_upper = Q(63,25)
assert near_upper**2 > 2*Q(22,7)*variance_factor
# Region 2: (eps+v)/u² <= 2^-13, so exponent <= -4096.
assert eps/u1**2+c*u2**2 == Q(1,2**13)
mid_upper = 2*u2/sigma*Q(factorial(4),4096**4)
assert mid_upper < Q(1,100000)
# Region 3: eps+v <= (257/256)*16*u^4. Substitute z=1/(4u) and
# bound 2*int_2^infty exp[-z²/(2*factor)] dz by factor*exp[-2/factor].
assert 1+eps/(c*u2**4) == variance_factor
tail_upper = variance_factor/exp_partial(2/variance_factor,7)
upper = Q(45,100)*Q(801,1000)*pair_upper*(near_upper+mid_upper+tail_upper)
assert upper < Q(49,50)
print(f"corrected core ratio > {float(lower):.15f} > 21/20 (exact comparison)")
print(f"matched quartic ratio < {float(upper):.15f} < 49/50 (exact comparison)")

# Scaled limit v_r=r² Phi(u/r), eps=r^4. Its broad coefficient for the
# corrected profile exceeds 1/4 by the same [1/2,1] calculation at eps=0.
assert Q(39,100)**2*2*Q(22,7) < 1 and Q(39,100)*broad > Q(1,4)
# Matched quartic: broad coefficient 2*normal_bar(2) < exp(-2)/sqrt(2pi).
assert Q(41,100)**2*6 > 1 and exp_partial(Q(2),4) == 7
assert Q(41,700) < Q(1,16)
print("scaled corrected and baseline broad factors >1/4 and <1/16 respectively")
