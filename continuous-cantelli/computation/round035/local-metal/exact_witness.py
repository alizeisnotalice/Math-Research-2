"""Exact finite-core density exclusions for two infinitely flat profile shapes."""
from fractions import Fraction as Q
from math import factorial


m, r, eps, sigma = Q(2), Q(1, 8), Q(1, 2**32), Q(1, 2**16)
target = m+eps
normal = Q(14, 25)
relative = Q(999, 1000)
assert sigma*sigma == eps and 3*sigma < r/8
assert r*r < m
assert normal**2 * 2 * Q(22, 7) < target

# v(u) = r² exp[(1-r/u)/4] for u>0. For u<=3σ,
# (r/u - 1)/4 > 512, so exp[(1-r/u)/4] < 8!/512^8.
assert (r/(3*sigma)-1)/4 > 512
v_ridge_upper = r*r*Q(factorial(8), 512**8)
assert relative**2*(1+v_ridge_upper/eps) < 1

# In the broad interval t=u/r∈[1/8,1], the derivative of
# log(t² exp((1/t-1)/4))=(2/t)-(1/(4*t²)) is nonnegative.
# The ratio is ≤1, hence u²≤v(u)≤r². Use e^-1/2>3/5.
assert relative**2*(1+eps/r**2) < 1


def lower_exp(x, degree):
    return sum((-x)**k/factorial(k) for k in range(degree+1))


assert lower_exp(Q(1, 2), 3) > Q(3, 5)
step = Q(1, 8)
near = 2*relative*step*sum(lower_exp((j*step)**2/2,31) for j in range(1,25))
broad_one_side = relative*Q(3, 5)*Q(7, 8)

# log(2) < 3/4 follows from the positive quadratic partial sum for exp(3/4).
assert 1 + Q(3, 4) + Q(3, 4)**2/2 > 2
for mode, a, broad_sides, source_loss, source_bound, threshold in (
    ("symmetric", Q(3, 2), 2, Q(4, 7), Q(14,25), Q(16,15)),
    ("right plateau", Q(5, 4), 1, Q(2, 5), Q(2,3), Q(27,25)),
):
    assert a*a > Q(3, 2)  # stronger than a²>2*log(2)
    loss = a*a*(1-1/target)/2
    # Narrow central region uses the same bound for both sides, including any plateau.
    assert loss+a*3*sigma+Q(9, 2)*eps < source_loss
    if broad_sides == 2:
        # Equal variances on the two sides permit the cosh(au) >= 1 sum bound.
        assert loss+r*r/2 < source_loss
    # In the plateau case, the broad interval is on x=a-u with source exponent ≥-a²/2.
    assert 0 < a-r < a and loss < source_loss
    assert lower_exp(source_loss, 3) > source_bound
    lower = normal*source_bound*(near+broad_sides*broad_one_side)
    assert lower > threshold
    print(f"{mode}, a={a}: lower ratio {float(lower):.15f} > {threshold} (exact comparison)")

print(f"common narrow integral lower bound = {float(near):.15f}")
print(f"broad integral lower bound per active side = {float(broad_one_side):.15f}")
