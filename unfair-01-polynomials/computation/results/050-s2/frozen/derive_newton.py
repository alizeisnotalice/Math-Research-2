"""Exact Newton identities and modulus-envelope data for the S={2} core."""
import json
import sympy as sp

x = sp.symbols("x")
b, a, c, d, e, f = sp.symbols("b a c d e f")
T, N = sp.symbols("T N")
# If P(x)=prod(1-z_i*x), then its roots' reciprocals z_i have
# monic polynomial y^7+b*y^6+a*y^5+c*y^4+d*y^3+e*y^2+f*y+1.
coeff = [1, b, a, c, d, e, f, 1]
s = {0: 7}
for k in range(1, 8):
    value = 0
    for j in range(1, k):
        value += coeff[j] * s[k-j]
    value += k * coeff[k]
    s[k] = sp.expand(-value)


def reduce_a(expr):
    """Return p,q with expr == p+q*a modulo a^2-T*a+N."""
    rem = sp.rem(sp.Poly(sp.expand(expr), a), sp.Poly(a**2 - T*a + N, a))
    return [sp.expand(rem.nth(0)), sp.expand(rem.nth(1))]

envelope = sp.Poly((1 + 2*x)**3 * (1 + x/sp.Integer(2))**3 * (1 + x), x)
out = {
    "newton_sums": {str(k): str(s[k]) for k in range(1, 8)},
    "reduced_p_plus_q_a": {str(k): [str(z) for z in reduce_a(s[k])] for k in range(2, 8)},
    "power_sum_bounds": {str(k): str(3*2**k + sp.Rational(3, 2**k) + 1) for k in range(1, 8)},
    "envelope_coefficients": [str(envelope.nth(k)) for k in range(8)],
    "envelope": str(envelope.as_expr()),
}
print(json.dumps(out, indent=2))
