# Exact calibration of finitely many probes with a linear output defect on the true Gaussian input

For every finite integer N≥2 and every finite list of distinct nonzero real frequencies z_i, there is a family of smooth bounded positive time graphs f_η for all sufficiently small η>0, equal to 1 outside a fixed compact set. With X standard normal, put M_f(z)=E exp(zX−z²f(X)/2). The graphs satisfy
```
M_{f_η}(z_i)=1 for every i,       M_{f_η}^{(k)}(0)=0 for 1≤k≤N,
f_η=1−η on [−1,1],               ||f_η−1||_(C^j) ≤ C_j η for each fixed j.
```
In particular E f_η=1. There is a common H_η=1+Cη≥sup f_η for which the Gaussian output defect obeys
```
liminf_{η↓0} η^(−1) W₁(law(X+sqrt(H_η−f_η(X))G), N(0,H_η)) ≥ 1/√(2π).
```
Constants and smallness thresholds may depend on the fixed finite probe list. This is a certified non-exact family for the original Gaussian input. Its absolute valley depth tends to zero; the theorem supplies no continuous nonconstant limit satisfying a growing list of probes.

**Finite-dimensional calibration.** Take a smooth compactly supported h₀ equal to −1 on [−1,1] with support in (−2,2). On (2,3) consider the functions v_j consisting of probabilists' Hermite polynomials He₀,…,He_(N−2) and exponentials exp(z_i x−z_i²/2). They are linearly independent on any open subinterval: analyticity extends a relation to the line, and repeated differentiation beyond the polynomial degree then gives independence of exponentials by their Vandermonde determinant. Take a smooth nonnegative compactly supported cutoff ψ inside (2,3) positive on an open interval, and set h_j=ψ v_j for these correcting functions (indexed starting at one, independently of h₀). The Gram matrix ∫ψ v_i v_j dγ is positive definite.

Write f=1+η h₀+Σa_j h_j. The constraints are real analytic in these finite parameters near zero, by boundedness of the perturbations and Gaussian exponential moments. The first variation of M_f(z) in direction h at f=1 equals −(z²/2)∫h exp(zx−z²/2)dγ. Thus the variations of the jets of orders k=2,…,N are nonzero multiples of inner products with He_(k−2). Evaluations at z_i give nonzero multiples of inner products with the listed exponentials. Order one always vanishes. The correcting Jacobian is consequently an invertible diagonal scaling of the Gram matrix. The implicit function theorem gives real analytic a_j(η)=O(η), and all stated spatial bounds follow. Increasing C if necessary gives H_η>sup f_η and f_η positive.

**Exact lower bound on the output defect.** This step uses just the valley and the upper envelope, independent of how the probe constraints were attained. Take H=1+ε with ε=Cη and y=0. The potential U of round 206 obeys U(H,0)=E|X+sqrt(H−f(X))G|−E|N(0,H)|, whether or not f is exact. Start the reverse Brownian motion at (H,0), and stop at the graph or height 1−η. Remaining in (−L,L), 0<L≤1, for time ε+η gives a positive terminal value at least m=inf_{|x|≤L}∫_{1−η}¹p_s(x)ds. Elsewhere any negative terminal value must be at a graph hit with f>1, has magnitude at most ε/√(2π), and requires the path to leave (−L,L) before time ε. These statements follow from the explicit U of round 206, which is nonnegative at heights ≤1. Caloric optional stopping and the reflection principle thus give
```
W₁ ≥ U(H,0) ≥ m S(L,ε+η)−4ε Φ̄(L/√ε)/√(2π).
```
Now L=η^(1/4) gives m/η→p₁(0), S→1 and the negative term is o(η). The displayed finite-η inequality itself can be used without an asymptotic approximation, with S≥max(0,1−4Φ̄(L/√(ε+η))). It also applies, with the corresponding shallower valley and increased envelope, to a uniform neighborhood of these graphs, regardless of probe values.
