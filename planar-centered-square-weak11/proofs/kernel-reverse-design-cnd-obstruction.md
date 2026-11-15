# Reverse-designed diamond potential: strict CND obstruction

Date: 2026-11-15 (Asia/Shanghai). Author workspace: `pcw-sol-kernel-certificates`.
Status: analytic obstruction; no new upper bound for the maximal operator.

## 1. Exact object and Fourier transform

Work in the diamond coordinates of `kernel-interface-audit.md`, so the unit averaging set is D={|u|+|v|≤1} and a feasible comparison kernel has cost (1/2)∫K. For 0<p<2 define U_p(u,v)=(|u|+|v|)^(-p). The Fourier convention is Ff(ξ)=∫exp(-i x·ξ)f(x)dx. Put a=|ξ_1|, b=|ξ_2|.

In tempered distributions, and as a positive function away from ξ=0,

    F U_p(ξ) = C_p (b^p-a^p)/(b²-a²),
    C_p = 2π/[Γ(p) sin(πp/2)] > 0.                 (1)

At a=b>0 the ratio is (p/2)a^(p-2), and at a=0<b it is b^(p-2). In particular the reciprocal, up to the harmless positive factor C_p, is

    φ_p(a,b)=(b²-a²)/(b^p-a^p),
    φ_p(a,0)=a^(2-p),  φ_p(0,b)=b^(2-p).           (2)

Proof of (1). The Laplace identity gives U_p=Γ(p)^(-1)∫_0^∞t^(p-1)exp(-t(|u|+|v|))dt, whose elementary Fourier factors are 2t/(t²+a²) and 2t/(t²+b²). Consequently

    F U_p = 4/Γ(p) ∫_0^∞ t^(p+1)/[(t²+a²)(t²+b²)]dt.

Substituting s=t² and q=p/2, then using

    s^q/[(s+A)(s+B)]
      = [B s^(q-1)/(s+B)-A s^(q-1)/(s+A)]/(B-A),
    ∫_0^∞ s^(q-1)/(s+A) ds = π A^(q-1)/sin(πq),

gives (1); all the latter integrals converge since 0<q<1. The limits cover a=b and one coordinate zero. The original distribution interchange is legitimate: for a Schwartz test, the absolute Laplace pairing is O(t^(p-1)) at t=0 and O(t^(p-3)) at infinity. The transform is locally integrable at ξ=0 because it is homogeneous of degree p-2 with p>0. This identifies the distribution with the displayed function, not merely its point samples.

## 2. A four-point necessary inequality excludes every 0<p<1

A real, even, continuous negative definite (CND) Markov symbol φ with φ(0)=0 must satisfy

    φ(x+y)+φ(x-y) ≤ 2φ(x)+2φ(y).                  (3)

Indeed use points (0,x,y,x+y) and coefficients (1,-1,-1,1) in the definition Σ_ij c_i c_j φ(z_i-z_j)≤0. The resulting quadratic form is exactly Q=2[φ(x+y)+φ(x-y)-2φ(x)-2φ(y)]≤0.

Take x=(a,0), y=(0,1), 0<a<1. Coordinate-sign invariance makes the two diagonals equal. For (2), direct cancellation gives the exact identity

    φ_p(a,1)-φ_p(a,0)-φ_p(0,1)
      = (a^p-a^(2-p))/(1-a^p).                   (4)

If 0<p<1, then a^p>a^(2-p), and the denominator is positive. Thus (3) fails strictly for every such a. Positive rescaling by C_p^(-1) does not change the sign. The reciprocal Fourier transform of U_p cannot be a CND Markov symbol for any 0<p<1. This is stronger than a sampling failure: the four-point contradiction is exact and covers the full interval of p analytically. It requires neither a differentiability theorem for stable spectra nor any assumption that the spectrum has a density.

A small rational witness is p=1/2, a=1/4. The normalized values are

    φ_p(a,1)=15/8,  φ_p(a,0)=1/8,  φ_p(0,1)=1.

The defect is 3/4, and the complete CND quadratic form is 3>0. This is an exact finite certificate.

For p=1, φ_1(a,b)=a+b and the unnormalized reciprocal is (|ξ_1|+|ξ_2|)/(2π), the valid symbol of two independent symmetric Cauchy coordinates. No assertion about 1<p<2 is needed below.

## 3. Independent integral and radius optimization check

For R>1 consider K=c(U_p-R^(-p))_+. It is integrable exactly in the present range p<2. Its minimal scalar enforcing K≥1 on all of D is c=(1-R^(-p))^(-1). Diamond polar integration, using |{||x||_1<r}|=2r², gives

    J_p(R)=(1/2)∫K
       = c p R^(2-p)/(2-p)
       = p R²/[(2-p)(R^p-1)].                   (5)

Writing T=R^p>1, the logarithmic derivative of T^(2/p)/(T-1) is 2/(pT)-1/(T-1). Hence the unique minimizer is T=2/(2-p), and

    inf_(R>1) J_p(R) = [2/(2-p)]^(2/p).         (6)

Thus the integral/radius formula proposed by Astra is correct. The boundary at ||x||_1=R contributes no extra integral; K vanishes continuously there. Equations (5)-(6) certify only C1/C2/C3 of this scalar candidate; generator feasibility remains a separate condition.

The cost in (6) is strictly increasing in p. Set t=p/2. Its logarithm is

    -log(1-t)/t = ∫_0^1 (1-ts)^(-1) ds,

whose integrand is strictly increasing in t for s>0. Its p↓0 limit is e, but every parameter in that potentially attractive branch 0<p<1 has been excluded by (4). For p≥1, the formal optimized cost is at least 4, with equality at p=1 (R=2,c=2). Therefore this reciprocal-potential scalar family cannot produce a cost below 4 through a CND generator, irrespective of the unsettled feasibility of p>1. In particular it cannot improve the inherited 3.590865312933 general upper bound.

This is an obstruction to the specified comparison-kernel route. It is neither a lower bound on C₂ nor a proof that arbitrary generators or arbitrary kernels have cost at least 4.

## 4. Convex mixtures have the same cost floor

The cost obstruction extends from a single hinge to every finite convex nondecreasing F:[0,∞)→[0,∞) with F(0)=0 for which K=F(U_p) is in L¹. No generator assumption is needed for the following purely integral statement.

The standard one-dimensional convex representation is

    F(s)=a s+∫_(0,∞) (s-t)_+ μ(dt),
    a=F'_+(0)≥0,  μ≥0.                               (7)

The measure is locally finite inside (0,∞); it may have infinite total mass. A nonzero a is impossible under L¹: F(U_p)≥a U_p and ∫_(||x||_1>1)U_p=∞ for p<2. Thus a=0. Write β=(2-p)/p>0. Every hinge H_t=(U_p-t)_+, t>0, has

    J(H_t)=(1/2)∫H_t=(p/(2-p)) t^(-β)=(1/β)t^(-β),
    H_t|_(∂D)=(1-t)_+.                              (8)

For 0<t<1, the ratio J(H_t)/(1-t) has its unique minimum at

    t=β/(β+1)=1-p/2,

and that minimum is exactly C(p)=[2/(2-p)]^(2/p). For t≥1 the following right side vanishes. Hence for every t>0,

    J(H_t) ≥ C(p)(1-t)_+.                           (9)

Tonelli, justified by nonnegative integrands, now gives

    J(F(U_p))=∫J(H_t)μ(dt) ≥ C(p)∫(1-t)_+μ(dt)
              =C(p)F(1).                           (10)

Since F is continuous at 1 and nondecreasing, K≥1 almost everywhere on D implies F(1)≥1. Therefore every such convex-composition majorant has J≥C(p), with equality achieved by the single optimized hinge. In particular, after excluding p<1 by the CND contradiction, any convex-composition route based on a reciprocal Markov potential U_p has cost at least 4.

At p=1 the elementary hinge inequality is particularly transparent:

    J(H_t)=1/t,
    1/t-4(1-t)=(2t-1)²/t≥0  for 0<t<1.

Equality is at t=1/2; normalization uses K=2(U_1-1/2)_+ and gives J=4. This rules out improvement through arbitrary convex mixtures and convex truncations in the stated family. It does not constrain an independently added nonradial S or nonconvex transformations.

## 5. Actual computation and dependency scope

`computation/kernel_reverse_design_cnd_check.py` uses only Python's exact Fraction arithmetic. It checks the full four-point quadratic form, independently checks identity (4), and checks the scalar integral/radius normalization at p=1 and p=1/2. The deterministic rational p list is 1/4,1/2,3/4,1,5/4,3/2,7/4; a=2^(-denominator(p)) makes every involved power rational. No random seed, floating tolerance, LP, large matrix, MATLAB, or new package is involved. The p>1 entries only record that this necessary four-point test no longer rejects them; they do not claim CND sufficiency.

The output records the actual interpreter/version, Git base, script hash, exact inputs, all quadratic forms, source-file hashes, and the CND test's logical limits. Reproduction, from the topic directory:

    python3 -B computation/kernel_reverse_design_cnd_check.py

The computation is a finite arithmetic cross-check. The proofs of (1)-(6) establish the continuum claims. R9's 1182-column certificate was read for scope only and was not recomputed. Frozen files were not changed. No commit/push was made in this delegated workspace; root owns integration.
