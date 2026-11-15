# R55 coordinate-Cauchy smoothing: independent formula and screen audit

2026-11-15. This attachment establishes analytic identities and conditional certification tools. It does not certify a new general C₂ upper bound. No frozen R51 files were changed, no MATLAB process was started, and no commit/push was made.

## 1. Normalization and all degenerate cases

Let ρ(x)=|x₁|+|x₂|, U=1/ρ and
p_t(x)=t²/[π²(t²+x₁²)(t²+x₂²)], t>0. The coordinate Cauchy generator is A₀=−|D₁|−|D₂| and A₀U=−2πδ₀. Since U=2π∫₀^∞p_s ds and p_s*p_t=p_(s+t), nonnegative Fubini gives

V_t:=U*p_t=(2/π)∫_t^∞ s²/[(s²+a²)(s²+b²)] ds,

where a=|x₁|, b=|x₂|. For a≠b this equals

V_t=(2/π)[b atan(b/t)−a atan(a/t)]/(b²−a²).

At a=b>0 its continuous value is

V_t=(1/π)[atan(a/t)/a+t/(t²+a²)].

At a=0<b it is 2 atan(b/t)/(πb). At a=b=0 it is 2/(πt). Thus V_t is finite and smooth everywhere; its apparent diagonal/axis divisions are removable. No singular contribution is discarded. Equivalently,

V_t=(2/(πt))∫₀¹[(1+(au/t)²)(1+(bu/t)²)]⁻¹ du.

For b≥a≥0, a+b>0, set h=t²+ab, z=t(b−a)/h and atanc(z)=atan(z)/z with atanc(0)=1. The cancellation-free closed expression is

V_t=2/[π(a+b)]·[atan(b/t)+(at/h)atanc(z)].

Indeed b atan(b/t)−a atan(a/t)=(b−a)atan(b/t)+a atan(t(b−a)/h). This avoids a small difference of nearly equal terms. An implementation should additionally guard extreme floating overflow/underflow; its algebraic identity is exact, its floating output is not directed.

For Q_t=U−V_t=2π∫₀^t p_s ds≥0 and a+b>0,

Q_t=2/[π(a+b)]·[atan(t/b)−(at/h)atanc(z)].

For a=0<b this reduces to 2 atan(t/b)/(πb). At a=b>0 it is [atan(t/a)/a−t/(t²+a²)]/π. At the origin Q_t has the original integrable positive singularity. If both a,b are large relative to t, the last closed form subtracts nearly equal terms; use its positive integral instead:

Q_t=(2t/π)∫₀¹ (tu)²/[((tu)²+a²)((tu)²+b²)] du.

## 2. Exact conditions for the heavy-tail C4 interface

Take a finite nonnegative atomic measure σ=Σa_jδ_(z_j), with A=Σa_j>0, Σa_jz_j=0, and common t>0. Let ν=σ*p_t, W_t=A U−U*ν and K_t=(W_t)_+. Then

W_t=W_σ+Σa_jQ_t(·−z_j),  W_σ=A U−U*σ.

The compact zero-mean tail proof gives W_σ∈L¹ and ∫W_σ=0; Q_t≥0 has ∫Q_t=2πt. Consequently W_t∈L¹ and

∫W_t=2πAt,    ∫K_t=(∥W_t∥₁+2πAt)/2.

The finite-translation decomposition also proves the generator identity without interchanging a noncompact ν with a singular U:

A₀W_t=2π(ν−Aδ₀).

Mollify W_t with a nonnegative compact mollifier. Away from its shrinking support at zero, A₀W_(t,ε)=2πν*η_ε≥0. For smooth convex increasing approximations Φ of positive part with Φ(0)=0 and 0≤Φ′≤1, the symmetric jump inequality gives A₀Φ(W_(t,ε))≥Φ′(W_(t,ε))A₀W_(t,ε)≥0 there. First pass ε→0 in L¹, then Φ→positive part in L¹. Pairings are continuous because A₀φ is bounded for compact smooth φ. Therefore A₀K_t≥0 on R²\{0} in distributions, including axes and positive-part interfaces.

Even σ gives even K_t. Nonnegative weights, matching singular mass A, zero finite-position mean, L¹, and the generator identity are substantive hypotheses. The heavy-tail measure ν itself does not possess the ordinary first moment used here. A nonzero finite-position mean generally destroys L¹, and unequal coordinate Cauchy times are not covered by this interface. ν may have density on D; separation from D is unnecessary.

For comparison one still needs m₀>0 with W_t≥m₀ almost everywhere on all D and a directed bound B≥∫K_t. Only then C₂≤B/(2m₀). Numeric sampling and full-domain adaptive quadrature do not meet those requirements. The signed integral alone gives ∫K_t≥2πAt, hence candidate cost at least πAt/m.

## 3. Full-D minimum and an elementary rigorous enclosure route

Differentiation of the positive V integral gives

|∂_aV_t|≤(4/π)∫_t^∞ a/(s²+a²)² ds≤3√3/(8πt²)=C/t²,

because max_(a≥0) a/(s²+a²)²=9/(16√3 s³). The same bound holds for ∂_b. Along an L¹-unit ray x=rv, 0<r≤1,

dW_t(rv)/dr≤−A/r²+AC/t².

Thus t²≥C=3√3/(8π), numerically t≥0.4546958717463487, suffices to reduce the full-D infimum to its boundary. The decimal is only explanatory; the exact inequality is the criterion. This condition covers the current screen's t=1, not .025, .1 or .3. It gives a boundary reduction, not the boundary minimum itself.

For any t, V_t≤2/(πt) gives W_t≥A/r−2A/(πt), so an inner diamond can be excluded once that lower bound exceeds the target m. On a box B away from zero, let r_max=max_Bρ(x), and let d_j,i be the distance of source coordinate z_j,i to the coordinate interval B_i. Since V_t decreases in each absolute coordinate,

inf_BW_t≥A/r_max−Σa_jV_t(d_j,1,d_j,2).

A finite box cover of D outside the excluded inner diamond, with directed evaluation, can certify m₀. This remains a proposed implementation until an actual candidate warrants it.

## 4. Stable far field and a genuine infinite-tail bound

For x₁,x₂≥0 set r=x₁+x₂ and d_j=ρ(x−z_j)−r. Exact algebra gives

W_σ=D₁/r²−Σa_j d_j²/[r²(r+d_j)],  D₁=Σa_jd_j.

For coordinate-reflection symmetric σ,

D₁=Σ_(z_j,1>0)2a_j max(z_j,1−x₁,0)+Σ_(z_j,2>0)2a_j max(z_j,2−x₂,0).

This computes D₁ as a nonnegative sum; outside both source coordinate extents D₁=0. It avoids cancellation among O(r⁻²) summands that hide an O(r⁻³) corner result. Add the positive Q integrals to recover W_t. Cancellation may still occur where the actual W_t is close to zero; no floating formula removes mathematical sign sensitivity there.

Here is a conservative exact full-tail budget. Put L=max_j max(|z_j,1|,|z_j,2|), S=max_jρ(z_j). For R>max(L,S),

∫_(ρ>R)K_t≤A[8L log(R/(R−S))+8t²/(R−S)].                 (T)

Proof: in the corner |x₁|>L, |x₂|>L, zero mean and convexity imply W_σ≤0. Everywhere with r>R>S, |W_σ|≤AS/[r(r−S)]. On each of the two cross strips, integration in its long coordinate yields at most 4AL log(R/(R−S)); adding strips gives the first term. Also ∫_(ρ>H)Q_t=2π∫₀^tP(ρ(X_s)>H)ds. A union bound on the two coordinates gives P≤(4/π)atan(2s/H)≤8s/(πH), whence this mass is ≤8t²/H. Every shift z_j has ρ(z_j)≤S, so use H=R−S and sum its masses. Finally (W_σ+Σa_jQ_j)_+≤(W_σ)_++Σa_jQ_j. This proves (T), including all cross strips and heavy tails. It may be conservative and has not been used to claim any numerical cost enclosure.

## 5. Independent read-only MATLAB screen audit

Audited `cauchy_heavy_tail_fixed_screen.m`, SHA256 `85b46ef12e307e72383876ad73cd29f583fb9616dcfba9f30ce40cac600272e1`, 4848 bytes. The final version includes the reflection-paired D₁ and quadratic remainder and asserts max source ρ<64, protecting the far branch from a Q origin singularity.

For a D4 source, the entire plane has eight congruent wedges x≥y≥0, so a half integral equals four times that wedge integral. Splitting at y=L is exhaustive. The maps x=y+Lu/(1−u), 0≤y≤L and y=L+Lv/(1−v) for the other piece have the stated Jacobians L/(1−u)² and L²/[(1−u)²(1−v)²]. Both range to infinity; there is no finite cutoff. The code comment saying “full quadrant is eight” is a wording slip: eight wedges cover the full plane; the actual factor4 is correct.

The grid and fmincon constraints sample the full triangular D4 portion of D, not merely its boundary; finite samples/local optimization do not certify its minimum. The full unbounded quadrature uses float tolerance, and shortQ's far branch is 8-point Gauss, so neither cost nor tail is a rigorous enclosure. The final output correctly records these limits. Exact D4, nonnegative weights and zero position mean of this actual input were checked independently below; the screen routine itself assumes its orbit input has those properties.

## 6. Actual bounded Python run and source record

`computation/cauchy_smoothed_potential_independent.py` actually ran with Python3.13.18, SymPy1.14.0, mpmath1.3.0 at 100 decimal digits, NumPy2.5.4, one process, no random inputs. Six symbolic identities passed, including the derivative primitive, axis/diagonal limits and compact cancellation identity. Seven independent positive-integral checks matched the closed V expressions to relative discrepancy below 10⁻⁸⁵. Ten floating stress cases gave stable-V max relative error 2.98e−16; the raw divided difference lost about 0.002716 relative accuracy on its worst selected near-diagonal case. Independently generated 8-point Gauss checks on six far-Q samples gave max relative error 4.01e−16. These are nondirected numerical tests, not uniform error certificates.

The original fixed LP JSON SHA256 is `8cb643aab0b7b732444f54abf74379afcc9130128d85c5924f5e76eabc1764d6` (24312 bytes). Its coefficients are interpreted as the exact rationals of binary64-parsed JSON values, matching MATLAB jsondecode; dividing by each orbit size (4 or 8) is exact. The exact positive orbit data are archived in `computation/inputs/cauchy-smoothed-potential-independent/fixed-positive-orbits.json`, SHA256 `99ca2dc483ff95e33b2b4a4ff1b99be2b891e6dc10f86fe8a8a315ffd06d1f44` (5913 bytes). Exact Fraction checks give 76 positive atoms, D4 symmetry, zero mean, A=1077997594621154451/576460752303423488, L=S=16<64. At five selected rational far points the reflection-paired first moment and compact identity agree exactly. Saved output also contains all 76 exact atoms, input origins and source hashes.

Read source heavy-interface note SHA256 `5502cf9c0b2617a9e6e993415a523e8a61bac3308f640b68015ce28f5b3bcc64`; read source compact/Kato note SHA256 `be663fccd68ec7e6fca91de29935fd34d740bb56e83b8b2501f3c2c0d310a339`. The C4 argument and formulas above were independently checked; the script cross-checks only scalar identities and fixed finite data, not the general theorem by sampling.

Reproduce from the topic directory:

    /Users/zzh/.venvs/main/bin/python3 computation/cauchy_smoothed_potential_independent.py

The final script/input/output hashes and whitelist are in the adjacent delivery manifest. Remaining mathematical work is a certified full-D m and full-plane positive-part integral for an actually useful candidate; the general exact constant remains unresolved.
