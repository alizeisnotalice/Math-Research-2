# Independent audit: every subunit diamond singularity fails stable C4

Date: 2026-11-15 (Asia/Shanghai). Proposed by `/root/astra_kernel`; this file is the Sol certificate auditor's analytic verification. Status: strict obstruction, with the assumptions below. It is not a claim about all comparison kernels.

## 1. Statement and normalization

Let 0<p<1, 0<α<2, c>0, R>0, S∈C_c²(R²), and

    K(x)=c(||x||_1^(-p)-R^(-p))_+ + S(x).

Let σ be a nonzero finite nonnegative angular measure on S¹ whose support spans R². Use the symmetric homogeneous stable generator

    A_σ f(x)=∫_(S¹)∫_R
      [f(x+tθ)+f(x-tθ)-2f(x)] |t|^(-1-α) dt σ(dθ).  (1)

Any conventional positive multiplicative normalization is harmless. The spectrum need not have a density or lack atoms. The conclusion is

    A_σ K cannot be a nonnegative distribution on R²\{0}. (2)

The proof has two cases that exhaust all α. No C1 nonnegativity of K is used; the contradiction applies even when the smooth correction has arbitrary sign. Compact truncation and S imply K∈L¹. The Fourier symbol of -A_σ is a positive multiple of ∫|ξ·θ|^α σ(dθ), strictly positive for ξ≠0 by the spanning assumption.

## 2. Subcritical origin: p+α<2

Choose χ∈C_c^∞, 0≤χ≤1, χ=1 near 0, and put χ_ε(x)=χ(x/ε). For a fixed smooth compactly supported φ, write

    K=cU_p+V,  U_p=||x||_1^(-p),  V∈L∞(R²).

Boundedness of V holds because it equals -cR^(-p)+S inside the diamond and -cU_p+S outside. Self-adjointness, homogeneity, and scaling give

    |<A_σ(cU_p),χ_ε φ>| ≤ C ε^(2-p-α),
    |<A_σ V,χ_ε φ>| ≤ C ε^(2-α).                    (3)

For completeness, after x=εz the first pairing is ε^(2-p-α)<U_p,A_σ g_ε>, where g_ε(z)=χ(z)φ(εz) ranges in a bounded C_c² set. These pairings are uniformly bounded: near the origin U_p is integrable (p<2), and for the far translated test functions the t-tail is integrable. For the second pairing, ||A_σ g_ε||_1 is uniformly bounded by the C² translation estimate at |t|≤1 and the L¹ translation estimate at |t|>1; scaling then gives ε^(2-α). Finite σ suffices for both estimates. The use of bounded V does not require V itself to be in L¹.

Suppose A_σK≥0 off 0. For φ≥0 the test (1-χ_ε)φ is nonnegative and supported away from 0. As both exponents in (3) are positive, passing ε↓0 yields <A_σK,φ>≥0. Thus A_σK is globally a nonnegative distribution, hence a positive Radon measure ν.

Let χ_R(x)=χ_0(x/R), where χ_0≥0 is smooth, compactly supported and equals 1 on the unit ball. The exact stable scaling gives ||A_σχ_R||_∞≤C R^(-α). Since K∈L¹,

    0≤<ν,χ_R>=<K,A_σχ_R>≤||K||_1 C R^(-α)→0.

Every fixed compact set lies where χ_R=1 for sufficiently large R, so ν must vanish. Hence A_σK=0. The continuous Fourier transform of the L¹ function K vanishes for ξ≠0 because the multiplier there is strictly positive, and it also vanishes at 0 by continuity. This forces K=0, contradicting its positive singular leading term and bounded S. This proves (2) in the subcritical case.

This argument does not infer global positivity merely from pointwise values; estimate (3) explicitly excludes any hidden negative origin-supported distribution.

## 3. Complementary case: p+α≥2

Now α≥2-p>1. Fix 0<u₀<R, and let v↓0 with x_v=(u₀,v). Near (u₀,0), the singular term has expansion

    c((u+|v|)^(-p)-R^(-p))
      =F₀(u)-a(u)|v|+v²H(u,|v|),
    a(u)=cp u^(-p-1)>0,

with smooth F₀,a,H on an appropriate fixed local rectangle. The model cusp has the exact one-dimensional integral

    ∫_R [|v+t|-2|v|+|v-t|]|t|^(-1-α)dt
       =4/[α(α-1)] |v|^(1-α).                     (4)

Indeed its integrand vanishes for |t|≤|v| and equals 2(|t|-|v|) above that threshold. Scaling tθ₂ in (4) gives the leading term

    A_σ K(u₀,v)
      =-4cp u₀^(-p-1)/[α(α-1)] · M₂ · v^(1-α)+O(1),
    M₂=∫_(S¹)|θ₂|^α σ(dθ)>0.                     (5)

Here the O(1) is genuinely uniform over the angular spectrum, including atoms. Locally, subtract the fixed cusp -a(u₀)|v|. The F₀ and v²H remainders have uniformly bounded second differences O(t²). For the mixed cusp term, expand a(u₀±tθ₁)=a(u₀)±a'(u₀)tθ₁+O(t²). The first-order pair is bounded using

    ||v+tθ₂|-|v-tθ₂||≤2|tθ₂|,

and is also O(t²), uniformly in θ and small v. Integrating t^(1-α) near 0 is finite since α<2. S contributes O(1) by the same C² estimate.

For far jumps, the only singularity lies at the origin. For any fixed bounded t-interval and any θ, projection onto θ gives

    ||x_v+tθ||_1 ≥ |t+x_v·θ|.

Therefore the line integral of ||x_v+tθ||_1^(-p) is uniformly bounded, because p<1 and ∫_I|t+b|^(-p)dt is bounded uniformly over translations b of a bounded interval I. The jump weight is bounded away from t=0 in this far part. The remaining t-tail and the removed model cusp tail are uniformly integrable because α>1. This proves the O(1) in (5).

The expansion is uniform when u₀ ranges over a small compact interval inside (0,R). Its strictly negative leading coefficient and 1-α<0 imply a negative open strip with v>0 sufficiently small. On that strip the generator integral is well defined: K is locally C² there, and p<1 gives the uniform far-line integrability just proved. Truncation and dominated convergence identify the distribution with its continuous pointwise representative there. Consequently a nonnegative test supported in the strip has <A_σK,φ><0, contradicting C4.

The spanning assumption ensures M₂>0. In particular the claim does not discard an axis as a zero-area set: the contradiction occurs in an open set adjacent to it.

## 4. Targeted audit of Astra's frozen author draft

Read in full: `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-kernel/planar-centered-square-weak11/proofs/kernel-subcritical-cusp-obstruction.md`, SHA-256

    8d8b296824a2e545bde73ebb0c95b4de39f6a6a7e47392458482c3ac8b4eb3b0.

The author uses t∈(0,∞), whereas (1) here uses t∈R. Thus every generator and cusp coefficient here is exactly twice the author's; author (3.2) has coefficient 2/[α(α-1)] and (4) here has 4/[α(α-1)]. There is no sign or normalization discrepancy.

Targeted findings:

* Author §2.2's weighted L¹ tail estimate is uniform over θ and permits arbitrary finite singular angular measures. Moving the support by y=z±tθ produces |z|≥t/2 for the large-t tail, so the claimed p-weighted and unweighted tail integrals converge.
* Author §2's cutoff limit legitimately extends off-origin positive distribution to the whole plane. Its final large-cutoff argument does not presume finite total mass; positivity and the L¹ bound force every compact mass to vanish.
* Author §3.1's uniform O(t²) conclusion is correct. The explicit F₀(u)-a(u)|v|+v²H decomposition in §3 above supplies a direct proof. A generic O((s-u₀)²) bound alone should not be read as a bound on its second difference; smoothness plus the displayed mixed-cusp cancellation is the operative fact.
* Author §3's far-line bound is valid for p<1 and uniform over all line orientations, including directions carrying spectral atoms. The resulting negative open strip suffices to contradict the distribution condition.
* Author §4's Brownian endpoint also holds: for A=Q:D² with Q positive definite, away from the origin the line v=0 carries -2cp u^(-p-1)Q₂₂ δ_(v=0). The other local terms and S have bounded densities, so the negative line measure cannot be compensated.

No blocking mathematical gap was found in the frozen author statement under its explicit assumptions. The proof above is an independent expansion and error estimate, rather than reliance on the author's finite arithmetic PASS record.

The author then adopted the explicit local decomposition in §3.1. The revised author file was read at the revised section and its actual SHA-256 was checked as

    02ed862985cce2ce8bea7aac85b6c3cb396e5b9a1bad1f3f1246116f4d263e6d.

The new §3.1 now explicitly proves the uniform O(t²) remainder with the mixed-cusp cancellation and avoids the loose Taylor wording. That revision resolves the sole requested clarification; this audit passes the revised statement. The remainder of the argument is the previously fully reviewed proof, reported unchanged by its author. Later revisions require separate hash identification.

## 5. Scope and arithmetic cross-check

The case split covers every α∈(0,2): either p+α<2 or α>1 and the axis expansion applies. It allows all finite angular spectra and arbitrary signed compact C² corrections. It does not cover p≥1, nonstable generators, corrections retaining a compensating axis singularity, or kernels with a different singular angular profile.

The proof is analytic, not a grid certificate. The companion exact Fraction script `computation/kernel_reverse_design_cnd_check.py` checks only the scalar identity 4[1/(α-1)-1/α]=4/[α(α-1)] at α=5/4,3/2,7/4; it does not numerically evaluate an arbitrary angular spectrum or prove uniform error estimates. Those estimates are given above. No R9 matrix or LP was rerun, and no new package or MATLAB process was started.

Read dependencies: `proofs/kernel-interface-audit.md`, especially the normalization/C4 distribution requirement and the older fixed-coordinate axis asymptotic; `proofs/kernel-family-floor.md`, used only to preserve the frozen R9 scope. Their hashes and this file's hash are recorded in the actual companion output. No commit or push was made; root owns review and integration.
