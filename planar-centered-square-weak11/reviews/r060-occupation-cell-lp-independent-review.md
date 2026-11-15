# R60 occupation-cell LP: independent geometry, source and result audit

2026-11-15. Read-only audit of frozen MATLAB `cauchy_occupation_cell_lp.m`, SHA256 `dab5992d41449cfeb047455d472d106ee05ce691e82f8f41f7e65bc241c9c1fd` (4897 bytes), and source proof `r060-positive-cell-occupation-interface-independent.md`, SHA256 `2674fedc81adc0ff781f9274438579e09a44d8a05265a1d10ce26b9ddf488e52` (10253 bytes). No source changes or MATLAB/LP launch were made by this audit. Earlier R55 and R60 rectangle/Jensen files are preserved.

**Outcome:** matrix direction, D4 geometry and exact half-mass formula are correct. Twelve completed outputs reconstruct. The lowest-cost midpoint candidate is rigorously C4-invalid on an explicitly certified positive-area set. Cellinf screens remain numerical candidates with costs above 4.7; no new general C₂ upper bound follows.

## 1. Geometry, source coefficients and matrix direction

For δ>0 and integers 0≤j≤i≤N the source representatives are δ(i,j), completed under both coordinate reflections and interchange, with duplicate points removed. These disjoint orbits contain each grid center δ(k,l), −N≤k,l≤N, exactly once. The physical square cells of side δ tile [−(N+1/2)δ,(N+1/2)δ]² up to zero-area boundaries. Thus J=(N+1)(N+2)/2 and the physical cell count is (2N+1)².

The vector a_j is a density per physical cell, not total orbit mass. If f=Σ_j a_jΣ_(z∈orbit_j)1_(cell_z), its mass is δ²Σ_j size_j a_j. A transition entry L_ij sums over source cells in orbit j and evaluates at target representative i. There is no division by orbit size in L, D, or f. Orbit size belongs in the mass objective. The implemented orientation `L(:,j)` is correct.

For the one-dimensional Cauchy density, the probability of a centered interval of width δ at distance d≥0 is

B_h(d)=(1/π)atan2(hδ,h²+d²−δ²/4).

It is decreasing in d since its derivative is

−2dhδ/[π((h²+d²−δ²/4)²+(hδ)²)]≤0.

Therefore for a target cell centered at representative r_i, the code's `cellinf` values

b_i=p_h(|r_i,1|+δ/2,|r_i,2|+δ/2),
L_ij=Σ_(z∈orbit_j) B_h(|r_i,1−z₁|+δ/2) B_h(|r_i,2−z₂|+δ/2)

are valid entrywise lower bounds over the entire target cell. Different summands need not attain their minima at a common corner: summing their individual lower bounds is conservative and remains sufficient. For first-wedge representatives r_i≥0 the b formula in the code is exactly the displayed one.

The midpoint version replaces those lower bounds by center values. This only imposes source positivity at centers. It does not certify the whole target cell; the actual completed outputs show the gap is material.

As an orientation check independent of the code, aggregated matrices satisfy size_i L_ij=size_j L_ji. This detailed balance follows from symmetry of each physical-cell pair contribution and the D4 action. High-precision reconstructions passed this identity for both midpoint and cellinf matrices. Midpoint row sums also equal the probability that a p_h increment from its representative lands in the entire square grid domain, verifying that source-orbit counts are neither omitted nor averaged.

## 2. Kernel, D constraint and exact objective

Write λ for `v(1)`, μ=λδ₀+f(x)dx, and Q_h=U−P_hU=2π∫₀^h p_s ds. The implemented kernel is

K=Q_h*μ=λQ_h+Σ_j a_jΣ_(z∈orbit_j)(Q_h*1_(cell_z))≥0.

For each physical cell, convolution equals

2π∫₀^h B_s(|x₁−z₁|)B_s(|x₂−z₂|)ds.

Thus `kmatrix` has the correct factor 2πh multiplying [0,1] Gauss weights; it does not need an additional area factor because B_s already integrates over the source interval. The exact full-plane half integral is

(1/2)∫K=πh[λ+δ²Σ_j size_j a_j].

This is precisely `cost`, including the free λ. No finite-tail integral is needed for this objective. The constraint block [−b,I−L] imposes a≤λb+La; [−q,−D] imposes sampled K≥1. Fine-grid rescaling by `dm` preserves homogeneous source constraints mathematically, but its D normalization remains sampled and nondirected.

`max(v,0)` after the floating LP changes the returned vector. Therefore LP success alone does not imply source feasibility of the retained vector. The code appropriately reports the reconstructed `lowerres`, and all final vectors require exact nonnegativity, directed source inequalities and full-D verification before becoming certified kernels. Tiny negative printed `lowerres` values cannot be silently promoted to exact feasibility.

## 3. Source proof: finite C4 and the density claim

The source proof's finite C4 argument is valid under its stated hypotheses. Nonnegative lower bounds b,L and a≤λb+La imply ν=λp_h+P_hf−f≥0 everywhere except irrelevant cell boundaries (and hence as a density measure). Outside the finite cell domain f=0, so positivity is automatic. Its total mass is λ. The identity

A K=2π(ν−λδ₀), A=−|D₁|−|D₂|,

follows from A Q_h=2π(p_h−δ₀) and finite-measure convolution. Thus K itself, without a positive-part repair, has distributional C4 away from zero. Even f supplies the required symmetry, K is nonnegative/L¹, and the exact integral above follows from Tonelli. For λ=0 the strictly substochastic finite matrix forces a=0, so a nontrivial feasible D-normalized kernel has λ>0.

The density proof separates three limits correctly. For fixed h, the discrete-observation killed process has μ_h=ΣT_h^nδ₀ and H_h=Q_h*μ_h. Its occupation continues within each h interval until the next observed exit; it is not prematurely identified with the continuously killed process. Geometric survival bounds from jumps larger than 2B make μ_h finite and provide a terminal surviving-mass settlement for every finite observation horizon.

Along nested dyadic observation times, τ_h decreases to the true exit time if the exit lies in the interior of Ωᶜ. The proof handles rectangle/polygon boundary hits using point polarity of one-dimensional Cauchy projections: its resolvent argument has finite potential away from a point and divergent potential at the point, so strong Markov forces zero hitting probability. Right continuity then identifies the observed exits, with an integrable survival envelope. Actual pathwise occupation domination H_h≥H yields the exact L¹ difference 2π(Eτ_h−Eτ), preventing mass escape.

For D compactly inside Ω, exit sources remain separated from D. The common U singularity cancels in H_h−H, and uniform continuity of shifted U on those separated source locations gives uniform D convergence. At fixed h, rational cell lower approximations, Neumann partial sums and the displayed substochastic error bounds give convergence to H_h when the cell mesh and π bounds are refined. The order “first choose h, then refine cells/π/Neumann tail sufficiently” is essential. This supports the existence/density statement. It provides no evaluated continuous-exit error for any of the twelve coarse screens, and no numerical kernel-cost improvement by itself.

The low cellinf costs below are not derived from the source proof's rational Neumann example; they are separate floating LP proposals. That distinction is maintained.

## 4. Actual twelve-case reconstruction

Read execution receipt SHA256 `ac5fb26dcbb4154c53dae11240873b0c582a2a036237376ccaabdd3c1ca936a9`. Every receipt-listed final result hash was checked before extracting a minimal reproducibility input. All twelve final outputs have successful runs; two earlier service-error attempts exited 1 and were retried separately. The receipt reports all fourteen owned process IDs exited; this audit did not independently inspect processes or start/stop any process.

The table lists printed exact-mass-formula costs after the code's sampled D rescaling. They are floating screening values.

| h | δ | midpoint cost | midpoint sampled source minimum | cellinf cost |
|---|---|---|---|---|
| .15 | .125 | 3.984887185270 | −2.8518614630 | 8.432387736387 |
| .15 | .0625 | 3.968208315983 | −1.6683443245 | 5.993784747860 |
| .3 | .125 | 4.306025472745 | −.3845707542 | 5.318204595696 |
| .3 | .0625 | 4.297698579323 | −.2234367067 | 4.718867226623 |
| .6 | .125 | 4.997518287098 | −.05256448625 | 5.180836669930 |
| .6 | .0625 | 4.989873551071 | −.02798159862 | 5.071068241608 |

N=16 for δ=.125 and N=32 for δ=.0625. Both use square halfwidth near 2. For every result, interpreting JSON numbers as exact binary64 rationals, this audit independently expanded its D4 orbits, reconstructed μ mass and πhμ mass, and matched the printed cost within 10⁻¹². It also recomputed the posted source-minimum point by summing the actual λp_h and all source-cell Cauchy probabilities; the twelve reconstructed residuals agree within 10⁻¹⁰. The six midpoint negative residuals are numerical witnesses; the strongest retained candidate has an additional exact certificate below.

All six cellinf outputs have floating conservative constraint minima slightly below zero (roughly 10⁻¹¹ or smaller), while their center/corner source samples are positive. Neither observation is a directed source certificate. Their lowest screening cost is 4.71886722662315 at h=.3,δ=.0625. These twelve cases produced no useful candidate below the inherited 3.590865312933 comparison target; this is not a lower bound for the entire occupation family.

## 5. Rigorous invalidity of the best midpoint candidate

For the h=.15, δ=.0625 midpoint candidate, source residual was printed at target orbit 5, center (1/8,1/16), corner (5/32,3/32). A corner is a density-interface point, so this audit moved it inward by the exact fraction ε=10⁻⁸ of the vector toward that cell center. The resulting point is strictly inside that physical cell, with positive distance to every side.

The retained binary64 λ and all density coefficients were converted to exact Fractions and checked nonnegative. An exact interval computation gives

ν(x)=λp_h(x)+(P_hf)(x)−a₅ < −1.668344309654557.

The displayed decimal summarizes an outward rational upper endpoint saved exactly in the JSON certificate. All interval operations are integer arithmetic on a 2¹⁹² lattice. Arctangent bounds use 80 alternating terms after reciprocal/π/4 reduction to an argument at most 1/2. Each multiplication/division is rounded outward, and the next-term remainder is included. π is enclosed by Machin's identity; its tangent identity was checked exactly. Thus this negativity does not rest on high-precision floating agreement.

The source smooth part has each partial derivative bounded by 2 μ(R²)/(π_lower²h³). Inside the target cell f is constant. Combining this bound with the exact distance to the cell sides proves ν<0 throughout an explicit L¹-radius ball of radius 1.5625e−10 around the certified point. This ball has positive area and excludes the origin. Therefore A K has negative density on an open set and the candidate strictly fails C4. This establishes invalidity of this specific midpoint kernel, not impossibility of the occupation construction or midpoint-generated proposals after repair.

## 6. Actual bounded tools, limitations and files

Final script `computation/r060_occupation_cell_lp_independent.py` ran successfully with Python3.13.18, mpmath1.3.0 at 80 decimal digits and NumPy2.5.4, one process. Deterministic exact geometry tests used N=1,2,4. High-precision matrix tests verified detailed balance, whole-square row sums and sampled entrywise cellinf domination. Five independent GL24 cell tests found relative errors up to +2.585e−6 just outside a source corner and −1.347e−7 at its center. These GL tests are nondirected numerical evidence; the exact negative-source certificate is separately implemented and labeled.

No full-D minimum was certified. The source positivity of cellinf proposals remains unverified with directed arithmetic. Exact mass identities and the negative-source certificate do not imply a general C₂ upper bound.

The script prefers integrated topic-relative MATLAB/proof paths and uses historical clone paths only if those files are absent. Its local minimal input archives all twelve retained coefficient vectors, orbit representatives and result-source hashes, so the rigorous best-midpoint certificate does not depend on disposable output directories. Reproduction from the topic directory:

    /Users/zzh/.venvs/main/bin/python3 computation/r060_occupation_cell_lp_independent.py

The adjacent delivery manifest records the final review/script/input/output paths, hashes and sizes. No commit or push was performed.
