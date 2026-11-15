# R62 V3: telescope interface density jumps before source-curvature bounds

2026-11-15. New independent V3 files only; frozen V1/V2 run sources are
unchanged. The objects, general positive occupation kernel, whole-cell
bilinear interpolation lemma and exact halfmass formula are those of
`r062-coupled-wholecell-source-interface.md`. This increment changes the
curvature representation, without changing the kernel family or D sampling.

## 1. Exact continuous derivative identity

Let m=2N+1 and partition [−(N+1/2)delta,(N+1/2)delta]^2 into square physical
cells C_uv centered at (u*delta,v*delta), −N≤u,v≤N. Let their densities a_uv
be nonnegative, D4-symmetric, and put all exterior densities equal to zero.
Write T_u(x)=integral_[u*delta−delta/2,u*delta+delta/2] k_h(x−z)dz.
Then

    (P_h f)(x,y)=sum_(u,v) a_uv T_u(x)T_v(y),
    T''_u(x)=k'_h(x−(u−1/2)delta)
              −k'_h(x−(u+1/2)delta).

For every vertical interface e=(k−1/2)delta and source row v, define

    d^x_(k,v)=a_(k,v)−a_(k−1,v),  −N≤k≤N+1.

The first term is the density on the right, the second on the left. The
convolution is smooth for h>0 and the sums are finite, so collecting the
shared endpoint terms is exact:

    (P_h f)_11(x,y)
      =sum_(k,v) d^x_(k,v) k'_h(x−e) T_v(y).         (1)

For horizontal interfaces use d^y_(u,k)=a_(u,k)−a_(u,k−1), top minus bottom,
and obtain

    (P_h f)_22(x,y)
      =sum_(u,k) d^y_(u,k) T_u(x) k'_h(y−e).        (2)

The boundary terms in (1)–(2) are essential and are included by the exterior
zero convention. Thus each direction has exactly m(m+1) interfaces. A
constant-density block has zero jumps at all its internal interfaces; only
the two exterior interface lines in each direction survive. This is exact
continuous cancellation, not a finite-grid approximation of a derivative.

## 2. Signed jumps have a linear upper-bound interface

Fix a target subrectangle B. Let Psi=[L,U] enclose an interface function,
either k'_h(x−e)T_v(y) or T_u(x)k'_h(y−e), throughout B. Introduce d+,d−≥0
with d+−d−=d. For every actual psi in [L,U],

    d*psi ≤ d+*U−d−*L.                              (3)

Indeed the right side minus the left is

    d+*(U−psi)+d−*(psi−L) ≥ 0.

It is not necessary to impose complementarity d+*d−=0. If both are positive,
(3) remains valid, though it may be weaker. This keeps the proposal problem
linear.

Let ox,oy be signed upper bounds for partial_11 p_h, partial_22 p_h on B.
For g=lambda*p_h+P_h f, use nonnegative auxiliaries eta_x,eta_y satisfying

    eta_x ≥ lambda*ox + sum_vertical(d+*U−d−*L),
    eta_y ≥ lambda*oy + sum_horizontal(d+*U−d−*L).   (4)

Then g_11≤eta_x and g_22≤eta_y. The same four coupled vertex rows as before,

    g(vertex)−a_target ≥ (ell_x² eta_x+ell_y² eta_y)/8,

certify the whole target subrectangle by the main R62 interpolation lemma.
The complete vertex value still uses the original positive cell
probabilities. Only the second-derivative calculation is telescoped.

The one-dimensional k' ranges use the exact extrema from V2. Probability
ranges use their maximum and minimum absolute distances from the source
center. For multiplication with k', which can be signed, all four endpoint
products are tested. Taking only a product of positive lower endpoints
would be wrong here.

## 3. Code-to-formula mapping

`computation/cauchy_occupation_jump_aux_lp.m` is a standalone new MATLAB
proposal. Its physical orbit map gives each grid cell the same a_j as its
D4 orbit. For MATLAB cell index k=1,...,m, centers are
(k−N−1)delta; for interface index k=1,...,m+1, edges are
(k−N−3/2)delta. Each vertical row of JD adds the right cell's orbit coefficient
and subtracts the left cell's, where present. The second block does top
minus bottom at each horizontal interface. Coefficients are densities;
there is no orbit-size factor in a jump.

The variable order is

    (lambda,a), dplus_(all vertical then horizontal),
    dminus_(same order), eta_x, eta_y.

The equality matrix is [JD,−I,+I,0], so dplus−dminus=JD*(lambda,a) has the
sign in (1)–(2). The curvature rows implement (4) as

    lambda*o + dplus*U − dminus*L − eta ≤ 0.

There are J=(N+1)(N+2)/2 source-density orbits,
B=J*subdiv² target subrectangles, E=2m(m+1) physical interfaces,
J+1+2E+2B variables, 6B source/curvature inequalities, and E jump equalities.
For N16,subdiv2 this is E2244,B612 and5866 variables. The adjacency equality
matrix is sparse; edge interval evaluation is vectorized. Ordinary doubles
are proposals only. All previously recorded frozen MATLAB/Python hashes
are unchanged.

The new exact checker `computation/r062_coupled_source_verify_v3.py` parses
only lambda and the original source densities. It reconstructs all physical
jumps with exact Fractions and omits exactly zero jumps. It does not trust
the LP's jump parts or curvature auxiliaries. It bounds (1)–(2) by signed
160-bit outward interval arithmetic and evaluates original coupled vertex
source values independently. Its coverage, negative-witness and repair
semantics are unchanged from V1/V2. The tiny self-test additionally asserts
that all internal constant-density jumps vanish exactly.

## 4. Targeted comparison and evidence boundary

Use exactly h=.15,delta=.125,N16,subdiv2 for the first V3 test, preserving
V2's kernel family, objective and D quadrature/sampling. Comparing the
actual costs isolates the effect of operator-level interior cancellation.
Root owns MATLAB launch, time bound, exact checker execution and prompt
process exit. No additional V3 grid was launched by this author.

Prior actual root results are source-certified but uncompetitive:

| Version | h,delta,N,subdiv | Floating sampled-D cost | Exact source | Full D |
|---|---|---:|---|---|
| V1 separate positive curvature | .15,.125,16,2 | 6.97777284372976 | PASS612 leaves, no repair | not certified |
| V2 signed curvature auxiliaries | .15,.125,16,2 | 6.51716115659593 | PASS612 leaves, no repair | not certified |

V2's exact halfmass/pi is32413700406210501/15625000000000000. Root receipt
records MATLAB PID51819 exited0, absent after completion,22.2628s wall;
its exact source checker took1.81964s. Candidate, source certificate and
receipt are copied byte-for-byte to
`computation/results/r062-coupled-source/v2-h015-d125-s2/`.

Root subsequently ran the V3 target comparison. The same-parameter
floating cost is6.326423184817601; exact source PASS covers all612leaves,
no repair,2.7094s. Root receipt records MATLAB PID55789 exit0, absent after
completion,26.3006s wall. Exact halfmass/pi is
322202087025565676571/160000000000000000000. Candidate, certificate and
receipt are copied byte-for-byte to
`computation/results/r062-coupled-source/v3-h015-d125-s2/`.
Thus operator-level cancellation improves this test from6.51716 to6.32642,
with source positivity fully certified, while the cost remains uncompetitive.
This does not establish a new general C2 bound or a family optimum.

Root also ran two isolated V3 subdiv4 comparisons, without changing source:

| h,delta,N,subdiv | Floating sampled-D cost | Root exit/absence | Exact source |
|---|---:|---|---|
| .15,.125,16,4 |6.180910217311804| PID56582 exit0/absent | PASS2448 leaves, no repair |
| .075,.125,16,4 |11.130723502016476| PID56583 exit0/absent | PASS2448 leaves, no repair |

Their exact source checks took10.3157s and20.9577s respectively. All five
actual candidates are now source-certified on6732complete terminal leaves,
with no unresolved or unchecked rectangle and no added origin repair.
Their actual saved outputs, source certificates and root receipts are under the corresponding
`computation/results/r062-coupled-source/v3-h015-d125-s4/` and
`v3-h0075-d125-s4/` folders. No new job or further grid was started here.
Full-D certification remains absent for every candidate.

The h comparison reveals a resolution issue, not a lower bound on the
occupation kernel family: f is constant on a complete source cell, while
g is continuous and must dominate both one-sided cell densities at shared
interfaces. Decreasing h at fixed delta increases delta/h; the central-cell
corner value of p_h becomes much smaller than its center peak. The R60
approximation theorem refines cells for fixed h before letting h decrease;
it does not promise improvement when h decreases on a fixed coarse cell
grid. This mechanism and the present observations motivate eventual
co-refinement or central adaptive cells, without claiming an executed
improvement or a general obstruction. No V4 is implemented.
