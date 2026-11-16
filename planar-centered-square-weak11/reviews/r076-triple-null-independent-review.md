# R76 independent review: actual joint-null corner payments

2026-11-16, Asia/Shanghai. Strict threshold 1, centered axis-parallel
common-radius squares. New R76 files only. I read the complete proof
`proofs/r076-triple-null-corner-payment.md`, its arithmetic checker, input
and saved result. The analytic proof passes independent review; its final
SHA256 is `d0354c1085ab5ec702e9d1f799ba96fe7de79938fdbcfe40782e423a1ac20936`.
The author's six algebra fixtures and three activity fixtures were not rerun.
An initial analytic/source review was followed by one new directed n=5
fixture execution, recorded below. No old certificate was changed.

## Exact objects and the role of the full grid caps

The internal measure is the full uniform S_3, total 1, not just its local
origin of mass m=1/9. The same global outside measure satisfies every
closed square cap. Its strict effective SW window P=(-1,7/18)² contains
two positive outside atoms; outside P there may be further sources.
The receiver is Gamma=(-1/2,-1/9)². Boundaries do not affect area, but
strict source capture remains essential.

The full caps rule out every strict full or residual witness of radius
r>=1/2. Every remaining witness centered in Gamma lies strictly in P.
Thus complete all-radius full/residual localization is valid, and the
only old source participating locally is the origin. The local three-
source full halo has seven labels; its residual has both singletons and
their pair. No residual label is replaced by an independent measure.

For either outside point q in P, the bbox of {0,q} has span less than 1.
A closed unit square contains these sources, so the ORIGINAL full-grid
cap implies m+b<=1. Consequently u=sqrt(m+b)/2<=1/2. This is a separate
use of the full caps beyond localization. It does not require the local
old atom to be a maximum block.

## Clipped mixed-pair payment

Set e=m=1/9, v=sqrt(b)/2. The mixed pair receiver in Gamma lies in the
old-centered square [-u,-e]², of side L=u-e, if u>e. A coordinate overlap
with a q-centered interval loses at most u-v when its radius decreases
from u to v. If q is inside the interval of length L<=u, the radius-u
interval contains it and the radius-v overlap is at least min(v,L).
If q is outside it, the overlap is a clipped affine function with slope
at most 1. These cases include all endpoints.

For outer coordinate overlap lengths A,B<=L and inner lengths a,b,
the inequalities a>=(A-u+v)_+, b>=(B-u+v)_+ imply an area increment at
most L²-(L-u+v)_+². The actual all-radius residual singleton is exactly
the radius-v receiver square. Hence the mixed pair's true contribution
after even this partial residual subtraction is bounded by

    (u-e)_+²-(v-e)_+².

The function F(s)=(sqrt(s)-e)_+² is continuously differentiable at e²,
with derivative 0 below that point and 1-e/sqrt(s) above it. On
[0,1/4] its derivative lies in [0,1-2e]. Since u²-v²=m/4, the bound is
at most (1-2e)m/4=7/324. This uses the SAME b in u and v. It is valid
for arbitrary q; no symmetry or averaging is used.

## Triple-null payment and the strict remaining branch

The old singleton contributes only (1/6-1/9)²=1/324 in Gamma. If the
mixed triple receiver misses Gamma, the only possible new full labels
are this singleton and the two old/outside pairs. Subtracting each
pair's own residual singleton gives a safe upper bound, because the
actual complete residual halo contains both of them. Summing areas is
an upper bound even if the remaining pieces overlap. Thus

    |[H(sigma_3+nu)\H(nu)] intersect Gamma|
       <=1/324+2(7/324)=5/108<1/18,

with gap 1/108. This is a two-label payment with a null triple, not a
claim that the same old mass can pay arbitrary many mixed labels.

For source bbox [a,b]×[c,d] of {0,q1,q2}, independently minimizing the
two one-coordinate farthest-distance functions over [-1/2,-e] gives

    delta=max((b-a)/2,(d-c)/2,b+e,d+e,
              -1/2-a,-1/2-c,0).

The center feasibility interval [b-r,a+r] meets the closed receiver
interval precisely at the corresponding one-coordinate inequalities.
The two coordinates are independent, so this max is the exact minimum
farthest-source radius over the closed Gamma square. With
L=m+b1+b2, the open triple receiver meets open Gamma iff L>4delta².
When strict, a minimizing closed center can be moved into the open
receiver with margin. At equality no point has distance strictly below
sqrt(L)/2; a closed-rectangle drawing may only add zero-area boundary
contacts. Equality is correctly classified as null.

An active triple is eligible. Only THEN do the global caps imply L<=1:
if L>1 and its bbox span is less than sqrt(L), a containing bbox-side
square violates max(1,side²). Empty triples may have mass exceeding 1;
the proof does not impose the upper bound on them. Accordingly delta>=1/2
forces the null branch for every globally admissible mass vector. A
valid full-grid sum bound b1+b2<=U with m+U<=4delta² does the same for a
whole fixed-position mass domain. This filter retains the common global
mass vector, rather than separately maximizing different witness labels.

## A sharper fixed-position mass filter

For two DISTINCT fixed outside positions, each complete finite-source
bbox cap has outside coefficients (0,0), (1,0), (0,1), or (1,1).
The old-only rows hold automatically since the old total is 1. Taking
the smallest right side in the other three classes gives the exact
standalone full-grid feasible polytope

    b1,b2>=0, b1<=V1, b2<=V2, b1+b2<=U.

The bounds are nonnegative because zero outside mass is feasible. Each
individual class exists for distinct outside positions by testing the
singleton's zero-side bbox. The common class exists by testing both
positions. The maximum total mass is exactly

    U_eff=min(U,V1+V2).

For T=U_eff, choose b1=min(V1,T) and b2=T-b1; this is feasible and attains
T. Thus m+U_eff<=4delta² is a sharper uniform triple-null sufficient
condition. If both outside masses must be strictly positive, it remains
a valid upper bound; zero-cap cases simply remove the two-positive branch.
Coincident locations should be merged and use the single-position case.

For an original input with additional sources outside P, deleting those
sources while retaining the full old grid enlarges the feasible mass
domain. Hence this standalone polytope is a safe outer domain for the
original input, which is sufficient for an upper-payment filter. It is
not permission to paste independently feasible corner inputs.

## Evidence scope

The author saved six exact rational parametrizations and three strict/
tie geometric fixtures. They test constants and bbox activity, not their
global cap qualification. The code and saved result were read, not rerun.
No position search or halo integration was launched by this reviewer.
A finite exact-cap/floating-area exploration,
if reported elsewhere, cannot establish the remaining triple-active
payment. Arbitrary active triples, the general two-outside corner budget,
the unrestricted uniform-grid 7/2 target and exact C2 remain open.

## General n: the actual condition and the truncated-radius repair

I read the complete final core proof
`proofs/r076-general-grid-joint-null-corner-payment.md`, SHA256
`d3d9b8848fc29523017c5fee6a365d929532b8f3f0bb7734af9a489497b0b33b`.
Its Sections 1–8 pass independent review. For n>=3, set e=1/n²,
K=floor(n/2), and m_P=K²/n². The actual joint-null hypothesis excludes
every strict witness capturing BOTH positive outside positions and ANY
old source. This is stronger than merely testing the origin-plus-two-
outside label when P contains several old grid points.

A pure internal witness in the negative corner captures an anchored
k×l prefix. With j=max(k,l)>=2 it would need both
r>(j-1)h+e and r<j/(2n), impossible since
2n²[(j-1)h+e-j/(2n)]=n(j-2)+2j>0. Thus its true corner halo is only
the origin singleton, of area A0=(1/(2n)-e)². This term must not be
replaced by the singleton halo of the collapsed mass m_P.

For a witness capturing exactly one outside atom, collapsing all local
old mass to the origin safely enlarges full capture on the negative
receiver orthant. The essential radii are

    u_i=min(1/2,sqrt(m_P+b_i)/2), v_i=sqrt(b_i)/2.

The actual witness has r<1/2, and the outside singleton cap gives
b_i<=1. Therefore 0<=v_i<=u_i<=1/2 and u_i²-v_i²<=m_P/4.
One must NOT infer m_P+b_i<=1 from the original caps: the original
sigma|P and q_i need not fit together in a unit square. Retaining the
radius cutoff repairs this potential modeling error without imposing
caps on the auxiliary collapsed measure.

The clipped annulus estimate proved above then charges at most
(1-2e)m_P/4 to each of the two single-outside choices. Under actual
joint-null their union plus the true pure-old contribution gives

    |D intersect Gamma|<=g_n=A0+(1-2/n²)m_P/2<m_P/2,

with exact gap (n-1)/n⁴ for even n and (2n-3)/(4n⁴) for odd n. No
independent choice of mass vector occurs in this union bound.

For each of the K² ACTUAL old prefixes T_(k,l), form delta_(k,l)
using its entire bbox with q1,q2. Actual joint-null is equivalent,
under full global caps, to all inequalities

    kl/n²+b1+b2<=4delta_(k,l)².

The strict reverse inequality creates an eligible source label; the
global caps then force its mass<=1 and hence permit a witness radius
below 1/2. Any additional captured source mass only strengthens the
witness. Equality gives no strict receiver. Neither testing only
k=l nor replacing all prefixes by the origin label is justified.
The full-grid mass-polytope filter uses U_eff in every one of these
K² inequalities and is safe for an entire fixed-position mass domain.

For one common globally feasible input with at most two positive
positions per P, and actual joint-null in every two-position corner,
the four exact local restrictions give the complete deletion bound

    Delta<=7/2-4/n³+1/n⁴                    (even n),
    Delta<=7/2-1/n+1/(2n²)-2/n³             (odd n).

Both are strictly below 7/2; n=3 gives 85/27. This structural theorem
retains the same global measure and complete residual halo. Active
joint witnesses and corners with more than two positions remain
outside its assertion. The author's 127 parity rows, 8,255 pure-prefix
rows and 635 truncated-radius identities are finite scalar checks,
not the all-n proof; I did not repeat that kernel.

## A new admissible n=5 fixture distinguishes the quantifiers

One directed exact run verifies why the full-prefix menu is necessary.
Take S_5 and outside points

    q1=(-4/25,1/5), q2=(1/5,-4/25), each mass 17/200.

`computation/review_r076_prefix_joint_fixture.py` independently enumerates
all 342 nonempty bbox closures of the 27 raw atoms. All square caps
pass; there are 155 eligible closures, the unique maximum is the full
S_5 mask 33554431, and the eligible hypergraph is connected. The strict
nonterminal gap is 371/27500>0.

The origin-plus-two-outside label has mass 21/100, delta=6/25 and
activity gap -51/2500, so it has no strict corner receiver. But the
2×2 old prefix with the same two outside atoms has mass 33/100,
delta=7/25 and activity gap 41/2500>0. Its actual full witness at
x=(-41/1000,-41/1000), r=57/200 has excess 51/10000. All three residual
labels fail at this center. The center is already in the old singleton
halo; this is NOT a payment or area counterexample. It proves precisely
that origin-triple-null cannot substitute for actual joint-null at n=5,
even in the unique-maximum, connected, nonterminal full-cap domain.

Python 3.15.0 Fraction, PID 52859, exit 0 and confirmed absent, about
0.037 seconds wall time. The exact input is embedded in the saved
script; raw atoms, all counts, both activity gaps and witness checks are
saved in `prefix-joint-fixture.json` with an actual process receipt.
No position scan or full halo area was computed.

## Producer screen: source and scope review only

I read the finite screen code, grouping code, scope note, grouped result
and two execution receipts. The all-bbox cap scaling and three outside
coefficient classes are correct. The vertex/interior points are feasible
and all original cap rows are rechecked. With zero masses, extra virtual
source locations only create redundant subset receivers; caching by the
actual positive measure is valid. Formal triple activity is correctly
separated by positive outside arity in the saved grouping.

Producer records report one 1,024-position-pair batch, 8,143 mass points,
581,126 exact cap checks and 5,392 distinct actual-measure floating area
evaluations. Grouping reads saved data and adds no geometry. For actual
positive arity two, 935 formal triples are active and 4,136 null. The
largest displayed area, about .03812305, is not a certified area bound
or a proof of the arbitrary-position budget. PIDs 42374 and 46453 are
recorded exited. I did not independently replay the batch, re-count all
saved cases or regenerate statistics. Version/source binding does not
enlarge this review scope.

## Final report and version binding

I read final `rounds/076.md` in full, SHA256
`a4dc75089325d8deea3fbb7006f1c73e08c76e000d13d4fe9b20f5ea2dd069c0`.
The mathematical statement and all evidence ranges pass review. The n=3
and general-n proofs, the actual prefix counter-fixture, and the finite
screen are kept distinct. The general scalar check source, result summary
and receipt were read; their counts are not claimed as independent runs.
The single actual independent run is the new n=5 fixture above.

Frozen proof/component/report versions are pinned in this review's final
manifest. Root's unified integrator handles payload/dependency hashes,
staged bytes and relative links; no duplicate generic binding program
or old-round computation was added. No commit or push was made here.
