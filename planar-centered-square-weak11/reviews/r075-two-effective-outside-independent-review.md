# R75 independent review: two-outside obstruction and certified corner payments

2026-11-16, Asia/Shanghai. Strict threshold 1, centered axis-parallel
common-radius squares. R74 files are unchanged. One new independent
three-source execution was made, followed by a separately authorized complete
continuous-parameter certificate replay; no MATLAB, new search batch, commit
or push.

## Complete actual southwest deletion

Let sigma=delta_0 and nu=4delta_(-2,0)+4delta_(0,-2). I reconstructed the
full halo from all seven nonempty source subsets and the residual halo
from all three outside subsets. The mixed outside pair must be retained.
For any finite positive atomic measure, the complete halo is the union,
over nonempty subsets of mass m and bbox [a,b]×[c,d], of

    (b-sqrt(m)/2,a+sqrt(m)/2) × (d-sqrt(m)/2,c+sqrt(m)/2).

Every strict witness yields such a subset. Conversely an interior point
of a nonempty receiver rectangle admits a radius larger than all source
distances and smaller than sqrt(m)/2. This accounts for all radii, not a
chosen witness menu. Open receiver boundaries are finite line segments
and have zero area.

Use positive receiver coordinates X=-x,Y=-y, u=sqrt(5)/2, v=sqrt(2).
Modulo the singleton residual rectangles, the full halo consists of

    O=(0,1/2)²,
    E_x=(2-u,u)×(0,u), E_y=(0,u)×(2-u,u),
    T=(1/2,3/2)².

The entire residual union in this orthant consists of

    L_x=(1,3)×(0,1), L_y=(0,1)×(1,3),
    R=(2-v,v)².

O has no residual intersection and contributes 1/4. The part of E_x
outside T survives residual subtraction exactly on
(2-u,1)×(0,1/2), of area (u-1)/2; E_y gives the transpose. Inside T,
the two L rectangles have disjoint interiors and total area 1/2. R lies
inside T, has area 4(v-1)², and overlaps each L rectangle in area (v-1)².
Hence the residual union within T has area 1/2+2(v-1)² and the deletion
within T has area 4v-11/2. The four remaining pieces have disjoint
interiors. Therefore

    |[H(sigma+nu)\H(nu)]∩{x<0,y<0}|
      =1/4+(u-1)+(4v-11/2)
      =4sqrt(2)+sqrt(5)/2-25/4 >131/250>1/2.

The strict lower certificate uses sqrt(2)>707/500 and sqrt(5)/2>559/500,
whose squared errors are 151/250000 and 19/250000. Thus the two-effective-
outside extension of the uncapped m/2 orthant lemma is false. The R74
single-effective-outside theorem is unaffected.

If the residual pair R is omitted, the erroneous area becomes
sqrt(5)/2-1/4, an overcount of 6-4sqrt(2). This is why singleton-only
residual subtraction cannot certify the claim.

## One independent exact execution

`computation/review_r075_three_atom_payment.py` reads the raw three
positions and masses, reconstructs every subset receiver, and computes
the actual difference using horizontal y slabs and merged x intervals.
Coordinates and areas are exact SymPy radical expressions; no producer
code or saved rectangle rows are imported. This differs from the
producer's vertical-event and Cartesian-midpoint algorithms.

The result agrees symbolically with the preceding formula. The same
single run verifies the residual-pair omission error and two receiver
clips. Python 3.13.18, SymPy 1.14.0, PID 9219, exit 0, PID absent after
wait; approximately 1.57 seconds wall time. Result and receipt are
`computation/results/round075-independent/three-atom-payment.json` and
`three-atom-payment-receipt.json`. The raw input SHA256 is
`173a66e72d79c26b999d3fa622cf23de0ab7b7e61da094ef56993cfe56d7f3ad`.

The producer's directed input, code, result and complete analytic note
were read. Their seven-full/three-residual scope and hand partition pass
independent review. The author's two 72-bit geometry paths were not rerun.
No additional d, random input, or numerical optimization was attempted.

## The orthant is not the grid corner, and unit caps are not all caps

For the unscaled three-source fixture, deletion inside [-1/2,0]² is only
the old singleton square, of area 1/4. Inside [-1/2,-1/9]² its area is
49/324. These two exact clips belong to the unscaled fixture; neither
is a uniform S_3 input.

The fixture's actual maximum eligible mass is 9, not its designated
internal mass 1. The side-2 square [-2,0]² captures mass 9>4, and even
unit caps fail at a mass-4 atom. Consequently this is not a full-cap
uniform-grid counterexample or a violation of the actual 7/2 target.

Scale lengths by 1/3 and masses by 1/9. The local three-source measure
now has total and maximum eligible mass 1 and therefore satisfies all
normalized square caps. Its designated old atom has mass 1/9, not the
maximum mass. Its whole-orthant deletion is the original area divided
by 9, still greater than (1/9)/2. Thus imposing local caps with upper
normalization 1 alone cannot rescue the general whole-orthant m/2 lemma
when the designated old mass is smaller than 1.

The scaled local THREE-source crop on Gamma_3=(-1/2,-1/9)² is a different
area. In original positive coordinates it clips to [1/3,3/2]². The old
piece has area 1/36, the two stems contribute (u-1)/3, and the central
piece still contributes 4v-11/2. Dividing by 9 gives

    4sqrt(2)/9+sqrt(5)/54-209/324 <143/4536<1/18.

The first strict upper uses sqrt(2)<10/7 and sqrt(5)<9/4. This is an
analytic clip calculation, not an extra numerical run, and concerns
only the scaled local three atoms.

Adding the eight missing S_3 grid atoms gives a full 11-source embedding
whose outside points are (-2/3,0),(0,-2/3), each of mass 4/9. All closed
unit caps do hold: a unit square containing both outside points can
capture only the old origin, so its mass is at most 1. A unit square
containing exactly one outside point can capture at most one old column
or row, so its mass is at most 4/9+1/3=7/9. With no outside point its
mass is at most 1. These cases include closed boundary captures.

However, [-2/3,4/9]² has side 10/9 and captures four old atoms and both
outside atoms, mass 4/3>100/81 by 8/81. It gives an eligible label heavier
than S_3. The full embedding is therefore inadmissible for the grid
target, despite passing every unit cap. Atomic outside mass cannot use
the R74 density<=1 unit-to-all-square lemma. Nor may its complete corner
halo be identified with the local three-source crop by the r<1/2
localization theorem: the required full caps fail.

The qualification code/result's 225 unit-center strata were source/read
reviewed, not independently rerun. The preceding direct analytic cases
independently verify their all-unit conclusion and the nonunit defect.
I read the full final qualification proof
`proofs/r075-two-outside-obstruction-and-global-cap-qualification.md`,
SHA256 `cb14a03352f0d84c7e933d27e3ed02d0f6e539734ead22c06f4dd8ad40b88ff1`,
including its general d>=2 partition and analytic limit 5/4-1/sqrt(2).
The stems and triple inclusion-exclusion, the rationalizations for the
limit, and the scaled qualification all pass. No other d was evaluated.
The general multiple-outside corner budget, uniform-grid actual 7/2
target and exact C2 remain open.

## Full-grid capacity reduction: a valid replacement interface

I also read the complete core proof
`proofs/r075-orthogonal-two-outside-full-grid-caps.md`, frozen SHA256
`7a0e6c14b12f1ddb610e7c63f061400254f38863b3c314ee1528fad7e2553cf3`.
Its Sections 1–4 pass independent mathematical review. For n>=3, full
uniform S_n, and equal outside masses b at (-a,0),(0,-a), 0<a<=1,
the exact all-square condition is

    b <= B_n(a) = (1/2) min_(1<=k<=n)
                    [max(1,(a+(k-1)h)^2)-k^2/n^2].

Necessity uses the closed squares [-a,(k-1)h]². A square capturing both
outside atoms captures an anchored k×l old prefix; its maximal coordinate
count j controls its side and old mass by the j×j row. A square capturing
only the left atom also captures an anchored prefix. If k is its old
column count, the double-outside k-row dominates its single-outside mass
because 2kl-k²<=n²<=n² max(1,(a+(k-1)h)²). Empty-old and no-outside cases
are covered separately. The restriction a<=1 is used in the k=1 bound;
the theorem is not asserted beyond it. All closed source ties are retained.

For e<a<=1, put K=1+floor((1-a)/h). Then 1<=K<n. The unit-side rows
decrease with k, while adjacent nonunit rows have difference

    2h R_k+h²-(2k+1)/n² >3/n²+h²>0.

Only rows K and K+1 can bind. The unit equality belongs to K, and a<=e
permits only b=0. I independently checked the n=3 five-piece boundary,
its two radical switch points, and B_3(2/3)=32/81. These are analytic
checks, not additional computation.

For admissible n=3 parameters, the four old-containing full receiver
rectangles minus the three complete residual rectangles, all clipped to
Gamma, give the exact area. b=0 and a=1 cause no endpoint gap. The
parameter-box certificate interface is also sound: enlarge every full
rectangle with upper radii/lower a, and subtract a common residual subset
formed with lower radii and the safe two a endpoints. Monotonicity of B_3
justifies pruning b_L>B_3(a_U) and reducing b_U to B_3(a_U). This validates
an upper-certificate mechanism, not an already completed area budget.
The completed parameter-cover execution is reviewed separately below; no
area conclusion is inferred from sampled parameters.

## Complete independent parameter-cover replay

Root explicitly authorized this distinct new verification after the
one-input counterexample audit. I did not regenerate the tree, choose new
parameters or rerun the producer's 961-point probe. The saved compressed
certificate SHA256 is
`36bd0b4bad73c7047601a66f96c6cf76a3ac77e1794ecdcee7c1b4e819cfa99c`.

`computation/review_r075_orthogonal_tree.py` independently traverses the
entire tree. All 71 nodes are reached once; the 35 binary splits exactly
cover their retained parent boxes, including common closed boundaries.
All seven removed upper-b strips are infeasible by the independently
recomputed monotone B_3(a_U). Each node's kept box is the exact cap-trimmed
original box. No pending, unresolved, orphaned or cyclic node exists.

For every retained node, I checked all four dyadic-root records by squaring
integers against the exact rational mass, including whether the upper
integer equals the lower integer at a perfect square. This checks 284 root
records. I reconstructed all four full outer and three residual common-inner
rectangles from box endpoints and these root bounds, clipped them to the
actual Gamma, and compared each saved rectangle exactly. A different
horizontal y-slab algorithm then recomputed all 71 area upper bounds,
using 572 slabs in total. It imports no producer code.

All 36 terminal leaves are certified below 1/18. Their worst upper is

    1563026333441017474428935886295754233898431 /
    28224380641890319913506143438820582531006464,

with positive target margin

    4994813330666965210294304749833684490817 /
    28224380641890319913506143438820582531006464.

This is a proof for the entire closed feasible parameter domain
1/9<=a<=1, 0<=b<=B_3(a), not a finite list of tested parameters. The
remaining 0<a<1/9 branch forces b=0 and has pure internal Gamma area
1/324. Its analytic treatment completes the family. The replay used
Python 3.15.0 Fraction, PID 26310, exit 0, PID confirmed absent, about
0.038 seconds wall time. Results and receipt are
`orthogonal-tree-replay.json` and `orthogonal-tree-replay-receipt.json`.

The other three effective corner windows contain no outside source. Under
the global caps their actual corner deletion is the pure old singleton
crop, each of area (1/6-1/9)²=1/324. The SW bound 1/18 thus gives total
corner deletion at most 7/108. The noncorner container has area 80/27, so
the **complete actual deletion** in this equal-mass orthogonal full-S_3
family obeys

    Delta <=80/27+7/108=109/36<7/2.

This conclusion does not require uniqueness, connectedness or a
nonterminal hypothesis. It is restricted to this family's equal outside
masses and positions; it does not settle arbitrary two outside positions,
unequal masses, or the unrestricted uniform-grid question.

I also read the full certificate/structural-class proof
`proofs/r075-certified-orthogonal-corner-payment.md`, SHA256
`6669be9c13302bbba84611999d44f177b7373018441aad2aad82c4c1ce55af67`.
It passes. It correctly combines the parameter certificate with R74
localization for a larger class: one already globally feasible full S_3
input may have in each effective P either an arbitrary single positive
outside position, or the corresponding D4 image of an equal orthogonal
pair, with different pair parameters in different corners. Additional
outside sources away from all P are arbitrary. Restricting nu to a P
while retaining full S_3 preserves the global caps and forces that pair's
parameters into the certified domain. Every corner costs at most 1/18,
so the same actual input has complete deletion

    Delta <=80/27+4/18=86/27<7/2.

This does not paste independently feasible pair inputs; it starts with
one common feasible nu. For pure pairs at a=1 the source-boundary rule
applies. At a>1, if full global capacities are separately retained, all
four P are empty of outside sources and the stated coarse full-area
bound also follows. No a>1 capacity formula or unequal-area conclusion
is imported into these claims.

## Unequal-mass capacity addendum

I read the complete short analytic addendum
`proofs/r075-unequal-orthogonal-outside-capacity-polytope.md`, SHA256
`6284834949379d97fbe2b142b2f914af8df8d749badf694814d83cf0069bdb17`.
It passes independent review without another execution. With
C_k=max(1,(a+(k-1)h)²), its exact all-n polytope is

    b1+b2<=U=min_k(C_k-k²/n²),
    b1<=V=min_k(C_k-k/n), b2<=V.

The individual necessary rows follow by containing the outside point
and the k×n old prefix; additional captured sources only strengthen
necessity. Sufficiency covers the same anchored captures as the equal
proof. Zero masses and closed ties are harmless. V also needs only the
last unit and first nonunit prefix, since the nonunit row difference is
strictly greater than (n+2)/n²+h². Taking the row minimizing V in
C_k-k²/n²<=2(C_k-k/n) proves U<=2V, so the equal slice reduces to b<=U/2.
This is a capacity theorem only. Neither the independent tree replay nor
the 109/36 area conclusion extends to unequal masses.

## Root symbolic evidence, read without re-execution

I read root's Wolfram source, final result, both receipts and execution
scope note. There were two attempts, with exactly one successful
mathematical run. The initial syntax failure returned launcher exit 0
but produced no result and no recorded kernel PID; its null result hash
is not a successful check. The final repaired source produced all five
piece-counterexample existence values False, the ordered-breakpoint
test True, and four joining identities True. Launcher PID 15038 and
kernel PID 15058 were confirmed absent after exit 0.

This final Wolfram 15.0.0 run checks only the n=3 scalar capacity
envelope on 0<a<1. Endpoint a=1 is covered analytically by the core
formula. It is not an all-n proof or a halo-area computation. Root's
component manifest pin is
`740beffb2ec3481b4110e13ec1163d6bf81377b6f19739af9e94d1676478f03e`.
I launched no Wolfram or MATLAB and do not count these attempts among
the two actual independent executions recorded above.

## Final report and archive scope

I read the full final `rounds/075.md`, SHA256
`f87ef11a5fa6dd5936034cb2a51aae71cfa20c0211afb562ee711c9d6aa25c9e`.
It passes mathematical and evidence-scope review: the uncapped
counterexample, the all-n capacity results, the n=3 equal-mass certified
areas and common-input structural corollary are kept distinct. The
producer midpoint/225-unit runs are not attributed to independent replay,
and the failed initial Wolfram attempt is not counted as a success.

The producer's final orthogonal certificate note was also read and passes;
its maximum sampled area is correctly separated from its uniform Phi.
The producer searched once for a tree; the independent audit traversed
that saved tree once. These are different executions. New independent
scripts, results and actual receipts, together with frozen source pins,
are listed in this review's component manifest. Root's unified integrator
handles full payload hashes, staged bytes and relative links; no generic
source-binding program or old-round execution was added here.
