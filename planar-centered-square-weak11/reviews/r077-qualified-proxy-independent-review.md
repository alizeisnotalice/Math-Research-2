# R77 independent review: a qualified-input mixed-proxy obstruction

2026-11-16, Asia/Shanghai. Strict threshold 1, centered axis-parallel
common-radius squares. This review keeps three objects distinct: the
original full S20 input, its synthetic local three-source proxy, and
the uncapped mixed-only example. No actual S20 corner deletion area
or general C2 bound is inferred. All old-round files are unchanged.

## Proofs read and scope

I read both final proofs in full and pass their mathematics:

- `proofs/r077-mixed-only-obstruction-and-unequal-axis-caps.md`, SHA256
  `3c71513ed72567445d4900b96f819f49d59481a7a6e9a1d1b4c1550cb4a211aa`;
- `proofs/r077-qualified-grid-total-mass-proxy-obstruction.md`, SHA256
  `279bef687f0e378466558b2a0490b26998d33925b4413fc41248a8140677218e`.

“Mixed-only” means the union of old/outside mixed-label receivers,
without ADDING an old-singleton receiver as a separate union term.
Any overlap of a mixed receiver with that singleton is still included.
It does not mean subtracting the old halo. Both proofs state this
definition and retain the same complete outside residual throughout.

## The uncapped construction and the full-grid cap theorem

For old mass 1 and outside masses 4 at (-28/25,0) and 8 at (0,-161/50),
the second old/outside pair receiver is empty: its mass 9 is below its
span squared. The first mixed receiver and the triple receiver have
disjoint Y ranges. The first encounters only the mass-4 singleton
residual; the triple encounters only the outside-pair residual. The
other residual singleton is above the triple. The pair residual is
contained in the triple and must be subtracted. Direct areas give

    78/25-53sqrt(5)/50-161sqrt(13)/100+161sqrt(3)/50
       >25849/50000>1/2.

The radical orderings and rational lower bounds in the proof are
correct, including the empty pair and positive residual-pair height.
No negative-width rectangle is treated as an area. This is a valid
uncapped auxiliary counterexample only; its atom masses violate unit
caps. No independent uncapped halo execution was made here.

The companion full-grid theorem for q1=(-a1,0), q2=(0,-a2), a1,a2>0,
nonnegative masses, and n>=2 is also correct. Retaining all old grid
positions, all square caps are equivalent to the n+n+n² prefix rows
plus the individual atom caps. In a single-outside capture, the old
points form an anchored k×l prefix. Its second-direction span is below
1, so the tightest old mass for its fixed first index is the k×n
prefix. A double-outside capture gives a k×l anchored prefix and its
actual maximum axis span. The outside pair cannot be captured by a
closed square without the old origin. For 0<a1,a2<1 the individual
atom caps are implied by the first prefix rows; for arbitrary larger
distances they must be retained separately.

The bbox necessity argument does NOT require a smallest containing
square to capture exactly the stated bbox closure. Its mass is at
least the closure mass, so any additional atoms only strengthen the
necessary inequality. Conversely the source set captured by any
closed square is its tight-bbox closure, with span no larger than the
square side. This proves sufficiency without a hidden exact-capture
assumption or a receiver-grid restriction.

## One independent exact S20 execution

The new script `computation/review_r077_s20_proxy.py` reads both the
producer's compact corrected case and its complete raw 402-atom input.
It verifies the original 400 old positions, their 1/400 masses and
spacing 21/400, and the two outside positions/masses

    (-71/200,0), b1=1/4; (0,-171/200), b2=1/2.

All 45,302 raw nonempty bbox closures are independently reconstructed
with integer coordinate masks. Their outside-mask counts are 44,100
old-only, 401 of each single-outside type, and 400 double-outside.
Every cap passes; 7,386 closures are eligible. The unique maximum is
the complete S20 of mass 1. An eligible closed label joining the old
origin and both outside atoms connects them to the eligible full old
block. The exact nonterminal margin is 8319673/5280000>0.

The separately derived 440 axis-prefix rows also all pass and their
three minimum bounds match the raw closure classes exactly:

    V1=7/20, V2=5281/6400, U=5729/6400.

The minimum strict outside-row slack is 1/10. Outside singletons have
mass strictly below 1. Thus an outside-containing label of span<1
has mass<1; one of span>=1 has mass<span² and cannot be eligible.
For any nonclosed eligible subset, closing its bbox preserves span
and increases mass, so it cannot evade this argument or create a
missed maximum tie. Proper old-only subsets have mass<1, while the
full old block has mass 1 and span<1. This proves uniqueness as well
as reproducing the raw enumeration.

The actual 100 local old-prefix joint-activity predicates are checked
without integrating the full halo. Exactly nine are active: (k,1),
k=1,...,9. Their captured old masses are at most 9/400. This is only
a full-witness activity statement, not surviving deletion area.

## Exact synthetic proxy area and its limited implication

The same run constructs only the local SYNTHETIC three-source
measure with old mass m_P=1/4 at the origin. It uses the three mixed
full masks 3,5,7, complete residual masks 2,4,6, and Gamma cut 1/400.
Full radii are capped at 1/2. No original 400-old-source halo is built.
Area is computed independently by exact-radical horizontal y slabs
and merged x intervals, without importing producer code. This differs
from the producer's directed vertical-event and Cartesian algorithms.

In positive coordinates put u=sqrt(2)/4, v=sqrt(3)/4, a=71/200,
d=171/200. The first mixed receiver is [e,u]² and the triple is
[e,1/2]×[a,1/2], disjoint because u<a. The second mixed receiver is
exactly the residual-pair receiver on Gamma. The second residual
singleton misses Gamma; the first cuts only the first mixed receiver,
and the pair cuts only the triple. The recovered exact area is

    43sqrt(3)/200-101sqrt(2)/1600-24641/160000
      >206483/1600000>1/8.

The simple lower certificate uses sqrt(3)>433/250 and sqrt(2)<283/200.
It also exceeds the n20 single-corner sufficient threshold 80003/640000
by more than 12951/3200000. These are PROXY conclusions, not actual
corner deletion conclusions, and do not yield a four-corner or C2
counterexample.

The original failed n20 input is a distinct retained test. Its outside
pair mass is 5/8, less than span²=(83/100)², so that pair's receiver is
empty. It must not be represented by a negative rectangle area. I did
not independently rerun its area or the uncapped production case.

This single independent execution used Python 3.13.18/SymPy 1.14.0,
PID 85854, exit 0 and confirmed absent after wait; about 0.696 seconds
wall time. The result and actual receipt are
`computation/results/round077-independent/s20-qualified-proxy.json`
and its `-receipt.json`. Input pins are raw402
`d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a`
and compact cases
`9317309e36df5cafac9dd12535f643e2998c456f2fc64087b18767750c1cc6fb`.

## Positive-area separation and failed cap preservation

The proof's cell E=[-49/100,-47/100]², area 1/2500, is a valid analytic
separation. Every collapsed source is strictly captured at common
radius 99/200 with excess 199/10000. No residual label can receive E:
any outside capture needs r>47/100, whereas their total 3/4 is below
4(47/100)². Thus E belongs to the mixed proxy and is outside its old
singleton halo.

For the ORIGINAL input, r<=47/100 captures no source, including at
equality because the squares are open. If 47/100<r<1/2, each upper
endpoint is below 3/100<21/400; at most the old origin and the two
outside atoms are captured. Their mass 301/400 is below 2209/2500,
so reception fails. At r>=1/2 the original global caps exclude strict
reception. Hence the entire E is outside the actual full halo. This
is a proof for a whole closed cell, not a floating point observation;
no additional full-halo execution was needed or claimed.

The independent run also checks the displayed unit square before and
after ACTUALLY replacing all local old points by m_P at zero while
retaining old points outside P. Its original mass is 339/400, but
the replacement mass is 409/400>1. Therefore original full-grid
qualification is not cap preservation of a synthetic full-grid
replacement. The isolated three-source proxy has total mass 1, but
its selected old mass 1/4 is not the maximum eligible mass. These
different normalizations cannot be interchanged.

## Other executions and final report

The author scalar sources/results were read, not rerun. Their exact
root bounds, expansion and E inequalities supplement the analytic
proof; they are not independent area executions. Production's three
directed cases and qualification output were source/result reviewed
and pinned. Only the corrected S20 case was independently recomputed.

Root's finite floating probe was not repeated. Its saved scope and
receipt record one actual 1,536-case batch after a wrapper syntax
error that occurred before launching a child. The negative batch
does not prove an upper bound or retrospectively become the source
of the later exact counterexample. This review makes no independent
claim about all its sampled areas or count reconstruction.

I read final `rounds/077.md` in full, SHA256
`2d5dcae99a843317c44957d9e254f66d8cc9268a76003c54fbd4bdf2ee9150de`.
Its mathematics and evidence ranges pass review. The new result is a
qualified-input obstruction to a particular total-local-mass mixed
proxy. Actual grid deletion, the active-joint payment problem, the
general 7/2 target and exact C2 remain open. Root's unified integrator
handles full byte and link checking. No MATLAB, old-round rerun,
commit or push was made here.
