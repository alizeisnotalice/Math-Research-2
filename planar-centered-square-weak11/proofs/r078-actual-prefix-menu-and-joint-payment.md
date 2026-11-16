# R78: complete actual prefix menus and a nine-prefix joint payment

2026-11-16. Strict threshold 1, centered axis-parallel squares with a common
radius. All source positions and masses below belong to one actual input.
No prefix is replaced by its total local mass, and no outside measure is
chosen separately for different receiver labels.

## 1. The finite actual-menu theorem

Let n>=2 and put

    e=1/n², h=(n+1)/n², w=1-e, K=floor(n/2),
    sigma_n=(1/n²) sum_(i,j=0)^(n-1) delta_(ih,jh),
    Gamma=(-1/2,-e)², P=(-1,1/2-e)².

Let nu be a finite positive atomic outside measure such that mu=sigma_n+nu
satisfies every closed square cap

    mu(Q)<=max{1,side(Q)²}.                             (1)

The four closed unit squares with coordinate intervals chosen from [-e,w]
and [0,1] each contain all internal mass 1. Hence nu vanishes on their union
[-e,1]², including its boundary; positive outside locations cannot coincide
with old grid sources.

Assume nu restricted to P consists of at most two positive locations q1,q2,
with fixed masses b1,b2; all other outside sources are outside P. For notation
suppose both are present and distinct; absent or merged locations only
decrease the menu size. Let I={1,2}. For 1<=k,l<=K define the actual old prefix

    S_kl={(ih,jh):0<=i<k,0<=j<l}.

For each B subset I, form T_klB=S_kl union {q_i:i in B}, with mass

    M_klB=kl/n²+sum_(i in B) b_i.                       (2)

For a nonempty pure-outside subset B, use T_B={q_i:i in B} with mass M_B.
For any such label T, let its actual source bbox be
[a_x,b_x]x[a_y,b_y], and define

    r_T=min{1/2,sqrt(M_T)/2},
    E_T=(b_x-r_T,a_x+r_T) x (b_y-r_T,a_y+r_T).           (3)

An interval with upper endpoint no greater than its lower endpoint is empty.
The same scalar r_T is used in both coordinates. Let

    U=union_(1<=k,l<=K, B subset I) E_(T_klB),
    R=union_(nonempty B subset I) E_(T_B).

Then the following are exact set identities on the open corner:

    H(mu) intersect Gamma = (U union R) intersect Gamma,
    H(nu) intersect Gamma = R intersect Gamma,
    Delta intersect Gamma = (U\R) intersect Gamma.    (4)

Thus at most 4K² old-containing prefix/subset labels plus three pure-outside
labels represent the complete actual full corner halo. The residual needs
all three pure-outside labels, including the pair. Pure old labels B=empty
are retained. Empty or coincident receiver rectangles may be removed, but
not inferred to be empty from midpoint tests.

## 2. Proof of completeness, including strict endpoints

First, every strict full or residual witness has radius r<1/2. For r>=1/2,
(1) bounds full mass by 4r²; residual mass is no larger. For x in Gamma and
r<1/2, the whole centered square is strictly inside P: each coordinate lies
between -1 and 1/2-e. Consequently sources outside P play no part in any
such witness, for full and residual halos alike.

At a negative receiver coordinate x_j, the lower endpoint x_j-r is negative.
Captured nonnegative grid coordinates therefore form an initial prefix,
with entry condition ih<x_j+r for open squares. The two axes share the same
r, and their Cartesian product is precisely S_kl if the captured old set is
nonempty. If either coordinate count is zero, there are no old sources and
only a pure-outside label is needed. Since x_j+r<1/2-e, each nonzero count
is at most K. The inequality i(n+1)<(n²-2)/2 has exactly the integer solutions
i=0,...,K-1, which proves this count without an asymptotic approximation.

Now take an actual strict open full witness at x. Its captured set consists
of one of those exact old prefixes and one outside subset B, or a nonempty
pure-outside subset. Its exact captured mass M_T exceeds 4r², and all its
sources have distance strictly below r from x. Hence

    max_(p in T)||p-x||_infinity < r < sqrt(M_T)/2,
    r<1/2,

so x belongs to E_T in (3). This proves the forward full inclusion. The same
argument without old sources proves the forward residual inclusion.

Conversely, if x in Gamma belongs to E_T, its farthest distance d_T is strictly
less than r_T. Choose one common radius s with d_T<s<r_T. The open square
captures every source in T and has mass at least M_T>4s². Any additional old
or outside atoms captured by that square only increase its mass. Thus T need
not be the **exact** captured set for this implication; forbidding additional
captures would incorrectly shrink the halo. This proves both reverse
inclusions and then (4). All labels still use the same original b1,b2 and
q1,q2 throughout.

A strict closed-square witness gives the same halo: increase its radius
slightly while keeping it below both 1/2 and sqrt(captured mass)/2. The
original boundary atoms are then strictly captured, and any extra atoms help.
This explains closed critical source events without changing strict receiver
inequalities. A point on one label's receiver boundary is not accepted by that
label merely because its closed rectangle contains it. If another strict
witness exists, the forward argument places it in another open menu receiver.

Degenerate source bboxes are allowed, for example a singleton. A **zero-width
receiver** in (3), including equality of its bounds, is empty as a strict
receiver, not a line of halo points. The same rule holds after clipping by
Gamma. If M_T>1, (1) prevents a nonempty radius-1/2 receiver: such a source
bbox would fit a closed unit square with mass greater than 1. The min in (3)
therefore also handles these labels safely. Receiver and crop boundaries may
be closed for Lebesgue area evaluation only after the exact open-set rule is
established; finitely many such boundaries have area zero. Source-space
atomic boundary mass was not discarded in the witness proof.

For completeness, the pure old part simplifies further. If j=max(k,l)>=2,
any Gamma receiver is at farthest distance greater than e+(j-1)h from that
prefix, whereas its mass radius is sqrt(kl)/(2n)<=j/(2n). But

    2n[e+(j-1)h]-j = j-2+2j/n >0.

Thus only the origin singleton can receive in Gamma without outside mass.
Its area is (1/(2n)-e)²=(n-2)²/(4n⁴), with the zero-width n=2 case empty.
The complete menu still retains all pure old rows before this simplification.

## 3. Application to the actual corrected S20 input

Use n=20, K=10 and the fixed actual outside data

    e=1/400, h=21/400,
    q1=(-71/200,0), b1=1/4,
    q2=(0,-171/200), b2=1/2.

The actual source file is
`computation/inputs/r078/n20-corrected-actual-input.json`, a byte-identical
copy of the qualified R77 input, SHA-256
`d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a`.
It contains all 400 original grid atoms, not a collapsed atom, and the two
outside atoms. Its full-grid square caps were proved and exactly checked
in the frozen R77 qualification package. This appendix binds that input;
it does not rerun its global qualification.

The theorem gives 400 old-containing prefix/subset rows (100 pure old,
200 old plus one outside, 100 old plus both), and three pure-outside rows.
These are the mathematical basis for the producer's actual 403-row menu.
Formula (4), not an implementation shortcut, justifies evaluating their
unions and subtracting the complete residual once.

The production result `computation/results/round078-actual-prefix/actual-prefix-area.json`
reports exact directed rational enclosures for the complete actual corner
deletion, approximately 0.007652366201647474, with two area algorithms agreeing.
Its joint-label union after complete residual is approximately
0.0029049224064398374. These numerical executions belong to that production
component and the independent review; this author supplies their mathematical
menu basis and the analytic joint payment below, not another implementation.

## 4. Why exactly nine fine joint prefixes survive

Put d=71/200, D=171/200 and V=sqrt(3)/4, the residual pair radius. A joint
prefix k by l has mass 3/4+kl/400 and vertical bbox span D+(l-1)h. For l>=2,
even its largest possible old mass, at k=10, cannot make it eligible:

    [D+(l-1)h]² > 3/4+l/40.                            (5)

Indeed the difference at l=2 is 3769/160000>0, and its increment from l to
l+1 is at least 2h(D+h)+h²-1/40>0. Thus all joint rows with l>=2 are empty,
without testing receiver samples.

For l=1 the vertical span D dominates the horizontal span d+(k-1)h even at
k=10, and every label is globally eligible. Define

    U_k=sqrt(300+k)/40, k=1,...,10.

In positive receiver coordinates X=-x,Y=-y, its cropped receiver is

    [e,U_k-(k-1)h] x [D-U_k,U_k],                       (6)

with open endpoints understood. The lower X cut is e since U_k>V>d.
All vertical intervals in (6) lie inside the positive Gamma crop.
The k=10 upper X endpoint is below e, so that receiver is empty. For k<=9,

    U_k>V>43/100>8h+e,

so (6) has strictly positive width. Hence the exact active joint prefixes
are (k,1), k=1,...,9, and no others.

Let x_k=U_k-(k-1)h, w_k=x_k-e and t_k=2U_k-D. Their widths decrease strictly:
U_k-U_(k-1)=(1/1600)/(U_k+U_(k-1))<h. Their vertical intervals expand strictly
about the same center D/2. These are nine different fine prefixes with nine
different true radii, not nine choices of a common total local old mass.

The two residual singleton halos miss all these joint receivers. The first
has positive Y at most 1/4, below every D-U_k>D-1/2>1/4. The second has
positive Y at least D-sqrt(2)/4>1/2, above every U_k<1/2. The full relevant
residual is therefore the single pair receiver

    R_0=[e,V] x [D-V,V],                               (7)

which has positive height because 2V>D. It is contained in the first joint
receiver (k=1). It is still the receiver of the same actual outside pair,
with its original total mass 3/4.

## 5. A union payment preserving all nine prefixes

Because the widths decrease and vertical intervals increase, adding row k
after rows 1,...,k-1 contributes exactly the two new vertical bands of total
height 2(U_k-U_(k-1)) and width w_k. Subtract (7) once, since it is contained
in the first row. Thus the actual joint union after complete residual has
the exact analytic area

    J_actual = (U_1-e)(2U_1-D)-(V-e)(2V-D)
       +sum_(k=2)^9 2[U_k-(k-1)h-e](U_k-U_(k-1)).       (8)

This accounts for overlap of the nine labels rather than summing their
individual surviving areas. Every increment in (8) is positive. In fact all
nine labels are needed for this joint union up to positive area: for an
intermediate k choose X between x_(k+1) and x_k and Y between U_(k-1) and
U_k. Earlier rows fail the Y test and later rows fail the X test. Those
patches also lie above the residual. For k=9 use X in (e,x_9), and for k=1
use X in (V,U_1), Y in (D-V,V), beyond all later x_k because
x_2<U_9-h<43/100<V. Therefore no
positive-surviving joint label was silently dropped.

There is a simple strict rational payment for (8):

    J_actual < 9/14080+161/68800
               =9019/3027200 <3/1000.                 (9)

To prove it, use U_9<11/25 and V>43/100, both checked by squaring positives.
The first difference in (8) is

    1/800-(D+2e)(U_1-V).

Since U_1²-V²=1/1600 and U_1+V<22/25, its strict upper bound is
1/800-(43/50)(1/1408)=9/14080. For the eight remaining bands,

    sum_(k=2)^9 2w_k(U_k-U_(k-1))
       =sum_(k=2)^9 w_k/[800(U_k+U_(k-1))]
       <[8(11/25-e)-36h]/[1600(43/100)]
       =161/68800.

The strict gap from the last rational sum to 3/1000 is 313/15136000>0.
Only original captured-prefix mass increments of 1/400 and their decreasing
horizontal widths have been charged. Neither the total old mass 1/4 nor the
outside pair's residual area is charged nine times.

This is a finite joint-active payment for this actual input, not a general
joint-active bound for arbitrary positions or masses. The staircase geometry
in (6) and the common residual containment in (7) are the precise structure
used by the payment. They are not assumed for a general prefix menu.
The separate core proof `proofs/r078-active-prefix-radius-contraction.md`
treats general radius contraction and the obstruction to deleting rows without
geometric justification; it is not a premise of the menu or payment proved here.

## 6. Actual scalar checks, source binding and limits

`computation/round078-kernel/check_joint_staircase_scalars.py` reads
`computation/inputs/r078-kernel/n20-joint-staircase.json` and the declared raw
402-source input. It checks the source hash and every original coordinate,
mass and S index against the formula above. Exact Fraction checks certify
the inequalities excluding l>=2 and k=10, positive widths for k<=9, radius
ordering, and every rational constant in (9). The actual Python 3.15.0 run
passed, exit 0; output and receipt are in
`computation/results/round078-kernel/checks.json`.

This scalar run does not rebuild the producer's 403 rows, approximate roots,
integrate any halo union, rerun global caps, or independently reproduce its
area endpoints. The exact identities and continuum menu proof above carry
the theory; production and independent geometry review have their own
source-bound execution records. All older frozen components are unchanged,
and this component creates no agent, uses no MATLAB and makes no commit.
The unrestricted joint-active corner budget and the general weak constant
remain open.
