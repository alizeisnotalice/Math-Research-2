# R76: the triple-null branch is below the half-mass corner threshold

2026-11-16. Threshold 1, centered axis-parallel common-radius squares. Put
e=m=1/9 and use the complete S3 grid {0,4/9,8/9}², each atom of mass 1/9.
Let the same finite outside measure nu satisfy every global closed square cap

    (sigma_3+nu)(Q)<=max{1,side(Q)²}.

Suppose its effective SW window P=(-1,7/18)² contains exactly two positive
outside atoms b1 delta_q1 and b2 delta_q2. Their positions are arbitrary in
P\[-1/9,1]² and need not have any symmetry. Outside P, nu may have further
atoms. Let Gamma=(-1/2,-e)² and Delta=H(sigma_3+nu)\H(nu).

The following resolves a genuine branch of the arbitrary-position problem:

    If the full receiver of {0,q1,q2} misses Gamma, then
    |Delta intersect Gamma|<=5/108<1/18.                (1)

No position or mass scan is used to establish (1). The remaining potential
violation must use a strict common witness capturing both outside atoms and
the internal origin.

## 1. Complete local labels and a clipped pair annulus

All strict full witnesses have radius below 1/2 by the global caps. For a
receiver in Gamma their source square lies strictly in P. The only internal
grid source there is the origin, of mass m=1/9. Thus locally there are all
seven nonempty full source labels and all three nonempty residual labels.
The latter are the two outside singletons and their pair. None is dropped.

For a fixed outside atom b delta_q, put

    u=sqrt(m+b)/2, v=sqrt(b)/2.

Origin and q have bbox span less than 1 because q lies in P. They fit a closed
unit square, so the full global caps give m+b<=1 and u<=1/2. The complete
mixed pair receiver is Q_0(u) intersect Q_q(u), with open squares. Its
contribution to Gamma after subtracting the outside singleton halo has area
at most

    B(u,v)=(u-e)_+²-(v-e)_+².                          (2)

For clarity, the residual subtracted in (2) is exactly H(b delta_q)=Q_q(v),
the true all-radius singleton halo. Intersecting Q_0(u) with Gamma gives the
square I² with I=[-u,-e] when u>e (otherwise it is empty). Its side is L=u-e.
For either coordinate let A_s=length(I intersect [q_j-s,q_j+s]). Since L<=u,

    0<=A_v<=A_u<=L,       A_u-A_v<=u-v.                (3)

If q_j lies in I, its radius-u interval covers I and its radius-v interval
overlaps length at least min(v,L). If it lies outside I, its overlap is
min(L,(s-d)_+) for a fixed distance d>=0 and is 1-Lipschitz in s. This proves
(3). The product difference is therefore at most

    L²-(L-(u-v))_+²=(u-e)²-(v-e)_+²,

which proves (2), retaining the same outside mass in u and v.

Let F(s)=(sqrt(s)-e)_+². Its derivative is 0 below e² and
1-e/sqrt(s) above e²; these meet continuously at zero at s=e². On
0<=s<=1/4, therefore, 0<=F'(s)<=1-2e. Since u²-v²=m/4,

    B(u,v)<= (1-2e)(u²-v²)=7/324.                     (4)

Receiver edges have zero Lebesgue area; closed limits are used only for
area evaluation, not to turn equality into a strict halo witness.

## 2. A true deletion payment when the triple misses the corner

The internal singleton's halo in Gamma has area

    (sqrt(m)/2-e)²=(1/6-1/9)²=1/324.                   (5)

The full local halo is the union of its seven exact receivers. Removing the
complete residual union leaves only the internal singleton, the two mixed
internal/outside pair receivers, and the triple receiver as possible new
contributions. If the latter misses Gamma, only the first three remain.
For each pair the actual residual H(nu) contains its own singleton halo, so
subtracting just that singleton gives the safe upper bound (4). Thus

    |Delta intersect Gamma|
       <=1/324+2(7/324)=15/324=5/108.

The gap below the requested m/2=1/18 is exactly 1/108. The two pair contributions
are explicitly counted twice; there is no assertion that an arbitrary number
of such label charges can be paid by one copy of m. The proof uses exactly
two pair labels and does not apply this sum when the triple is active. Full
residual credit from the outside pair, other outside labels, or sources
outside P can only reduce this upper bound.

## 3. An exact algebraic test for the remaining branch

Write the bbox of {0,q1,q2} as [a,b]x[c,d], and put

    delta_Gamma=max{(b-a)/2,(d-c)/2,b+e,d+e,
                    -1/2-a,-1/2-c,0},
    L=m+b1+b2.

Then exactly

    the triple receiver intersects Gamma
       iff L>4 delta_Gamma².                           (6)

Indeed delta_Gamma is the minimum of the farthest source distance over the
closed receiver square [-1/2,-e]². In one coordinate the source-containing
center interval [b-r,a+r] meets [-1/2,-e] precisely when
r>=max{(b-a)/2,b+e,-1/2-a}; take the maximum of the two coordinate values.
If sqrt(L)/2 is strictly greater than this minimum, a minimizing center can
be moved slightly into the open Gamma while retaining the strict distance
margin. Equality gives no strict receiver. This proves (6) with all ties.

Whenever (6) holds, the triple is a genuinely eligible label. The global
square caps then imply L<=1: otherwise its source bbox fits a square of side
s<sqrt(L), whose mass L exceeds max{1,s²}. Thus every possible violation of
the half-mass bound must satisfy the necessary conditions

    4 delta_Gamma² < m+b1+b2 <=1.                       (7)

The upper bound in (7) is conditional on triple activity. It is not asserted
for arbitrary empty triples. In particular delta_Gamma>=1/2 ensures the
triple-null branch for every globally admissible pair of masses. More
generally, if a full-grid capacity calculation supplies b1+b2<=U and
m+U<=4delta_Gamma², all masses in that position pair are covered by (1).

This is a decision filter for a remaining finite geometry, not a replacement
for full-grid capacity checks. The triple-active region of (7) is the
unresolved part: there the common outside-pair residual and global caps
involving internal atoms outside P must still be used jointly. No
symmetrization of arbitrary outside positions has been assumed.

## 4. Minimal exact arithmetic checks and scope

The script `computation/round076/check_triple_null_branch.py` reads
`computation/inputs/r076/triple-null-branch.json`. A Python 3.15.0 Fraction run
actually checked (4)-(5)'s exact constants, six rational parametrizations of
u²-v²=m/4 with u<=1/2, and three exact bbox activity examples including an
equality case in (6). It performed no position search, full-grid cap
qualification, halo-area integration, or large source enumeration. The
examples test algebraic activity only, not admissibility. Output and receipt
are in `computation/results/round076-triple-null/checks.json`.

The continuum geometry and all-radius argument are the proof above, not
inferences from these nine checks. This branch does not yet prove the
arbitrary-two-outside bound 1/18, and supplies no counterexample to it. It
instead identifies an exact common-triple witness as a necessary feature
of any remaining counterexample. All R74/R75 frozen files remain unchanged;
no MATLAB, subagent, commit, reset or push is used by this component.
