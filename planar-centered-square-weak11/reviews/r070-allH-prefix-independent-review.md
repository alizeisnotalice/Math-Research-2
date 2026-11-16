# R70 independent review: all-H compatibility and complete tied prefixes

2026-11-15. Reviewer: kernel-sources worker, independently assigned by root.
This is separate from this worker's author packet for the balanced-square
payment class; that six-file packet remains frozen and is not reviewed here.
The object here is the all-H necessary condition supplied by root/core and
its exact reduction to distance prefixes. Both statements PASS.

## 1. The necessary condition, with all strict boundaries retained

Let A be a finite positive atomic input at strict threshold1 with centered
axis-parallel common-radius squares. Let S be any maximum-mass eligible
source subset, L=m(S)>0. Let T be an eligible mixed subset meeting S and
write U=T intersection S, m=m(T)<=L. Fix x in its OPEN outer receiver E_T(m):

    ||x-q||_infinity<sqrt(m)/2 for every q in T.          (1)

For ANY nonempty H subset S\U, put M=m(H)>0 and
c_H(x)=max_(p in H)||x-p||_infinity. The independently checked conclusion is

    M+m>L  implies  c_H(x)>sqrt(M+m)-sqrt(m)/2.         (2)

Proof. Since H misses U and H subset S, H is disjoint from T. Thus
m(H union T)=M+m; double-counting captured sources is excluded by this
hypothesis. Suppose M+m>L but the conclusion of(2) fails. Every coordinate
span of H union T is controlled by one of three endpoint-pair types:

* H-to-H: its span is at most d(S)<sqrt(L)<sqrt(M+m).
* T-to-T: eligibility (or(1)) gives span<sqrt(m)<sqrt(M+m).
* H-to-T: by (1) and the assumed non-strict reverse of(2),

      ||p-q||_infinity <= c_H(x)+||x-q||_infinity
                       < c_H(x)+sqrt(m)/2 <= sqrt(M+m).

Therefore H union T is strictly eligible with mass M+m>L, a contradiction.
H itself need not be eligible and need not be bbox-closed. No restriction on
source rank or number of labels was used. The exact equality
c_H=sqrt(M+m)-sqrt(m)/2 is FORBIDDEN because(1) is strict. By contrast,
M+m=L does not contradict maximality and MUST NOT activate(2).

## 2. Exact prefix equivalence, independent of geometric realizability

This reduction is an algebraic finite-set lemma. Let P be any nonempty
finite set of positively weighted points, x any point, and m>0,L real.
Write d_p=||x-p||_infinity. For every distinct attained distance r, form

    P_r={p in P:d_p<=r}, M_r=m(P_r).

These are COMPLETE distance prefixes: every source tied at r is included.
There are at most |P| such nonempty sets. The following are equivalent:

(A) For all nonempty H subset P, if m(H)+m>L, then
        max_(p in H)d_p>sqrt(m(H)+m)-sqrt(m)/2.
(B) The same implication holds for every complete prefix P_r.

Proof. (A) implies(B) because each prefix is a subset of P. Conversely,
if(A) fails at H, set r=max_(p in H)d_p. The prefix P_r contains H and
has THE SAME attained maximum distance r. Positive weights give
M_r>=m(H), while z -> sqrt(z+m)-sqrt(m)/2 is increasing. Hence

    M_r+m >= m(H)+m > L,
    r <= sqrt(m(H)+m)-sqrt(m)/2
      <= sqrt(M_r+m)-sqrt(m)/2.

Thus(B) fails at P_r. This proves the equivalence, including all equal-
distance cases and all equality boundaries. Applying this lemma with
P=S\U proves the requested compression of(2). It compresses the all-H
necessary conditions for a FIXED x,U,m; it does not assert that satisfying
them realizes any eligible mixed label or proves any area bound.

## 3. Exact rational comparison and an independent evaluation path

For rational r>=0,M>0,m>0, put q=M+3m/4-r². The original radical comparison
is exactly equivalent to a signed rational-square test:

    r <= sqrt(M+m)-sqrt(m)/2
       iff q>=0 AND q²>=r²m.                           (3)

Indeed adding sqrt(m)/2 and squaring preserves order between nonnegative
quantities and gives r sqrt(m)<=q. Squaring this last comparison is valid
only AFTER requiring q>=0. Equality holds iff q>=0 and q²=r²m. The test
therefore uses no fixed-bit approximation and retains the forbidden
geometric equality. The separate trigger M+m>L is also evaluated exactly.

The exhaustive-subset code uses(3). To avoid merely comparing two copies
of the same predicate, the prefix path evaluates the ORIGINAL two radicals
using adaptive exact rational enclosures. At precision k, integer isqrt
constructs lower/upper dyadic bounds for each square root. If r is below
the resulting lower bound for the radical difference, it returns true; if
above the upper bound, it returns false. Otherwise it doubles k. Exact
algebraic equality is detected first by the equality criterion just stated.
For non-equality, converging enclosures must eventually separate; there is
no precision cap or fixed-precision fallback. The maximum k observed in
this run was32, an outcome rather than an acceptance threshold.

## 4. Actual checks, their scope and reproduction

Standalone source: `computation/round070/review_allH_prefix.py`.
Minimal input: `computation/inputs/r070-allH-prefix-review.json`.
Output and receipt: `computation/results/round070-allH-prefix-review.json`
and `round070-allH-prefix-review-receipt.json`.

The final actual run used Python3.13.18, Fraction and integer isqrt only,
seed70071,96 generated small S fixtures with rank3..8, positive rational
masses summing to L, random proper nonempty U, a<m<L, and rational x.
Every nonempty H subset S\U was enumerated. The code checks not only the
existence equivalence, but that EVERY individual violating H expands to a
violating complete prefix of exactly the same radius. A separate explicit
tie fixture forces two equal-distance sources into one prefix. Five scalar
cases include exact radical equality and close values on both sides; three
further cases verify the strict mass trigger at and on both sides of
M+m=L. These generated fixtures test the algebraic equivalence; no ambient
A or realizable mixed T is claimed for them.

Separately, the frozen R69 balanced-eight-source input is read from root
integration at commit71e6a09. The code enumerates all255 source subsets,
finds the13 eligible sets and the unique maximum S, and checks all-H
necessity and prefix equivalence at9 strictly interior rational receivers
for EACH of its four actual mixed eligible labels (36 receivers total).
The witness half-side is chosen rationally strictly between half the bbox
span and sqrt(m)/2, so no floating witness-boundary acceptance is used.
This is a local necessary-condition check, not a full halo/W area run.

Final totals:1703 exhaustive subset rows,347 complete prefix rows,
61 cases with tied distances,32 algebraic cases with a violation, and all
36 real receivers compatible. All checks PASS. PID48490 exit0 and was
checked absent after wait; script0.034531917015556246s, wall0.10304404198541306s.
An earlier preliminary smoke run (before the three explicit mass-trigger
cases were added) also passed; the archived receipt is the final version.
No MATLAB, Wolfram, LP or random search for counterexamples was launched.

From the topic directory, reproduce with an existing Python3 interpreter:

    python3 -B computation/round070/review_allH_prefix.py --input computation/inputs/r070-allH-prefix-review.json --balanced-input computation/inputs/r069-balanced-core-obstruction.json --output /Volumes/data/CodexTemp/T2/r070-kernel/allH-prefix-recheck.json

The independent review proves only the necessary condition and exact
pointwise prefix compression, plus the listed finite checks. A general
3.5 receiver-area bound and the original C2 problem remain unresolved.
