# Gibbs Propositions 1--2: base-cover audit

Scope: this note audits only the reduction from arbitrary planar sets of diameter
at most one to the slanted Pal polygon `P(sigma)`, at the fixed parameter
`sigma = 1.5494 degrees`.  It does **not** audit the later `CS`, `ES`, `AS`,
`AH`, or `EH` cuts.  The source PDF is Philip Gibbs, *An Upper Bound for
Lebesgue's Covering Problem*, arXiv:1810.10089v1, whose local copy has SHA-256
`df420d8754fd888a2adc8a1392d92f4d2f9ef4789adf263437c0afffad48712f`.

## Result and quantifiers

Let `S` be a nonempty bounded subset of `R^2` with `diam(S) <= 1`.  Then there
is **one plane isometry** `g` (translation composed with one rotation; a
reflection is permitted by the problem convention but is unnecessary here) such
that

`g(S) subset P(1.5494 degrees)`.

The rotation/translation is selected once for the whole completed body and is
then restricted to `S`.  It is not a separate rotation or reflection for
different points or pieces of `S`.  This is the quantifier needed to establish a
universal cover:

`for every S, there exists one g, such that g(S) is contained in P`.

Consequently `P(1.5494 degrees)` is a cover in the full-isometry convention.
Since the construction below never uses a reflection, it is also a cover if only
orientation-preserving isometries are permitted.  This latter observation does
not authorize later reflection-dependent deletions under the restricted
convention.

## What the Gibbs source states

On PDF p. 3, Definition 2 calls a plane-diametrically-maximal unit-diameter set
an "orbiform".  Proposition 1 says that every bounded unit-diameter planar set
is contained in one, citing Vrecica (1981); Corollary 1 reduces covering all
sets to covering all orbiforms.  PDF p. 4, Definition 3 defines `P(sigma)` as
the regular inradius-`1/2` hexagon `H` after two sides of a congruent hexagon
rotated by `30 degrees + sigma` cut the `C` and `E` corners.  Proposition 2
asserts that it covers all unit-diameter sets, citing Baez--Bagdasaryan--Gibbs
(2015), rather than proving it.

Thus neither proposition is an independent proof in Gibbs.  The needed
completion theorem and the elementary cut argument are set out below.

## Completion: imported classical theorem, and the exact hypothesis

We use the planar completion theorem in this form:

> **Completion theorem.** Every nonempty bounded `S subset R^2` with
> `diam(S) <= d` is contained in a compact convex body `K` of constant width
> `d`.

For the present application take `d=1`.  This is the required external classic
result, not something established by the sentence in Gibbs alone.  The exact
source cited by Gibbs is S. Vrecica, *A note on curves/sets of constant width*,
> Publications de l'Institut Mathematique **29** (1981), 289--291 (Gibbs PDF
> p. 21, ref. [16]).  The 2015 paper, p. 8, explicitly states the same use:
"any subset of the plane with diameter 1 can be extended to a curve of constant
width 1."

For clarity, the standard Zorn part does not itself prove constant width.  On
the poset of supersets of `S` having diameter at most one, a chain union remains
of diameter at most one, so Zorn produces a complete set `K`.  Maximality gives
the useful identity

`K = intersection_{x in K} closed_ball(x,1)`:

the right side may be adjoined to `K` without increasing its diameter.  Hence
`K` is closed and convex (the right side is an intersection of closed convex
balls), and it is compact because it is bounded.  A complete bounded set cannot
have diameter strictly below one: a point just outside a boundary point and
sufficiently close to it could then be adjoined without making the diameter
exceed one.  The genuinely geometric remaining assertion, that a complete
planar compact convex set is a constant-width-one body, is the planar
complete-set/constant-width theorem supplied by the cited result.  This also
deals directly with `diam(S)<1`; no invalid attempt is made to append a single
distance-one point to `S`.

## Lemma: a width-one body fits in the regular hexagon

Let `K` be a compact convex body of constant width one.  Write `h(theta)` for
its support function and put

`c(theta) = h(theta) - 1/2`.

Constant width gives `h(theta)+h(theta+pi)=1`, hence `c(theta+pi)=-c(theta)`.
For a proposed orientation `theta`, let the three unit normals of the three
opposite-side pairs of `H` be

`n_0=(1,0), n_1=(cos(pi/3),sin(pi/3)), n_2=n_1-n_0`.

A translation
`t` fits the oriented body into the width-one strip perpendicular to `n_i` iff

`n_i dot t = -c(theta+i*pi/3)`.

The equations for `i=0,1` determine `t`.  Since `n_2=n_1-n_0`, the third holds
exactly when

`F(theta)=c(theta)-c(theta+pi/3)+c(theta+2*pi/3)=0`.

But continuity of the support function and oddness of `c` give

`F(theta+pi/3) = -F(theta)`.

The intermediate value theorem therefore produces an orientation satisfying all
three equations.  All three pairs of strip inclusions then hold, so the one
rotation and one translation place `K` in `H`.  This is a self-contained proof
of the regular-hexagon placement used below.

## Lemma: the two slanted corner cuts retain a copy

Fix `0 <= sigma < pi/6`; in particular this includes `1.5494 degrees`.  Let
`H'` be a copy of `H`, concentric with `H`, rotated by `pi/6+sigma`.  Its six
supporting side-lines split `H minus H'` into six disjoint open corner
triangles `A,B,C,D,E,F`, cyclically labelled as in Gibbs Figure 3.  Each
opposite pair, `(A,D)`, `(B,E)`, `(C,F)`, lies in the two exterior open
half-planes of a pair of opposite sides of `H'`.  Those side-lines are distance
one apart.  Hence points in interiors of opposite triangles have distance
strictly greater than one.

For completeness, the asserted triangle description follows directly from
`pi/6 <= pi/6+sigma < pi/3`: at each vertex of `H` exactly one rotated-hexagon
side has that vertex in its exterior side (at the right endpoint it becomes a
degenerate triangle).  That line meets the two adjacent sides of `H`; the six
resulting corner triangles have disjoint interiors.  The projection on its unit
normal proves the strict distance claim above.

Place `K` in `H` by the preceding lemma.  It cannot meet interiors of both
members of an opposite pair, since its diameter is one.  Thus the set of corner
triangles whose interiors meet `K` contains no opposite pair.  Any such subset
is contained in one of the following eight allowed triples:

`ABC, BCD, CDE, DEF, EFA, FAB, ACE, BDF`.

Here the first six are consecutive triples and the last two are alternating. A
multiple of the `60 degrees` rotation symmetry of `H` sends every consecutive
triple into `FAB`, and sends each alternating triple into `FBD`.  Apply that
one rotation to the entire placed `K`.  It remains inside the fixed `H` and
avoids the interiors of both `C` and `E`.  The cut boundary is retained, hence
it lies in

`P(sigma)=H minus (interior(C) union interior(E))`.

Restricting the same isometry to `S subset K` proves the stated result.  This
is the content of Reduction 1, pp. 7--8, of Baez, Bagdasaryan and Gibbs,
*The Lebesgue Universal Covering Problem*, J. Comput. Geom. 6(1) (2015),
288--299; its locally downloaded PDF has SHA-256
`60fe8dcc58a04f2ba420e341ab0864ea74411acee98fb329f53a549488e385a2`.

## Fixed-parameter coordinate cross-check

The script `check_p_sigma.py` uses 80-decimal-digit `mpmath` arithmetic.  It
intersects the two Definition-3 cut lines with the adjacent sides of the
coordinate hexagon and compares the remaining shoelace area with the analytic
base-area formula already recorded in `proofs/009-gibbs-base.md`.  Its JSON
output is retained as `check_p_sigma.output.json`.

This is only a coordinate/branch check; the proof above, not numerical output,
establishes containment.  The later five deletions and Gibbs's final claimed
area are expressly out of scope.
