# Fixed-angle audit of the Sprague cuts

## Scope and source

This note independently fills the elementary geometry needed for Gibbs's
Definitions 4--7 and Propositions 3--6 at the *single exact input*
\(\sigma=1.5494^\circ=7747/5000\) degrees.  It does **not** reprove the
prior normalization statement (Gibbs Proposition 2), the orbiform reduction
(Corollary 1), any of the AH/EH reflection cuts, or the final numerical area
claim.

The reviewed source is Philip E. Gibbs, *An Upper Bound for Lebesgue's
Covering Problem*, arXiv:1810.10089 (2018), PDF pp. 6--8 (file pages 6--8),
SHA-256
`df420d8754fd888a2adc8a1392d92f4d2f9ef4789adf263437c0afffad48712f`.
Those pages state the contact assertions and the three propositions, but do
not give the support-strip argument or the endpoint-distance inequalities
below.  In particular, the phrase “further than unit distance from that
point” in Proposition 3 needs a statement uniform in the unknown contact
point.

Coordinates and labels use the project's audited convention:

\[
n_j=(\cos(j\pi/3),\sin(j\pi/3)),\qquad
m_j=(\cos(\pi/6+\sigma+j\pi/3),\sin(\pi/6+\sigma+j\pi/3)),
\]

and every support line is \(u\cdot x=1/2\).  Definition 3 cuts only two
corners of the original hexagon, so here (up to the label convention)
\[
 P=\{n_j\cdot x\le1/2\ (0\le j<6),\quad m_5\cdot x\le1/2,
 m_3\cdot x\le1/2\}.
\]
It is **not** the intersection with all six sides of the rotated hexagon.
The relevant support lines/cross-sections used below are

\[
\begin{array}{c|c}
\text{segment}&\text{normal}\\ \hline
F_3F_2&m_2\\
B_3B_2&m_0\\
B_1C_3&n_0\\
D_1E_3&n_4\\
D_1C_2&n_5.
\end{array}
\]

In the shorthand of the task, `CB`, `DE`, and `DC` mean respectively
\(B_1C_3\), \(D_1E_3\), and the \(n_5\)-face \(D_1C_2\); the contact
ultimately needed on the last face is \(D_1G\), not merely all of DC.

## Two elementary lemmas

**Support-strip lemma.**  If a compact convex constant-width-one body \(K\)
is contained in \(\{-1/2\le u\cdot x\le1/2\}\), then it meets both boundary
lines.  Indeed, its maximum minus its minimum in direction \(u\) is one,
while both extrema are in an interval of length one, so they equal \(1/2\)
and \(-1/2\).  If \(K\subset P\), a resulting contact belongs to the
corresponding face of \(P\).

**Endpoint lemma.**  For a segment \([a,b]\), if
\((z-a)\cdot(b-a)\le0\), then \(a\) is a nearest point of the segment:

\[
 |z-(a+t(b-a))|^2=|z-a|^2-2t(z-a)\cdot(b-a)+t^2|b-a|^2\ge |z-a|^2
 \quad(0\le t\le1).
\]

The first lemma supplies contacts on the **original-hexagon** faces
\(B_1C_3,E_1F_1,D_1E_3,D_1C_2\) (and their opposite faces).  It applies to
Gibbs's **orbiforms**, i.e. the constant-width objects in Propositions 3--6;
it is not a claim about arbitrary diameter-one sets.

Here is the missing justification of the two cross-section contacts stated
pictorially in Gibbs pp. 6--7.  The C-corner cut gives
\(h_K(m_5)\le1/2\).  Since \(m_2=-m_5\) and \(K\) has width one,
\(h_K(m_2)=1-h_K(m_5)\ge1/2\).  On the other hand, a support-strip contact
\(q\in K\cap B_1C_3\) satisfies \(m_2\cdot q<1/2\).  (It is enough to
check the two endpoints of this segment; the strict interval check is in the
certificate.)  Convexity now gives a point of \(K\) on
\(m_2\cdot x=1/2\).  Its membership in the original hexagon forces it into
the section \(H\cap\{m_2\cdot x=1/2\}=F_3F_2\).

Likewise the E-corner cut gives \(h_K(m_3)\le1/2\), hence
\(h_K(m_0)=1-h_K(m_3)\ge1/2\).  The forced contact on the original side \(E_1F_1\)
has \(m_0\cdot x<1/2\), so convexity gives
\(K\cap B_3B_2\ne\varnothing\).  Thus the two claims
\[
K\cap F_3F_2\ne\varnothing,\qquad K\cap B_3B_2\ne\varnothing
\tag{C}
\]
are now proved; neither assertion incorrectly treats \(F_3F_2\) or
\(B_3B_2\) as a face of \(P\).

## What the interval certificate proves

Run from this directory:

```sh
/Users/zzh/.venvs/main/bin/python -B lebesgue-universal-cover/computation/certify_sprague_fixed.py \
  > lebesgue-universal-cover/computation/028-sprague-fixed.json
```

`certify_sprague_fixed.py` uses copied `rational_trig.py` (SHA-256
`97e0076b911f51a9c532e14cacb9fb6e3ae7e68707a3c5af77048f9c6623e87f`) and
outward rational intervals, rounded after every operation to denominator
\(10^{24}\).  The actual run exits zero.  It establishes the following
facts without floating-point tests:

Write \(\operatorname{Fan}(A;\gamma)\) for the union of the closed segments
from \(A\) to points of an arc \(\gamma\).  The phrase “closed shape bounded
by” in Definitions 4--6 is made unambiguous here by the following exact
bounded-component decompositions (not by the interval bounding boxes):

\[
\begin{aligned}
C_S&=\operatorname{Fan}(C_2;\widehat{KG}_{F_3}),\\
E_S&=\operatorname{Fan}(E_2;\widehat{IJ}_{B_3})\cup
     \operatorname{Fan}(E_2;\widehat{JH}_{C_3}),\\
A_S&=\operatorname{Fan}(A_1;\widehat{QX}_{E_3})\cup
     \operatorname{Fan}(A_1;\widehat{XP}_{G}).
\end{aligned}
\]

Here every hat is the selected minor arc.  Cutting the displayed Jordan
boundary at the diagonal from the apex to its circle--circle crossing gives
exactly the two fans in the latter two identities; the diagonal is traversed
once in each direction and is not part of the outer boundary.  Thus the fan
decomposition covers the complete defined region, rather than merely a
collection of boxes around portions of it.  The segment-incidence, selected
root, minor-arc, and opposite-root-outside checks in the certificate are the
coordinate verification of this branch choice.  More specifically, it proves
the selected \(J\) and \(X\) strictly satisfy all six original-hexagon
halfplanes and both cut halfplanes.  It also proves
\[
\begin{array}{c}
[E_2I,E_2J],\ [E_2J,E_2H],\ [E_2I,E_2H]\text{ have positive cross products},\\
[A_1Q,A_1X],\ [A_1X,A_1P],\ [A_1Q,A_1P]\text{ have negative cross products}.
\end{array}
\]
Thus in each case the diagonal ray is strictly between the two outer rays
through a total angle strictly below \(\pi\); the two visible fan interiors
are disjoint and meet only on the diagonal.  This supplies the geometric
content of the two union identities.

1. `K,G,I,H,P,Q` lie on the six source segments, and the selected roots
   `J,X` are strictly in \(P\): all eight defining halfplane slacks have
   positive lower endpoints.  The other `J` has \(x>1/2\), while the other
   `X` lies below the bottom of the original hexagon.
2. The five circular pieces are minor arcs.  The script certifies their
   endpoint dot products, all segment-Voronoi endpoint signs, and the
   strict external-disk margins at their fan apexes.  Exact zero endpoint
   signs are supplied algebraically by the tangent/projection definitions,
   not inferred from a rounded zero.  The expected tangent endpoint is fixed
   before checking: K, I, Q are endpoint 0 of CS, ES-B, AS-E respectively;
   H and P are endpoint 1 of ES-C, AS-G.
3. The five fans have exhaustive interval boxes (the minor arc is covered by
   radial projection of its chord on 128 rational parameter cells).  The
   boxes separate \(C_S,E_S,A_S\), so their closed regions are pairwise
   disjoint at this fixed angle.  This is a continuum interval cover, not
   sampling.
4. On the DC subsegment, it certifies

\[
 |C_2-F_3|^2-1>0,\quad (G-F_3)\cdot(C_2-F_3)-1>0,
\quad (z-F_3)\cdot(F_2-F_3)<0\quad(z=G,C_2).
\]

The first two signs imply by the quadratic identity in the endpoint lemma
that every point of \((G,C_2]\) has distance strictly more than one from
\(F_3\); the last signs put the whole segment in the Voronoi half-plane of
\(F_3\) relative to \([F_3,F_2]\).

For clarity, the certificate's arc-to-region step uses the following small
analytic observation.  If a minor unit-circle arc with endpoint vectors
\(u_0,u_1\) is swept from an exterior apex \(A\), and
\((A-c)\cdot u_i\ge1\), then the whole minor arc has the same inequality:
the set \(\{u\in S^1:(A-c)\cdot u\ge1\}\) is a connected arc of length
less than \(\pi\).  Every non-boundary point of the fan is then strictly
outside the unit disk about \(c\), since

\[
 |(1-t)(A-c)+tu|^2-1=(1-t)^2(|A-c|^2-1)
   +2t(1-t)((A-c)\cdot u-1)>0\quad(0\le t<1).
\]

For the Voronoi functional the certificate uses 256 rational radial-chord
cells for each minor arc.  Every ordinary cell has a strictly negative
outward interval upper bound.  A cell containing a tangent endpoint is
handled by the exact derivative of \(L(c+w/|w|)\), whose interval has the
strict sign directed into the negative half-plane.  Thus a rounded endpoint
zero is not treated as a sampled sign.  Consequently the certificate gives
the exact inclusions

\[
\begin{aligned}
C_S^\circ&\subset \operatorname{Vor}_{F_3}([F_3,F_2])
                 \cap\{|z-F_3|>1\},\\
E_S^\circ&\subset
 [\operatorname{Vor}_{B_3}([B_3,B_2])\cap\{|z-B_3|>1\}]\\
 &\hspace{11mm}\cup[\operatorname{Vor}_{C_3}([B_1,C_3])\cap\{|z-C_3|>1\}],\\
A_S^\circ&\subset
 [\operatorname{Vor}_{G}([G,D_1])\cap\{|z-G|>1\}]\\
 &\hspace{11mm}\cup[\operatorname{Vor}_{E_3}([E_3,D_1])\cap\{|z-E_3|>1\}].
\end{aligned}
\]

These are the all-points statements missing from the source's pictorial
argument; they are stronger than checking distances at finitely many drawn
vertices.

## Propositions 3--5 and the DC endpoint issue

Let \(K\subset P\) be a unit-width orbiform.

For Proposition 3, take \(q\in K\cap F_3F_2\).
For \(z\in C_S^\circ\), the first inclusion and the endpoint lemma imply
\(|z-q|\ge|z-F_3|>1\), contradicting \(\operatorname{diam}K=1\).

For Proposition 4, obtain a contact point \(q_B\in B_3B_2\) and
\(q_C\in B_1C_3\).  Each \(z\in E_S^\circ\) is farther than one from the
entire first contact segment or from the entire second contact segment, by
the displayed union and the endpoint lemma.  It is therefore farther than
one from \(q_B\) or \(q_C\), again impossible.

For Proposition 5 one must not silently replace the support contact
\(K\cap[D_1,C_2]\ne\varnothing\) by a contact on \([D_1,G]\).  If a contact
were in \((G,C_2]\), item 4 above and the already forced contact on
\(F_3F_2\) would give two points of \(K\) at distance greater than one.
Hence some \(q_G\in K\cap[D_1,G]\) exists.  Together with
\(q_E\in K\cap[D_1,E_3]\), the third displayed inclusion proves the
Proposition 5 exclusion exactly as in Proposition 4.

This also addresses the closure subtlety relevant to the later contact
argument.  Merely knowing \(K\cap C_S^\circ=\varnothing\) does *not* remove
the boundary segment \([G,C_2]\) from
\(S=\overline{P\setminus(C_S\cup E_S\cup A_S)}\).  The separate strict
distance argument above is what excludes a DC contact in its relative open
part.  Equivalently, a convex body with nonempty interior that has such a
support contact enters the adjacent interior of \(C_S\); the distance proof
is the more direct formulation and does not rely on an unstated boundary
convention.

## Lemma 1 and Proposition 6

The interval boxes make the three closed fan regions disjoint; each is the
closure of its interior.  Let \(U=C_S\cup E_S\cup A_S\), and let
\(S=\overline{P\setminus U}\).  Do **not** replace this closure by
\(P\setminus\operatorname{int}U\): that equality can fail when a cut meets
the boundary of \(P\).  Instead, Propositions 3--5 give
\(K\cap\operatorname{int}U=\varnothing\).  A unit-width convex body is
full-dimensional, so choose \(z\in\operatorname{int}K\).  Since
\(U=\overline{\operatorname{int}U}\), the open set \(\operatorname{int}K\)
cannot meet \(U\), for otherwise it would meet \(\operatorname{int}U\).
Thus \(\operatorname{int}K\subset P\setminus U\).  For every \(q\in K\),
the points \((1-t)q+tz\) lie in \(\operatorname{int}K\) for \(0<t\le1\)
and converge to \(q\) as \(t\downarrow0\).  Hence \(q\in S\), proving
Gibbs Lemma 1 at \(\sigma=1.5494^\circ\).

Combining this lemma with the **source's** Proposition 2 and Corollary 1
proves the stated fixed-angle instance of Proposition 6.  The latter two
results were read for scope only, not independently re-certified here;
accordingly this audit does not claim a new unconditional universal-cover
theorem beyond the source's reduction chain.

## 本项目第28轮衔接

第26轮 proofs/026-base-universal-cover.md 已独立证明直径不超过1的集合可补全为定宽1凸体，并将其放入这里的 P。因此本轮引理无须再依赖来源的 Proposition 2 与 Corollary 1，即可推出 S 含有每个直径不超过1的集合的等距像。后续 AH/EH 切除与最终凸覆盖仍须另行证明。
