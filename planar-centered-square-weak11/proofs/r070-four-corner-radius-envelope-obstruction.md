# R70: a finite obstruction to the radius-aware source-only area budget

Use the notation of `r070-captured-mass-radius-and-prefix-compatibility.md`.
This proof concerns M(S), including its full internal-radius clamp. It
does not by itself test the stronger common-scalar prefix region.

## 1. Four equal corners

Normalize L=1. Put equal masses 1/4 at the four corners of [0,d]²,
where d=99/100. The full set S is eligible and has mass 1. Let

    C=[-1/2,d+1/2]²,
    rho=(sqrt(7)-1)/2,
    K=[d-rho,rho]².

Except for boundaries, the four open side-one squares centered at the
corners have union C: the two one-dimensional intervals overlap because
d<1, and the support is a Cartesian product. For every x in C outside K
there is a singleton U with r_U(x)<1/2. Its mass is a=1/4, so the clamp
parameter q=max{1/4,4r_U(x)²} is below 1. The function

    f(q)=sqrt(3/4+q)-sqrt(q)/2

is increasing for q≥1/4, and f(q)≤f(1)=rho. Since x∉K means
c_S(x)>rho, the exact scalar criterion puts x in M(S). Thus

    |I∪M(S)| ≥ |M(S)| ≥ |C\K|
      = (1+d)²-(sqrt(7)-1-d)²
      = (199/50)sqrt(7)-7 > 7/2.                       (1)

The last strict comparison has the rational certificate

    7(199/50)²-(21/2)² = 791/1250 > 0;

both sides before squaring are positive. No numerical area integration
is needed. This already disproves a uniform 7/2 bound for this necessary
source-only envelope. Its area need not equal a true deletion area.

## 2. Embedding in the actual remaining domain

The same S can be placed in a connected, nonterminal input with a unique
maximizing block. Add just

    q=(33/32,99/200),  m(q)=1/1000.

Every proper S subset together with q has mass at most 751/1000<1.
The only subset containing q with mass at least 1 is the full input;
its span is 33/32 and (33/32)²>1001/1000, so it is not eligible.
Therefore S is the unique maximizing eligible block, with L=1.

The pair consisting of q and either right corner is eligible: its span
is 99/200 and (99/200)²<251/1000. Since S itself is eligible, the eligible
hypergraph of the whole input is connected. The whole input is
nonterminal, since

    (33/32)² > (35/33)(1001/1000).

Both spans of S are 99/100>3/4. Its equal-side support-square union has
area (199/100)²>7/2. In the R68 support-tree criterion an adjacent edge
has cost d and a diagonal edge cost 2d-d²>d. Every three-edge spanning
tree therefore has cost at least 3d, attained by an adjacent-edge tree.
Consequently kappa(S)=297/100>5/2.

Thus the obstruction to the M(S) area budget persists in the actual
fat, connected, nonterminal domain beyond the support-tree payment
criterion. It is not a counterexample to a bound for W′ or D. The
production R70 four-corner pilot evaluates the S-only input and does not
claim to have computed the halo of this five-source embedding.

## 3. Evidence and remaining obligation

The companion Fraction checker enumerates all 31 nonempty subsets of
the five-source embedding, confirms the domain assertions, and checks
the displayed rational strict inequalities. This small domain check is
separate from the production directed-cell area enclosure for M(S).
The latter is a finite computation; the stronger uniform-grid limit is
proved separately in the independent R70 source-only obstruction
attachment. None of these statements gives a lower bound for the
original constant C₂.
