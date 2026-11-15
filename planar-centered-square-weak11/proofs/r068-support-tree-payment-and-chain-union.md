# R68: support-tree payment and the exact monotone-chain union

2026-11-15. This supplement gives a general structural class of maximum
blocks paid by a constant at most3.5. It does not prove the unrestricted
W/Wprime selection proposition. It also supplies a strict obstruction to
making this support-only tree criterion universal.

## 1. Fixed objects and the R67 interface

Let A be any finite positive atomic input with distinct positions, threshold1
and centered axis-parallel common-radius squares. Let L>0 be its maximum
eligible subset mass and S an eligible block of massL. Put s=sqrt(L),
write w,h for the bbox widths of S, and let

    Q_p=p+(-s/2,s/2)^2,  Z(S)=union_(p in S)Q_p.

Strict eligibility of S gives w,h<s. The source weights enter this
supplement through L and the maximum-block hypothesis; a point's mass is
not counted anew for every label.

Use the full-mass erosion unions W_A(S),W'_A(S) of the published R67
`proofs/r067-support-rank-and-full-erosion-union.md`. Its proved containment
is

    D_A(S) subset W'_A(S) subset W_A(S) subset Z(S). (1)

For clarity, a full eligible closure T has mass m<=L, outside mass
b=m(T\S), and its Wprime label is E_T(m) minus the receiver using the
actual bbox of T\S and massb. Every point of E_T(m) is within sqrt(m)/2<=s/2
of every captured point of S, so the last containment in(1) is independent
of label count. This supplement does not assume W=Wprime.

## 2. General tree lemma

For p,q in S define the nonnegative geometric edge cost

    c_s(p,q)=s*(|p_x-q_x|+|p_y-q_y|)
              -|p_x-q_x|*|p_y-q_y|.                (2)

Since both coordinate differences are less than s,

    |Q_p intersect Q_q|
      =(s-|p_x-q_x|)*(s-|p_y-q_y|),
    |Q_p\Q_q|=c_s(p,q).                            (3)

For EVERY rooted spanning tree G on the support of S,

    |Z(S)| <= L+sum_(edge in G)c_s(edge).            (4)

Proof. Add the root square first, paying L. Order the remaining vertices
so each parent precedes its child. Before adding a child v, its parent's
square is already in the previous union. Hence

    Q_v\previous_union subset Q_v\Q_parent(v).

Its new area is bounded by the corresponding cost in(3). Summing these
finite geometric increments proves(4). This is an overlap-credit bound
on the true spatial union, not a sum of erosion-label masses.

Define kappa_s(S) as the minimum of the finite tree costs. The immediate
selection reduction is

    exists maximum block S with kappa_s(S)<=5L/2
        implies |W'_A(S)|<=|W_A(S)|<=7L/2.          (5)

The external sources and all eligible labels are otherwise unrestricted,
subject only to S remaining a maximum block. This includes arbitrarily
large ranks and thick bboxes. Thus after this reduction an unhandled input
must have kappa_s(S)>5L/2 for EVERY maximum block, including all ties.

For an explicit tree on rational data, no floating sqrt is needed to
certify(5). Let P=sum_edges(|dx|+|dy|), Q=sum_edges|dx*dy|. Its cost is
sqrt(L)*P-Q, and because both sides are nonnegative,

    sqrt(L)*P-Q<=5L/2
      iff L*P^2 <= (5L/2+Q)^2.                     (6)

This rational test certifies any proposed tree. Computing an approximate
MST is optional; its global optimality is unnecessary for a positive
payment certificate. The examples below use L=1 and exact Kruskal weights.

## 3. Exact monotone-chain formula: a uniform constant below3

Suppose the support can be ordered p_1,...,p_n with both coordinates
nondecreasing, after an optional reflection of either coordinate.
Coordinate increments dx_i,dy_i may be zero, but the points are distinct.
Then

    |Z(S)|=L+s*(w+h)-sum_i dx_i*dy_i <3L.           (7)

Proof. In the unreflected coordinates, for j<i, the intersection
Q_i intersect Q_j has the same lower corner p_i-(s/2,s/2), and upper corner
p_j+(s/2,s/2). These intersections are nested as j increases. Their largest
is Q_i intersect Q_(i-1). Therefore the increment when adding Q_i is
EXACTLY c_s(p_i,p_(i-1)), rather than merely bounded by it. The coordinate
increments telescope to w and h, yielding(7). Since w,h<s and all products
are nonnegative, the strict upper bound3L follows. Reflections preserve
areas and the square geometry.

This chain class can have arbitrary rank and both widths greater than
3s/4; it is not the previous thin-axis or rank-at-most-three reduction.
The tree criterion(5) also covers supports which are not monotone chains.

### Centered-star corollary

There is another simple nonchain class inside(5). If S has at most four
sources and one is the geometric center of its bbox, root the star at that
source. Each edge has |dx|<=w/2<s/2 and |dy|<=h/2<s/2. The expression
s*(dx+dy)-dx*dy is increasing in each nonnegative coordinate on this range,
so every edge costs less than3L/4. There are at most three edges. Thus

    kappa_s(S)<9L/4,  |Wprime|<=|Z(S)|<13L/4.

No chain ordering or restriction on the external sources is required.
This is a corollary of the general tree lemma, not a separate computational
example or a claim about four arbitrary support points.

## 4. Exact positive examples inside the remaining R67 domain

All examples have L=1, so s=1. Every position and mass below is rational.
The saved input and check output give the complete subset census; no
floating LP, Wprime computation or numerical area approximation is used.

### Thick chain

S has four points (0,0),(3/10,3/10),(3/5,3/5),(9/10,9/10), each mass1/4.
Add the outside point(11/10,9/10) of mass1/100. S is the unique maximum
block, its two widths are9/10>3/4, and its rank is4. The full input is
eligible-connected and nonterminal, with exact margin

    d(A)^2-(35/33)m(A)=229/1650>0.

The chain's three edge costs are51/100 each, hence

    |Z(S)|=253/100<3.

This is a general Wprime/W/deletion payment certificate for this maximum
block, not a computation of those regions.

### Thick nonchain paid by a tree

S has points O=(0,0),P=(2/5,3/5),Q=(3/5,2/5),R=(19/20,19/20), all mass1/4.
The outside point is(11/10,19/20), mass1/100. Again S is the unique maximum,
connected/nonterminal margin is229/1650, widths19/20>3/4 and rank4. No
coordinate reflection makes this support a monotone chain: O and R force
the coordinate directions to agree, while P,Q remain incomparable.

The tree OP,PQ,PR has costs19/25,9/25,283/400, totaling731/400. Consequently

    |Wprime|<=|Z(S)|<=1131/400=2.8275<3.5.

Independent rational coordinate-cell integration in the code gives the
actual support-union area1063/400. Its coarse bbox area is1521/400>3.5,
so the new tree certificate pays a genuinely larger structural class than
the simple bbox-area test. None of these values asserts an exact Wprime
area.

In both five-source examples S is uniquely maximum because the full set
has mass101/100 but bbox span squared121/100 and is ineligible, while any
proper outside-containing subset has mass at most76/100<1. The external
point joins R by an eligible pair, and S itself is an eligible hyperedge.
All domain conditions are also checked by complete31-subset enumeration.

## 5. A strict limitation: the tree criterion is not universal

Take the two-row support

    S={(7j/24,0),(7j/24,7/8):j=0,1,2,3}.

Give seven points mass1/400 and the corner(7/8,7/8) mass393/400. Add the
outside point(11/10,7/8) of mass1/800. All masses are strictly positive.
S is the unique maximum block of mass1: the full input has mass801/800
and is ineligible, while any proper outside-containing subset has mass
at most799/800. S has width and height7/8, rank8, and the full input is
connected via the heavy corner. Its nonterminal margin is1303/8800>0.

For these eight support positions, the minimum tree cost is exactly21/8.
More generally, with n equally spaced points per row, horizontal spacing
epsilon=w/(n-1), and 0<epsilon<=h<1, any tree has k>=1 cross-row edges. Each
horizontal edge costs at least epsilon, and each cross-row edge costs
h+(1-h)|dx|>=h. Thus total cost is at least

    (2n-1-k)*epsilon+k*h
      =2w+h+(k-1)*(h-epsilon)>=2w+h.

Adjacent horizontal row paths and one vertical edge attain2w+h. Here it
is21/8>5/2, so every tree fails the sufficient test(5). The two full row
bands also give the exact support-union area

    |Z(S)|=(1+7/8)^2=225/64>7/2,

with excess1/64. This supplies an explicit R67-domain obstruction to a
universal support-only or tree-selection rule. It does NOT refute the
Wprime selection candidate: Wprime can be a strict subset of Z, and its
area is not computed here. Grid/cycle overlap and captured-mass information
remain needed outside the new structural class.

## 6. Actual bounded checks and scope

`computation/round068/check_support_tree_payment.py` uses Python3.13.18,
SymPy1.14.0 and exact Fraction arithmetic. Two symbolic zero identities
check the edge algebra and chain telescoping expression. Three fixed
fixtures are checked by complete eligibility enumeration(31,31,511 subsets),
hypergraph connectivity, unique maximum/domain tests, exact Kruskal cost,
all coordinate-reflection chain tests, and a separate rational compressed-
cell rectangle-union algorithm. The final run took about.0272s internally,
exit0; launched Python PID14178 was checked absent afterward. A preliminary
run used the same mathematical inputs/results before output-label/PID
metadata cleanup; there was no expanded search or additional input family.

Inputs and all rational outputs are saved; the receipt records the exact
command, software/source hashes and actual process exit. These finite
checks validate the examples and algebra, while the general union bounds
are proved in§2–3. No MATLAB, new Agent, floating LP, large computation,
commit or push was initiated. Frozen R65–67 files and other workspaces were
not modified. General Wprime selection, a general3.5 upper bound and the
precise C2 matching problem remain open.
