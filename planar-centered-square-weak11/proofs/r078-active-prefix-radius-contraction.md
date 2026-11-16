# R78: active original prefixes contract the near faces of the joint receiver

2026-11-16. This note keeps the original uniform grid. It gives a general-n geometric restriction and an explicit residual-strip bound for joint-active labels of two negative-axis outside atoms. It does not establish the full corner payment for this family.

## 1. Hypotheses and actual joint receivers

Let n>=3, e=1/n², h=(n+1)/n², K=floor(n/2), and let sigma_n put mass e at each point of {0,h,...,(n-1)h}². Fix actual outside atoms

    q1=(-a,0), q2=(0,-d), a,d>0,
    b1>0, b2>0, B=b1+b2, v=sqrt(B)/2.

The complete original input, possibly with other outside atoms outside the corner source window, satisfies all closed-square caps with maximum eligible mass1. Assume these are the only two positive outside positions in the southwest effective window P=(-1,1/2-e)². Thus a,d<1. Write Gamma+=[e,1/2]² in receiver coordinates X=-x,Y=-y.

For an actual anchored k-by-l old prefix, 1<=k,l<=K, define

    p_k=(k-1)h, p_l=(l-1)h,
    r_(k,l)=sqrt(B+kl e)/2.

Call the joint label active when its strict full receiver intersects Gamma. Original caps imply that an active label has B+kl e<=1, hence

    v<r_(k,l)<=1/2.                                 (1)

Indeed a strict receiver makes the label eligible; a mass above1 would violate its enclosing-square cap. An inactive label contributes no area. For an active label its exact receiver in the corner has the following rectangle representation, with open source-distance inequalities:

    J_(k,l)=Gamma+ intersect
          (a-r,r-p_k) x (d-r,r-p_l).                 (2)

The area is unchanged by using closed endpoints. This follows directly from the actual bbox [-a,p_k]x[-d,p_l], not from replacing the old mass by a single atom. If a label's witness also captures more old or outside mass, that only strengthens the full witness; the complete actual-prefix menu remains valid.

The true pair residual receiver is

    R_12=Gamma+ intersect(a-v,v)x(d-v,v).            (3)

It may be empty. It is part of the complete residual H(nu), which also contains the two singleton receivers and any outside-P contributions. Keeping only R_12 to remove points gives a safe upper bound on surviving joint area. These formulas concern one common actual outside pair, with no independent receiver choices.

## 2. Radius contraction beyond the first row or column

For every active joint label,

    k>=2 implies r_(k,l)-p_k<v,
    l>=2 implies r_(k,l)-p_l<v.                     (4)

We prove the first implication, the second being symmetric. If l>=2, activity at a center with Y>=e forces

    r>p_l+e.

Consequently

    r-v=kl e/[4(r+v)]
        <kl e/[4(p_l+e)]
        <(k-1)h.                                   (5)

For the last strict inequality, use k,l>=2 and

    kl<=4(k-1)(l-1), h²>e.

In fact kl e<4(k-1)(l-1)h²<4(k-1)h(p_l+e). No assertion is made that a selected collapsed old mass and the outside pair obey the original capacity rows.

If l=1, concavity of the square root yields

    r-v=[sqrt(B+k e)-sqrt(B)]/2<=sqrt(k)/(2n).

For k>=2, k<4(k-1)² and h>1/n, so sqrt(k)/(2n)<(k-1)h. This proves(4), including the one-row boundary. The case k=1 has no such contraction assertion: those labels can cross X=v. Similarly l=1 labels can cross Y=v. At n=3, K=1, there are no k>=2 or l>=2 labels, and the statements are vacuous rather than an additional n3 payment theorem.

The scalar implications in(4) are not restricted to negative-axis positions. For any exterior label of total mass B>=0, if adding the actual k-by-l prefix has a strict receiver in Gamma and r=sqrt(B+kl e)/2, the same argument applies: activity gives r>p_l+e and r>p_k+e. Neither the global caps nor the outside geometry is needed for this scalar step. The negative-axis assumptions are needed for the exact receiver/residual face formulas(2)-(3), and the global caps supply r<=1/2 for their corner use.

The strict inequalities in(4) follow from actual activity and original prefix spacing. They are stronger than merely bounding the total internal mass by m_P, and hold for unequal axis distances and unequal positive outside masses.

## 3. Only two far faces remain for multirow and multicolumn labels

For k,l>=2, equation(4) places J_(k,l) strictly inside the near-side square

    e<=X<v, e<=Y<v.

Any point of this receiver not in the true residual pair must therefore satisfy

    X<=a-v or Y<=d-v.                               (6)

Away from the boundary lines, these are the only two possible exposed faces. In particular the higher-prefix joint labels cannot create an independent contribution beyond X=v or Y=v; such near-face contributions can only come from k=1 or l=1 respectively.

This gives a uniform geometric area envelope for the union of ALL joint labels with k,l>=2, after subtracting the true complete residual. Define

    t=(v-e)_+,
    A=min(t,(a-v-e)_+), D=min(t,(d-v-e)_+).

If at least one joint label is active, (1) gives v<1/2, so no additional clipping at the outer side of Gamma is needed. Otherwise the joint contribution is zero. Since the two strips in(6) lie in the same square of side t, their union has area

    t(A+D)-AD.                                     (7)

Thus

    |[(union_(k,l>=2) J_(k,l))\H(nu)] intersect Gamma|
         <=t(A+D)-AD.                              (8)

The overlap AD is subtracted once. The outside pair and its residual are common to every prefix; no independent strip optimization or duplicated source mass is used. Bound(8) may be loose and is not by itself a half-mass or1/8 payment. If no multirow/multicolumn label is active, its left side is zero regardless of the right side.

## 4. A conditional reduction to the first-row/first-column menu

If the actual parameters also satisfy

    a<=v+e and d<=v+e,                              (9)

then A=D=0. Every k,l>=2 joint receiver lies inside the residual pair up to receiver-region boundary lines, and contributes zero surviving area. The exact deletion-area calculation may therefore discard those joint labels while retaining all pure-old labels, all single-outside labels, the complete residual and the joint rows

    k=1 or l=1,                                   (10)

of which there are2K-1. This reduces the joint part of the menu from K² to2K-1 for the conditional parameter class(9). It does not reduce the entire full menu to those rows.

The qualification in(9) is essential. For the R77 S20 input, v=sqrt(3)/4, a=71/200 and d=171/200; its second inequality fails. The theorem does not claim that the residual removes every multirow/multicolumn label for all negative-axis parameters, nor does it exclude first-row active contributions. At equality in(9), a line on the closed corner boundary can remain outside the open residual receiver; it has zero area. On the interior of Gamma the residual containment is pointwise.

## 5. A qualified example where the far-face joint labels are indispensable

The conditional reduction in Section4 cannot be made universal, even for a symmetric negative-axis pair. Fix n=20 and

    a=d=1/2, b1=b2=37/100, B=37/50.                 (11)

Keep the complete original400-point grid. The exact full-grid capacity interface from the R77 unequal-axis lemma applies. In this equal-distance case, the joint cap at fixed j=max(k,l) is strongest for k=l=j, and its bound is

    U_j=max(1,[1/2+(j-1)h]²)-j²/400.

The unit rows decrease through j=10. The nonunit rows increase from j=11, since their forward difference is

    h-2h²+(2j+1)(h²-e)>0.

Thus U=min(U_10,U_11)=min(3/4,1197/1600)=1197/1600. The analogous individual bound is V=1/2: its unit rows decrease through10, and its nonunit forward difference h+(2j-1)h²-1/20 is positive. Consequently

    b_i<V, B<U, U-B=13/1600>0.                     (12)

All original square caps hold. Strict slack for outside-containing rows and proper-old mass below1 give unique maximum S20. Each outside atom connects to the origin by a strict eligible edge, and old grid-neighbor edges connect the grid. The whole input has nonterminal gap

    (599/400)²-(35/33)(87/50)=698811/1760000>0.      (13)

Now consider the closed receiver cell

    E=[-61/1000,-59/1000]², |E|=1/250000.           (14)

At every point of E, the common open radius89/200 captures the actual8-by-8 old prefix and both outside atoms. The greatest outside distance is441/1000<89/200; the old maximum coordinate7h=147/400 is also strictly captured, and the next coordinate8h is not. The full mass is9/10, with strict margin1079/10000.

Any square at a point of E capturing either outside atom has radius strictly above439/1000. The entire residual mass is37/50, below4(439/1000)²=192721/250000. Hence E avoids the complete residual. The same radius threshold excludes every thin joint label, meaning k=1 or l=1: its old mass is at most10e, so its total mass is at most153/200<192721/250000. Any single-outside label has mass at most37/100+m_P=31/50, also below that threshold. Pure-old corner witnesses consist only of the origin singleton, whose radius1/40 cannot reach E. Radii at least1/2 are excluded by the original full caps.

It follows that E lies in true deletion but in no full receiver of the menu obtained by dropping all k,l>=2 joint labels. The removed labels are therefore indispensable: deletion area decreases by at least1/250000. This is a strict counterexample to universal first-row/first-column menu pruning, not to any corner area payment. It also proves a positive-area effect rather than merely showing that a high-prefix full label is active.

The production certificate separately evaluates the exact403-label original-prefix menu for this one input, using two directed72-bit union algorithms. It finds true corner deletion approximately.016617159266386683, while the deletion with all k,l>=2 joint labels removed is approximately.004715488935277193; their certified positive difference is approximately.01190167033110949. The exact endpoints, not these decimals, carry the finite-area conclusions. All440 grid-prefix capacity rows plus the three pure-outside rows are checked. The complete403-label approach avoids a fresh45302-closure computation. The positive cell(14) has the analytic proof above; the production run also excludes all322 retained labels over that cell by exact bbox-to-cell distance inequalities.

Here a=d=1/2>v+e, so both far-face widths in Section3 are positive. This is precisely the regime not covered by(9). Neither this example nor its positive removal loss violates the half-local-mass or1/8 single-corner budgets; both displayed actual deletion areas are well below those thresholds.

## 6. Actual scalar checks and scope

The author ran `computation/round078/check_prefix_contraction.py` once with Python3.15.0/Fraction. At n in {3,4,5,9,20,33,65,129}, it verifies5247 sufficient two-dimensional contraction comparisons and123 one-row squared comparisons. It also checks the explicit monotone cap-minimum constants, full-cell witness and exclusion inequalities, area1/250000 and nonterminal gap for(11). PID6040 exited0 and was confirmed absent; kernel time was about.017 seconds.

These are finite exact regressions of the displayed analytic proof, not a proof of all-n validity by sampling. The run does not rebuild the403-label catalog, integrate a halo, execute the producer's440-row capacity check or repeat the independent review's calculations. The source, minimal input, result and receipt are archived together. No general position search is part of this note.

## 7. Consequence for the active-branch strategy

The exact original-prefix menu must still account for the single-outside and first-row/first-column joint pieces, and for the two far-face strips when(9) fails. Equations(4)-(8) identify where the surviving higher-prefix joint contribution can occur and provide a shared-residual geometric bound. They do not complete the general-n two-outside budget, and they do not permit using total m_P in place of each kl/n².

The separate R78 menu proof establishes complete old-prefix/full-label coverage and handles the S20 finite staircase calculation. This note does not repeat that calculation or claim a general theorem from one input. Its remaining question is whether the far-face terms and the single-outside contributions can be paid together by the complete capacity constraints and residual, without adding them independently at their separate maxima.
