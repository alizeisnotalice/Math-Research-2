# R69: balanced corners refute the core-or-one-fifth alternative

2026-11-15. Strict threshold1, finite positive atomic inputs and centered axis-parallel common-radius squares. The heavy-core upper bound remains valid. The stronger proposed selection alternative is false: failure of every heavy-core budget does not force a maximum block with net mixed increment at most L/5. This note verifies every maximum-block and core quantifier in a single eight-source example. It is not a counterexample to the original3.5 deletion or erosion-union selection problem.

## 1. Input, unique maximum and the remaining support domain

Let S consist of the four corners

    p0=(0,0), p1=(9/10,0), p2=(9/10,9/10), p3=(0,9/10),

each of mass1/4. Add four residual sources, each of mass749/1000:

    q0=(-3/4,-11/50), q1=(28/25,-3/4),
    q2=(33/20,28/25), q3=(-11/50,33/20).

The rotation (x,y) -> (9/10-y,x) cyclically permutes both lists. S has mass1 and two widths9/10, so is strictly eligible and fat. It is the unique maximum eligible source subset, L=1. The following exhaustive classification proves this assertion.

* With no residual source, the only eligible subsets are the four singletons and the full S. Every pair/triple has squared span81/100 and mass at most3/4.
* With one residual source and one S source, only the corresponding nearest pair{pi,qi} is eligible. Its mass999/1000 exceeds its squared span9/16. Every other pair has span at least28/25 and mass999/1000.
* One residual source with any two or more S sources has mass at most749/1000+k/4. For two S points its span is at least28/25, whose square784/625 exceeds1249/1000. For three or four S points its span is at least33/20, whose square exceeds the maximum mass1749/1000.
* Any set containing two residual sources has squared span at least(187/100)^2. If it has exactly two, its total mass is at most2498/1000, below that square. If it has three or four residual sources, its span is12/5 and its mass is at most3996/1000, again below its squared span.

Residual singletons themselves are eligible. There are therefore exactly13 eligible subsets: eight singletons, S, and the four nearest mixed pairs. All are bbox-closed in the full input. None of the mixed sets has mass1, so S is uniquely maximum. S and its four nearest mixed pairs connect the full eligible hypergraph. Moreover

    d(A)^2-(35/33)m(A)=(12/5)^2-(35/33)(999/250)=837/550>0.

The minimum support-tree cost of S is27/10, since each of its three tree edges costs at least9/10 and a side path attains that value. Its full side-1 support union is C=(-1/2,7/5)^2, of area361/100. Thus the old thin, low-rank, tree and filled-support conditions do not pay this support.

## 2. All fifteen core choices fail, and Phi_* is exact

Use the full refined heavy-core functional from this round:

    Phi_S(H)=Gamma_(sqrt(3m(H))/2)(H)
             +min{9(E+beta), |Z_sqrt(E)(S\H)|},
    E=1-m(H),
    Gamma_r(H)=|C|-|C intersect K_H(r)|+|Z_1(H) intersect K_H(r)|.

There are fifteen nonempty cores H. Rotation and reflection reduce their geometry to the following five possibilities. All areas retain C intersect K; in particular an adjacent two-corner core has one side of K outside C.

* One corner: r=sqrt(3)/4<1/2, so K_H is contained in its side-1 support square. Thus Gamma=361/100.
* Two adjacent corners: r=sqrt(6)/4, r^2=3/8. The only part of (C\Z_1(H)) intersect K is a rectangle of area(2r-9/10)(r-1/2). Consequently Gamma=241/100+19r/10. Since r>3/5, Gamma>355/100.
* Two opposite corners: with the same r, the two missing rectangles have total area2(r-1/2)^2. Thus Gamma=236/100+2r>356/100.
* Three corners: r=3/4, K=[3/20,3/4]^2. Its uncovered corner square has side1/4, so Gamma=361/100-1/16=1419/400.
* All four corners: Z_1(S)=C up to boundaries, so Gamma=361/100 and E=0.

These formulas prove Gamma>=1419/400>7/2 for every core, already sufficient to certify failure of all fifteen budgets because the low term is nonnegative. They also determine the minimum of the full refined functional exactly. For a proper core, E>=1/4 and V=S\H is nonempty. Its support union contains at least one square of area E, while9(E+beta)>=E. Therefore

    Phi_S(H)>=Gamma+E>=1519/400>361/100  (H proper).

For H=S the low term is zero, giving

    Phi_*(S)=361/100>7/2.                              (1)

No nearby-source approximation, fixed-bit equality decision or unsuccessful sampled core is used to assert(1). Every core is covered by the symmetry classification.

## 3. A net mixed patch strictly larger than one fifth

Let I=H(S), R=H(A\S), P=H(A)\(I union R). The complete internal halo is the union of four open squares centered at the pi with half-side1/4 and the full-S receiver(2/5,1/2)^2. These five rectangles are pairwise disjoint, and

    |I|=4(1/4)+(1/10)^2=101/100.

The residual halo consists only of the four single-source squares. Since4(433/1000)^2>749/1000, let Rplus be the union of the four closed squares at qi with half-side433/1000. Every internal rectangle is disjoint from every rectangle in Rplus, by direct coordinate separation. Hence

    I intersect R = empty.                             (2)

Set rho=999/2000. The two-atom mass999/1000 strictly exceeds4rho^2=998001/1000000. The open rectangle

    J0=(-999/2000,-501/2000) x (-999/2000,559/2000)

consists entirely of strict witnesses of this one common half-side rho for{p0,q0}. Let J1,J2,J3 be its successive rotations. The four rectangles are pairwise disjoint, and all avoid I: for example J0 has right edge-501/2000<-1/4 and upper edge559/2000<13/20, separating it from the internal singleton and central rectangles.

Each Ji intersects only its own residual outer square, and that intersection has width365/2000 and height57/80. Each Ji itself has width249/1000 and height779/1000. Removing Rplus therefore leaves a subset F of P with exact area

    |F|=4[(249/1000)(779/1000)-(365/2000)(57/80)]
        =255759/1000000>1/5.

By(2), the same number is a strict lower bound for the net mixed increment:

    |P|-|I intersect R| >=255759/1000000> L/5.          (3)

This is a whole-region rational certificate, not a numerical difference of full halo computations.

## 4. What is refuted, and what remains true

There is only one maximum block. Equation(1) verifies every core in its premise, and equation(3) verifies the failure of the proposed net-mixed conclusion. Thus the candidate saying that all-core failure should imply some maximum block with net increment<=L/5 is false, even in the connected, nonterminal, fat support-tree remaining domain. Both the unrefined and refined core functional fail here.

The general heavy-core theorem is unaffected: it supplies an upper bound which is simply too large for this input. Nor is this an example of large W or large true deletion. The complete eligible-label classification in section1 gives the direct bound

    Wprime subset W subset I union union_i E_{pi,qi}(999/1000).

Every mixed receiver has area

    (sqrt(999/1000)-3/4)(sqrt(999/1000)-11/50)
        <(1-3/4)(1-11/50)=39/200.

Consequently

    |Wprime|<=|W|<101/100+4(39/200)=179/100<7/2.        (4)

No exact full-halo area is needed for(4). The balanced example shows that a general method must use more of the actual internal or mixed-label geometry than the worst compact bound plus a fixed one-fifth allowance, even after a heavy-core screen. It supplies no counterexample to the original selection problem and no new general C2 bound.

## 5. Actual exact check and reproducibility

Root ran `check_balanced_core_obstruction.py` with Python3.15.0 and Fraction arithmetic. It enumerated all255 nonempty source subsets, found exactly the13 labels in section1 and the unique maximum S15, and checked connectivity, nonterminal status and the support-tree cost. It integrated the four rational patch rectangles after subtracting the residual outer squares by an exact coordinate-cell algorithm. For all15 cores it separately evaluated Gamma at rational outward upper radii434/1000,613/1000,3/4,867/1000 according to core size. Their squared values bound3m(H)/4 from above; since Gamma is nonincreasing in radius, these are strict LOWER certificates at the optimal permitted radius. All fifteen computed lower bounds are at least1419/400. The exact Phi_* formula in section2 additionally follows from the analytic nonnegative low-term argument above.

The run exited0, script0.005539seconds and wall0.037407seconds, with owned PID34763 checked absent. Source SHA `d7ab753cc15e01dcd7dabc290ee986513f5aeebcb7239b2a15616defb2e2e6a6`. Minimal input, all255 subset rows, all15 core lower bounds, exact output and process receipt are retained. There was no complete halo evaluation, parameter optimization or new random search.

To reproduce, copy the standalone script into a fresh directory under `/Volumes/data/CodexTemp/T2` and run `python3 -B` on that copy. It writes its two JSON outputs beside the copy, requires mounted data, and rejects disabled assertions. No old workspace or third-party source module is imported. The independent review separately binds its actual verification and the final proof version.
