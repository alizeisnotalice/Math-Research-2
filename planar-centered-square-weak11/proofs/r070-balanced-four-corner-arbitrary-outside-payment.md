# R70: a balanced four-corner class with arbitrary outside sources

2026-11-15. Strict threshold1, finite positive distinct atomic sources in R²,
centered axis-parallel common-radius squares. This is a sufficient class for
W and hence Wprime, with no restriction on the number or geometry of outside
sources. It is not an unrestricted3.5 theorem or a new C2 upper bound.

## 1. Statement and scope

Let S be a maximum-mass eligible source subset of a finite input A. Write
L=m(S)>0, s=sqrt(L). Suppose S consists of exactly four sources at the
corners of a square of side ts, each with mass L/4, where

    sqrt(3)/2 <= t <= 9/10.                              (1)

Translation, reflection and exchange of the coordinate axes are allowed;
the square remains axis-parallel. All remaining sources may be arbitrary,
provided S retains its maximum eligible mass property. Then

    |Wprime_A(S)| <= |W_A(S)|
       <= L F(t) <= L(32+sqrt(2))/10 < (67/20)L < (7/2)L,
    F(t)=15/4-5t+5t²+(sqrt(2)-1)(1-t).                   (2)

We prove the stronger bound for the union J_A(S) of all outer receiver
rectangles of eligible labels T meeting S. The R67 definitions give
Wprime subset W subset J_A(S). Labels need not be bbox-closed in our
argument. Thus arbitrary source rank and arbitrarily many labels in A are
covered, although the selected internal block S has exactly four sources.
There is no sum of repeatedly captured label masses.

The R69 balanced eight-source obstruction has t=9/10 and equal-quarter S,
so belongs to this class. Its all-core failure does not prevent (2). More
generally the proof does not use that input's four outside masses, locations
or its13-label classification. For t=9/10, the support-tree cost is2.7L and
the full side-s support union has area3.61L; neither R68 support-only
certificate reaches3.5L there.

## 2. A general captured-subset carrier

Normalize L=s=1. For any nonempty source set U and mass budget m, write

    E_U(m)=intersection_(p in U)(p+(-sqrt(m)/2,sqrt(m)/2)²),
    K_S(r)=intersection_(p in S)(p+[-r,r]²),
    c_S(x)=max_(p in S)||x-p||_infinity.

Fix any mixed eligible T meeting S. Put U=T intersection S, a=m(U),
b=m(T\S)>0, m=a+b. Maximality gives m<=1, and U is proper. The union
S union T has mass1+b>1 and cannot be eligible, so its bbox span is at
least sqrt(1+b). Both S-to-S and T-to-T spans are strictly less than1,
whereas sqrt(1+b)>1. Hence an endpoint pair attaining that span must be
a cross pair p in S, q in T. For x in E_T(m),

    sqrt(1+b) <= ||p-q||_infinity
               <= ||p-x||_infinity+||x-q||_infinity
               < c_S(x)+sqrt(m)/2.

It follows that

    c_S(x)>f_a(m):=sqrt(1+m-a)-sqrt(m)/2.               (3)

This argument uses the actual captured mass a, rather than charging a once
for every label. Since a<=m<=1, elementary differentiation gives

    rho(a):=min_(a<=m<=1)f_a(m)
       =sqrt(3(1-a))/2,          0<a<=1/4,
       =1-sqrt(a)/2,            1/4<=a<=1.             (4)

Indeed f'_a has the sign of3m+a-1, so the minimizer is(1-a)/3 if a<=1/4,
and a otherwise. At the endpoint m=a, allowing b=0 only enlarges the
mixed-label outer bound. Equation(3) is strict, including on the boundary
of the closed K. Since x in E_U(1),

    E_T(m) subset E_U(1)\K_S(rho(a)).                  (5)

Thus, for ANY finite maximum block S, writing I=Halo(S),

    J_A(S) subset I union union_(empty!=U proper subset S)
                      [E_U(1)\K_S(rho(m(U)))].        (6)

No assumption that U is itself eligible is made. The union in(6) is finite
and depends only on S; it covers every actual label and every outside input
for which S is maximum. It is a necessary carrier, not a sufficient
realizability condition for mixed labels. A tighter version retains the
constraint m>4d_U(x)² in(3); near-full-width examples show that even that
S-only carrier does not yield a universal3.5 bound. That separate obstacle
is not a counterexample to (2) or to Wprime<=3.5L.

## 3. Exact union geometry for the balanced square

Now let S={(0,0),(t,0),(t,t),(0,t)}, masses1/4 each. All proper subsets of
size2 or3 have bbox span t and mass at most3/4, so none is eligible under
(1), including at the strict equality endpoint t²=3/4. The exact internal
halo I is the disjoint union (up to boundaries) of four side-1/2 squares
centered at the corners and the central square

    G=(t-1/2,1/2)².

In particular |I|=1+(1-t)². Let

    C=(-1/2,t+1/2)²,
    r1=rho(1/4)=3/4, K1=[t-3/4,3/4]²,
    r2=rho(1/2)=1-sqrt(2)/4, K2=[t-r2,r2]².

All square boundaries are null for receiver area; endpoint changes below
do not change any area statement.

For singleton U the union of E_U(1) is exactly C: the two unit intervals
on each axis overlap since t<1, and all four products are present. Thus
the union of the four singleton carriers in(6) is C\K1, of area

    (1+t)²-(3/2-t)²=5t-5/4.                           (7)

Inside K1 each corner's internal side-1/2 square contributes a corner
square of side1-t; G also has side1-t. They are disjoint since t>3/4.
Consequently

    |I intersection K1|=5(1-t)².                       (8)

For a diagonal pair or any triple U, its bbox is that of S, so
E_U(1)=G subset I. These captured sets contribute no additional area.
For a bottom adjacent pair,

    E_U(1)=(t-1/2,1/2) x (-1/2,1/2).

Its part in K1 but outside K2 is exactly the bottom strip

    (t-1/2,1/2) x [t-3/4,t-r2),

of area(1-t)(3/4-r2). The other adjacent pairs give the corresponding
left, top and right strips. All four strips are disjoint and avoid I:
their longitudinal interval is(t-1/2,1/2), whereas corner-square intervals
end at1/4 or start at t-1/4; their transverse interval is separated from
G because r2>1/2. Opposite strips are separated as well; adjacent strips
are separated since t-1/2>t-r2 and r2>1/2. Hence their total new area is

    4(1-t)(3/4-r2)=(1-t)(sqrt(2)-1).                  (9)

Equations(7)-(9) give the EXACT area of the coarse carrier union on the
right side of(6), namely F(t). This exactness concerns the stated carrier,
not J_A(S) or the full W area of any outside input.

Finally F'(t)=10t-4-sqrt(2)>0 on(1): sqrt(3)>17/10 and sqrt(2)<3/2 already
imply5sqrt(3)-4-sqrt(2)>0. Therefore its maximum is at t=9/10, where

    F(9/10)=(32+sqrt(2))/10<67/20<7/2.

This proves(2). Rescaling restores the factor L.

## 4. Actual small exact checks and reproduction

The standalone Python script `computation/round070/check_balanced_four_corner_payment.py`
reads `computation/inputs/r070-balanced-four-corner-payment.json`. It performs
exact algebraic coordinate-cell geometry at t=sqrt(3)/2,7/8,9/10, checks
all15 internal subsets and all14 proper nonempty captured sets per fixture,
checks the five pure internal labels, singleton complement, four disjoint
new pair strips, diagonal/triple containment and the exact F(t) identity.
It separately checks the derivative identity and the rational radical
comparisons used in the uniform bound. No outside input is generated and
no full halo, LP, optimization or random search is executed.

The actual run used Python3.13.18 and SymPy1.14.0, PID46298, exit0, and
PID absence was checked after wait. All three fixtures passed; script time
0.2406792500114534s, process wall0.482558459014399s. Source, input, complete
output and receipt are pinned by the accompanying integration manifest.
The continuum proof is sections2-3; finite fixture checks do not replace it.

From the topic directory, run (choose any fresh output file):

    /Users/zzh/.venvs/main/bin/python3 -B computation/round070/check_balanced_four_corner_payment.py --input computation/inputs/r070-balanced-four-corner-payment.json --output /Volumes/data/CodexTemp/T2/r070-kernel/payment-recheck.json

General3.5 maximum-block selection, the near-full-width arbitrary-outside
four-corner family, and the precise C2 value remain open. This worker does
not commit or push; root owns publication and remote verification.
