# R79: all thin joint prefixes share one residual-pair band

2026-11-16. This is a uniform-grid structural bound, separate from the general-block profile obstruction. It keeps the original internal prefixes, a common outside measure and the complete residual. It does not prove the unrestricted two-outside corner target.

## 1. Actual negative-axis setting

Let n>=3, e=1/n², h=(n+1)/n², K=floor(n/2), m_P=K²e. The original sigma_n is the mass1 uniform n-by-n grid. In its southwest source window suppose the two actual positive outside positions are

    q1=(-a,0), q2=(0,-d), 0<a,d<1,
    b1,b2>0, B=b1+b2, v=sqrt(B)/2.

The same original full measure satisfies all closed-square caps for maximum eligible mass1. Other outside sources may lie outside this corner's window. The old origin and both outside atoms fit in a closed unit square, so B+e<=1 and v<1/2. Actual corner witnesses have r<1/2 and only use these local sources.

Use positive receiver coordinates Gamma+=[e,1/2]². The actual joint label for a k-by-l old prefix has radius

    r_kl=sqrt(B+kl e)/2

and, when active, receiver

    J_kl=Gamma+ intersect(a-r_kl,r_kl-(k-1)h)
                            x(d-r_kl,r_kl-(l-1)h).

Open endpoints express strict witnesses; closing finitely many receiver boundaries is permissible for area only. Call a joint label thin if k=1 or l=1. These2K-1 labels use old mass at most K e, not the whole local mass K²e.

## 2. One band contains their union after the true residual

For any0<=s<=1/2, define the cropped pair rectangle

    P_s=Gamma+ intersect(a-s,s)x(d-s,s).

This is the receiving rectangle of the two outside positions at a common radius bound s. The sets P_s increase with s. The true pair residual is P_v and is contained in the complete H(nu). Set

    u=min(1/2,sqrt(B+K e)/2).

Then v<=u<=1/2. Every thin actual joint receiver lies in P_u: its old mass is at most K e, its actual radius is below1/2, and dropping the old source-position constraints can only enlarge its receiver. Therefore

    (union_thin J_kl)\H(nu) subset P_u\P_v.         (1)

This is one common band for the union, not2K-1 separate area payments. It keeps the same B,a,d and the same true pair residual.

For fixed a, the one-dimensional side length of P_s is

    ell_a(s)=[s-max(e,a-s)]_+.

It is continuous and piecewise linear. Where it is positive, its derivative is1 or2, and0<=ell_a(s)<=(s-e)_+. The same holds for d. Hence the area p(s)=ell_a(s)ell_d(s) is absolutely continuous and, almost everywhere,

    0<=p'(s)<=2ell_a(s)+2ell_d(s)<=4(s-e)_+.

Integrating from v to u gives

    |P_u\P_v|<=2[(u-e)_+²-(v-e)_+²]
               <=2(1-2e)(u²-v²)
               <=(1-2e)K e/2=:T_n.                 (2)

The middle inequality is the clipped-square estimate: F(z)=(sqrt(z)-e)_+² has derivative at most1-2e on0<=z<=1/4. The final inequality uses u²-v²<=K e/4, including truncation at1/2. Empty side intervals, v<e and equality of an interval's endpoints are all covered by the positive-part formula.

Thus the entire thin joint contribution has a bound of order1/n. No deletion-area limit or Borel transfer is being asserted: this is an explicit finite-n set inclusion and exact integral inequality.

## 3. Combining with the higher-prefix far-face envelope

The R78 contraction lemma proves that every active joint label with k,l>=2 has both near faces strictly below v. Its surviving area is confined to the common far-face strips inside[e,v]². With

    t=(v-e)_+,
    A=min(t,(a-v-e)_+), D=min(t,(d-v-e)_+),
    F_far=t(A+D)-AD,

the union of ALL higher-prefix joint receivers after complete residual has area at most F_far. Splitting the actual joint menu into thin and higher prefixes and using(1)-(2) therefore yields

    |(union_alljoint J_kl)\H(nu)|<=T_n+F_far.        (3)

The groups may overlap, so this is an upper bound, not an additive identity. The common residual is not optimized independently in either term. Its exact union can give a smaller bound. When a,d<=v+e, F_far=0 and the joint area is at most T_n. This condition is not valid for all globally feasible negative-axis inputs; the R78 far-face example shows why its omission is substantive.

## 4. What this does and does not pay for the full corner

Pure-old corner deletion has area at most

    A0=(n-2)²/(4n⁴).

For each single-outside group, the R76 clipped-annulus argument, with radius truncated at1/2, gives upper bound(1-2e)m_P/4. It is valid for the actual original input even if the auxiliary collapsed single-outside input would not preserve full caps. Applying the union bound to these two single-outside groups and the complete actual joint group gives

    |Delta intersect Gamma|
       <=A0+(1-2e)m_P/2+T_n+F_far.                  (4)

This coarse combination can charge overlapping old contributions more than once. Its usefulness must be checked against the desired budget rather than inferred from each separate term.

In the conditional class a,d<=v+e, equation(4) becomes

    Q_n=A0+(1-2e)e K(K+1)/2.                        (5)

For odd n it simplifies exactly to

    Q_n=1/8-1/(8n²)-1/n³+5/(4n⁴)<1/8.              (6)

The strict sign follows from -n²-8n+10<0 for n>=3 after multiplication by8n⁴. Thus odd uniform grids with such a negative-axis pair do have a full actual corner payment below1/8, allowing joint-active witnesses. This is a conditional structural class, not arbitrary two-outside geometry.

For even n, however,

    Q_n=1/8+1/(4n)-3/(2n³)+1/n⁴.                   (7)

Its excess above the sufficient per-corner threshold1/8+3/(4n⁴) is

    (n³-6n+1)/(4n⁴)>0 for even n>=4.                (8)

So this particular combination of individually safe bounds does not close even the zero-far-face even-n branch. Equation(8) is an obstruction to the coarse accounting estimate, not an actual input violating the target. Positive far-face terms can create additional loss for either parity. Further improvement must couple the single-outside and joint pieces or use finer original old profiles, rather than add separately maximized upper bounds.

For completeness, suppose n is odd and the SAME globally feasible input satisfies in each of its four D4 corner frames either: at most one effective outside position; or exactly a negative-axis pair obeying that corner's a,d<=sqrt(B)/2+e. In the one-position case the R76 bound is no greater than Q_n. Then four copies of(6), plus the geometric noncorner area3-3/n⁴, give

    Delta<=7/2-1/(2n²)-4/n³+2/n⁴<7/2.               (9)

The input must satisfy all original full-grid caps jointly. No assertion is made that four separately feasible pairs can be pasted together. This remains a uniform-grid, restricted-corner-multiplicity result and does not supply the universal maximum-block condition needed for general C2.

## 5. Evidence and scope

The eight n rows in the author's `computation/round079/general-profile-result.json` verify the exact parity formulas(5)-(8), together with the specified constants of the separate general-block fixture. There is no new source-position or mass search in this note. The all-n band inequality is the analytic proof in Sections1-2, not a numerical extrapolation. The original-profile quadratic envelope and its cell computations are separate components and are not rerun here.

The general-block issue is distinct: arbitrary maximum blocks may have concentrated local mass and unequal corner costs. The companion core note `r079-general-block-profile-direction-and-obstruction.md` states the safe two-sided profile directions and gives a qualified finite example. Even completing every uniform-grid negative-axis branch would leave that broader area-payment obligation open.
