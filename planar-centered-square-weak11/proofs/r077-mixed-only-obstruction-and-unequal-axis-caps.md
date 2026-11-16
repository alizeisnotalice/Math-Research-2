# R77: a mixed-only half-mass obstruction and unequal-axis grid caps

2026-11-16. Strict threshold 1, centered axis-parallel common-radius squares.
All preceding frozen components are unchanged. This note first disproves an
unrestricted **mixed-only** bound, then records exact full-grid constraints
for a separate unequal-axis family. Neither result is an actual uniform-grid
corner counterexample.

## 1. The mixed-only question and a strict counterexample

Let the internal measure be alpha=delta_0, and take two distinct outside
atoms

    beta_1=4 delta_(-d,0), d=28/25,
    beta_2=8 delta_(0,-D), D=161/50.

Write A_1,A_2,A_12 for the complete receivers of the three labels containing
the internal origin and at least one outside atom. Let

    R=H(beta_1+beta_2), O=(-infinity,0)².

The candidate inequality

    |[(A_1 union A_2 union A_12)\R] intersect O|<=1/2

is strictly false. Its exact left side for this example is

    78/25-53sqrt(5)/50-161sqrt(13)/100+161sqrt(3)/50
       >25849/50000=1/2+849/50000.                     (1)

The internal singleton receiver A_0 is not added to this union. Points of
A_0 which also belong to an A_i are naturally included by the stated
mixed-only definition; A_0 is not subtracted either.

Here is the complete true-residual geometry proving (1). Use positive
receiver coordinates X=-x,Y=-y on O, and put

    a=sqrt(5)/2, s=sqrt(2), V=sqrt(3), U=sqrt(13)/2.

The first mixed receiver in O is

    A_1=(d-a,a) x (0,a).

The second mixed receiver is empty: its mass is 1+8=9 but its source span
D exceeds 3. The triple receiver in O is

    T=A_12 intersect O=(0,U) x (D-U,U).                (2)

The full residual in O has all three nonempty labels:

    R_1=(d-1,d+1) x (0,1),
    R_2=(0,s) x (D-s,D+s),
    R_12=(0,V) x (D-V,V).                              (3)

These are exact receivers of the singletons and pair, so every residual
radius is accounted for. The interval in X for the pair receiver has been
clipped at zero because d<V. The triple has the same clipping because d<U.

The geometric orderings needed below are

    1<a<d<s<V<U,   D>U+s,   D<2V.

In particular T lies strictly above A_1 and R_1, and strictly below R_2.
The pair residual R_12 lies inside T and above A_1. Thus A_1 encounters only
R_1, while T encounters only the pair residual. They are separate components
of the mixed deletion, with no source choices changed between them.

The first contribution is

    |A_1\R|=(2a-d)a-(a-d+1)
       =3/2+d-(1+d)a=(131-53sqrt(5))/50.               (4)

The second contribution is

    |T\R|=U(2U-D)-V(2V-D)
       =1/2-D(U-V).                                   (5)

Adding (4)-(5) gives (1). No disjointness of hypothetical annulus payments
is being assumed: the actual receiver pieces and the complete residual have
just been computed, and the separation in Y proves the required additivity.

## 2. Strict rational verification and scope

The ordering and excess use only elementary positive rational square-root
bounds:

    a<d=28/25,
    U<1803/1000,
    V>433/250,
    s<283/200.

Each follows by squaring positive quantities. The last two upper bounds give
U+s<3218/1000<D, and the lower V bound gives D<2V. Formula (4) decreases in
a, so

    |A_1\R|>3/2-d²=307/1250.

Also U-V<71/1000, hence

    |T\R|>1/2-(161/50)(71/1000)=13569/50000.

Their sum is the strict lower bound in (1). These checks certify the proposed
failure without depending on floating-point signs or random sampling.

This counterexample has no global maximum-block caps: each outside atom
alone exceeds the unit mass cap. It also does not automatically fit a grid
corner source window. If its old mass is matched by spatial/mass scaling to
the full local grid mass m_P=floor(n/2)²/n², then sqrt(m_P)>=1/3 for n>=2.
The scaled second outside distance is at least D/3=161/150>1, outside P_SW.
This observation only excludes this particular normalization; it does not
prove a localized or globally capped mixed-only bound for other parameters.

The mechanism is nevertheless decisive for the unrestricted proposed lemma:
one first-pair contribution can be near m/4 while a disjoint triple contribution
exceeds m/4. Removing the explicit old singleton from the R75 formula does not
restore the m/2 bound. Any successful corner proof must retain additional
radius, source-window, global-cap or original internal-position information.

## 3. Exact all-n square caps for two unequal negative-axis positions

This separate interface keeps the full uniform grid. Let n>=2,

    h=(n+1)/n², w=1-1/n²,
    sigma_n=(1/n²) sum_(i,j=0)^(n-1) delta_(ih,jh),
    nu=b1 delta_(-a1,0)+b2 delta_(0,-a2), a1,a2>0, b1,b2>=0.

For k,l=1,...,n define

    V1_k=max{1,(a1+(k-1)h)²}-k/n,
    V2_l=max{1,(a2+(l-1)h)²}-l/n,
    U_kl=max{1,max(a1+(k-1)h,a2+(l-1)h)²}-kl/n².

Then the complete family of closed square caps

    (sigma_n+nu)(Q)<=max{1,side(Q)²}                   (6)

is exactly equivalent to

    b1<=1, b2<=1,
    b1<=min_k V1_k, b2<=min_l V2_l,
    b1+b2<=min_(k,l) U_kl.                            (7)

For effective negative-axis positions 0<a1,a2<1, the first two atom caps are
already implied by k=1 and l=1. Thus the displayed n+n+n² prefix inequalities
suffice there, including unequal distances and unequal masses. No symmetrization is
used.

To prove necessity and sufficiency, enumerate actual tight source-bbox
closures, not subsets of a candidate receiver grid. For finite positive
sources, all square caps are equivalent to each bbox closure of span s
having mass at most max{1,s²}: a smallest containing square gives necessity,
and the captured set of any square is such a closure with span no greater
than that square's side.

An old-only closure has mass at most 1 and is automatic. An outside singleton
without old points requires b_i<=1. A closure containing q1 and old points but
not q2 has k prefix columns and l prefix rows, so its mass is b1+kl/n² and its
span is max{a1+(k-1)h,(l-1)h}. Since the second term is at most w<1, its cap
is max{1,(a1+(k-1)h)²}, independent of l. The strongest old mass at fixed k
is attained at l=n, exactly giving V1_k. The q2 case is symmetric without
interchanging its distance or mass with q1's.

A closure containing both outside points necessarily contains the origin.
Its old atoms are precisely k by l prefixes, with mass b1+b2+kl/n² and span
max{a1+(k-1)h,a2+(l-1)h}; these are exactly U_kl. This proves both directions
of (7), including every boundary tie. The pure outside pair cap is not
forgotten: it is weaker than the k=l=1 inequality, because any closed square
containing both outside points also captures the origin. If a_i may exceed
1, the individual atom caps b_i<=1 must still be included separately; this
is why they appear in the general statement (7).

This theorem qualifies the full input before any synthetic collapse. It does
not assert that replacing the original prefix masses by a single mass m_P
at the origin preserves true corner deletion. In particular an area computed
for that collapsed measure remains a proxy until compared with the actual
full-grid halo and the same residual.

The separate core proof `proofs/r077-qualified-grid-total-mass-proxy-obstruction.md`
owns the qualified S20 proxy obstruction; no actual S_n deletion area is
computed or inferred in this appendix.

## 4. Minimal actual check and next decision

`computation/round077/check_mixed_only_counterexample.py` reads
`computation/inputs/r077/mixed-only-two-outside-counterexample.json` and checks
the exact rational radical bounds, receiver separation inequalities and
strict excess 849/50000. The Python 3.15.0 Fraction run actually passed,
exit 0. Its input/script/version/timing receipt and output are in
`computation/results/round077-mixed-only/checks.json`.

This author did not repeat the root's finite floating probe or the compute
worker's complete directed-area union. The analytic full residual partition
above is the proof; the minimal check only certifies its stated scalar
inequalities. No new full-grid halo area, random scan, MATLAB run, subagent
or commit is part of this component. The next decision is to retain the full
prefix positions and their global constraints when addressing genuine joint-
active corner witnesses; an unrestricted mixed-only half-mass lemma is no
longer an available shortcut.
