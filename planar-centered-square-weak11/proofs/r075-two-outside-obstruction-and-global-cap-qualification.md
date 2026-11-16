# R75: two outside atoms break the unrestricted half-mass payment

2026-11-16. Strict threshold 1, centered axis-parallel common-radius squares.
This note leaves every R74 file unchanged. It supplies an exact obstruction
to extending the **no-cap** R74 collapse argument to two outside positions,
and identifies specific global grid constraints missing from that relaxation.
It does not give a counterexample to the original uniform-grid corner budget.

## 1. A three-atom family with its complete residual halo

Put sigma=delta_(0,0). For d>=2 define

    nu=d² delta_(-d,0)+d² delta_(0,-d),
    O=(-infinity,0)²,
    a=sqrt(d²+1)/2, b=d/sqrt(2), c=sqrt(2d²+1)/2,
    v=d/2.

Let D=H(sigma+nu)\H(nu). Every nonempty subset of these three sources is
strictly eligible. Thus the full halo is the union of all seven receiver
rectangles, and the residual halo is the union of the two residual singleton
squares **and their pair receiver rectangle**. No residual radius or source
label has been dropped.

It is convenient to use positive receiver coordinates X=-x, Y=-y on O. The
full rectangles that can contribute after residual subtraction are

    A=(0,1/2)²,
    L=(d-a,a) x (0,a),     L^T,
    T=(d-c,c)².

Here A is the internal singleton, L and L^T are the two internal/outside pair
receivers, and T is the triple receiver. The true residual in this orthant is

    R_L=(v,3v) x (0,v),     R_L^T,
    R_pair=(d-b,b)².                                    (1)

These formulas follow from intersecting the appropriate source-centered
open squares of radius sqrt(label mass)/2. They account for all strict
centered-square witnesses, not just a chosen radius. Receiver boundaries
have zero Lebesgue area; source capture remains closed/open as appropriate
to the usual slight-radius-expansion argument.

For d>=2 the inequalities v<a<=b<c<d, d-a>=1/2 and d-c>=1/2 hold. Hence A
is disjoint in its interior from T, L and L^T, and also from the residual.
The part of L outside T and outside the residual is exactly the stem

    (d-a,v) x (0,d-c),                                 (2)

of area (a-v)(d-c); the other stem is its transpose. These stems and A have
disjoint interiors. Inside T, inclusion-exclusion with the complete residual
(1) gives

    |T\H(nu)|=2[(c-v)²-(b-v)²]
              =1/2-2d(c-b).                            (3)

For example, each residual singleton cuts a square of side c-v from T, while
its intersection with R_pair is a square of side b-v. The two singleton cuts
meet only on a boundary. The residual pair rectangle itself has area
4(b-v)² and is contained in T. These are the complete terms in (3).

Therefore the exact whole-orthant deletion formula is

    |D intersect O|=1/4+2[(c-v)²-(b-v)²]
                         +2(a-v)(d-c)
       =3/4-2d(c-b)+2(a-d/2)(d-c).                     (4)

The third mixed label T is essential. Summing only two separate two-atom
payments would miss this coupled source label, while subtracting only the
residual singletons would overstate the deletion by ignoring R_pair.

## 2. The exact d=2 counterexample

At d=2, a=sqrt(5)/2, b=sqrt(2), c=3/2. The three disjoint contribution types
in (4) are

    internal singleton:            1/4,
    two exterior stems together:   a-1,
    triple after true residual:    4b-11/2.

Thus

    |D intersect O|=4sqrt(2)+sqrt(5)/2-25/4
                   approximately 0.524888238242275.     (5)

It is strictly greater than the proposed m/2=1/2. A short rational certificate
is sqrt(2)>24/17 and sqrt(5)>20/9, verified by squaring positive rationals.
Substitution in (5) gives

    |D intersect O|>311/612=1/2+5/612.                 (6)

Omitting R_pair from the subtraction in T would incorrectly add
6-4sqrt(2) to (5). Formula (5) does include it. The root/compute component
has separately checked d=2 by all seven full and all three residual labels
with two exact directed-area algorithms; this note does not rerun that union
or attribute its execution to this worker.

For context, rationalizing a-v and c-b in (4) gives

    2(a-v)(d-c)
      =[d-sqrt(2d²+1)/2]/[sqrt(d²+1)+d],
    2d(c-b)=d/[sqrt(2d²+1)+sqrt(2)d].

Consequently

    lim_(d->infinity) |D intersect O|=5/4-1/sqrt(2)>1/2.

This is an analytic limit of (4), not a scan of d. No new d=4 or d=8 numerical
area evaluations were run. It proves that an unrestricted northeast internal
mass collapsed to one atom cannot retain the R74 m/2 payment with two outside
atoms. It says nothing by itself about the original full-grid square caps.

## 3. A scaled S3 embedding passes every unit cap but fails a nonunit cap

Scale the d=2 construction by 1/3 in space and 1/9 in mass. Its local internal
atom has mass 1/9 at the origin, and its outside measure is

    beta=(4/9)delta_(-2/3,0)+(4/9)delta_(0,-2/3).        (7)

Now restore the other eight internal atoms to form the **complete** S3 grid,
with coordinates {0,4/9,8/9}² and masses 1/9 each. Every closed unit square Q
satisfies

    (sigma_3+beta)(Q)<=1.                               (8)

Here is an all-position proof, including ties. If Q contains no outside atom,
its mass is at most the internal total 1. If it contains only the left outside
atom, its first coordinate interval contains -2/3 and has upper endpoint at
most 1/3. It captures at most the internal column x=0, so internal mass is at
most 1/3 and full mass at most 7/9. The bottom outside atom is symmetric. If
Q contains both outside atoms, each coordinate interval contains -2/3 and 0
and has upper endpoint at most 1/3; the only internal point captured is the
origin. Its mass is exactly 8/9+1/9=1. Since 1/3<4/9, no boundary convention
can admit an additional internal row or column in this argument.

By inclusion in a unit square, (8) also gives the required small-square cap
for every side R<=1. Nevertheless the global square cap fails at

    Q_*=[-2/3,4/9]²,  side(Q_*)=10/9.

It captures four old grid atoms and both outside atoms, so

    (sigma_3+beta)(Q_*)=4/3,
    side(Q_*)²=100/81,
    mass-side²=8/81>0.                                 (9)

The captured six-atom label is strictly eligible with mass greater than 1.
Thus S3 is not a maximum eligible mass-1 block in this embedding. This is not
a qualified grid counterexample, despite all unit and smaller square caps
passing. In particular, if the two outside masses in (7) are replaced by a
common mass M, this one necessary global cap already requires M<=32/81,
whereas (7) uses 4/9=36/81.

The obstruction has a precise geometric location: the SW effective source
window for S3 has upper coordinate 1/2-1/9=7/18, while Q_* reaches 4/9>7/18.
The extra internal atoms involved in (9) cannot appear in a valid r<1/2 SW
corner witness, but their positions still impose global capacity constraints
on the same outside mass vector. A local model cannot discard them from the
caps merely because they are absent from its small-radius witnesses.

## 4. Whole orthant versus the actual corner: a second distinction

The whole-orthant value (5) is also not a value for the finite corner region.
For the **scaled local three-atom measure only** (the origin atom and (7),
without restoring the other eight S3 atoms), consider

    Gamma_3=(-1/2,-1/9)².

In the original unscaled positive coordinates this becomes
[1/3,3/2]², followed by area division by 9. The internal singleton contribution
there is 1/36, the two stems together contribute (a-1)/3, and the triple
contribution is unchanged at 4b-11/2. Hence the scaled local three-atom deletion
in Gamma_3 is exactly

    [1/36+(a-1)/3+4b-11/2]/9
       =4sqrt(2)/9+sqrt(5)/54-209/324.                  (10)

Using sqrt(2)<10/7 and sqrt(5)<9/4 gives the strict upper bound

    (10)<143/4536<1/18=m/2.                            (11)

Thus clipping even this relaxed local construction does not retain its
whole-orthant violation. Equation (10) is **not** the corner deletion of the
full S3 embedding: that larger input fails (9), so its witnesses cannot be
restricted to r<1/2 using the original maximum-block caps. Additional internal
sources may then influence its full halo. The two distinctions are separate:
correct corner clipping and retaining all global square caps are both needed.

The scaled local three-atom measure itself has total mass 1, so it automatically
satisfies every cap mass(Q)<=max(1,side(Q)²). This is still not the target grid
model: its selected internal measure has only mass 1/9. Its complete three-atom
label has mass 1 and span 2/3 and is eligible, so that local internal atom is
not a maximum eligible block. Restoring the required full internal S3 is
precisely what produces the nonunit capacity violation (9).

## 5. Actual small exact certificate and remaining direction

The minimal input is
`computation/inputs/r075/scaled-counterexample-unit-qualification.json` and
the checker is `computation/round075/check_scaled_unit_qualification.py`.
It partitions all unit-square center positions at every source entry/exit
coordinate, retaining singleton and open interval axis atoms. The 225 product
cells have constant captured source sets, so one exact rational representative
per cell proves (8) for **all** centers. Outside the saved coordinate hull,
one coordinate interval contains no source and mass is zero. The result saves
every cell, representative, captured indices and exact cap slack. This is a
complete finite certificate, not a finite sample offered as a general cap test.

The same Fraction run checks (9), the rational square-root bounds in (6) and
(11), and the positive excesses. It does not recompute any complete halo union:
(4), (5) and (10) are the analytic rectangle decompositions above. Actual
Python 3.15.0 execution passed, exit 0; input/script/version/timing receipt and
full 225-cell table are in
`computation/results/round075-unit-qualification/qualification-checks.json`.

From the topic directory:

    python3 -B computation/round075/check_scaled_unit_qualification.py \
      --input computation/inputs/r075/scaled-counterexample-unit-qualification.json \
      --output /Volumes/data/CodexTemp/T2/r075-unit-qualification-recheck.json

The unrestricted two-effective-outside collapse extension is therefore
strictly false. The natural remaining direction keeps the original internal
grid positions in the global nonunit caps while retaining the complete
residual union in a genuine corner calculation. No finite or Borel transfer
is used here, and no violation of the qualified uniform-grid corner budget
or its general 3.5 deletion target has been proved. This component
runs no MATLAB, spawns no agent, and makes no commit or push.
