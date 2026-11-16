# R75: a three-atom counterexample to the uncapped half-mass collapse bound

Let the background be Lebesgue measure on R2 and let H be the strict
threshold1 halo for centered axis-parallel common-radius squares. Put

    sigma=delta_(0,0),
    nu=4 delta_(-2,0)+4 delta_(0,-2),
    O={x<0,y<0}, D=H(sigma+nu) minus H(nu).

Then the complete actual deletion area is

    |D intersect O|=4sqrt(2)+sqrt(5)/2-25/4
                  approximately .524888238242275 >1/2.

This disproves only the proposed **uncapped** collapse auxiliary bound
with one internal atom and two outside atoms. The full input has eligible
mass9, and sigma is not a maximum block. The closed square[-2,0]^2 has
mass9 and area4, violating the normalized maximum-block square cap.
It is not a uniform-grid-cap counterexample, a violation of3.5, or a C2 bound.

## Complete full and residual geometry

For every nonempty source subset B of mass m and bbox[a,b]x[c,d], its
receiver rectangle is

    (b-sqrt(m)/2,a+sqrt(m)/2)
      x(d-sqrt(m)/2,c+sqrt(m)/2).

The complete finite halo is the union of these rectangles. Every actual
strict witness supplies one captured subset, and every point of such a
receiver rectangle has a slightly larger-than-capture radius strictly
below sqrt(m)/2. Subsets need not be bbox closed to be included: redundant
subsets only contribute valid receiver rectangles. Here all7 nonempty
full subsets and all3 nonempty residual subsets are retained. In particular,
neither the full triple nor the residual outside pair is omitted.

Write u=sqrt(5)/2 and v=sqrt(2), so1<u<v<3/2. Within O, the full rectangles
not already wholly residual are

    A=(-1/2,0)^2,
    E_x=(-u,u-2)x(-u,0),
    E_y=(-u,0)x(-u,u-2),
    T=(-3/2,-1/2)^2.

The actual residual is the union of

    L_x=(-3,-1)x(-1,0),
    L_y=(-1,0)x(-3,-1),
    R=(-v,v-2)^2.

A is disjoint from all three residual rectangles and contributes1/4.
The squares R and T satisfy R subset T. Both T intersect L_x and
T intersect L_y have area1/4 and have disjoint interiors. Their respective
intersections with R each have area(v-1)^2. Since |R|=4(v-1)^2,

    |T intersect H(nu)|=1/2+2(v-1)^2=13/2-4v,
    |T minus H(nu)|=4v-11/2.

E_x lies within T in its x coordinate and extends above T only over
y in(-1/2,0). On that extension, only L_x removes area, leaving precisely

    (-1,u-2)x(-1/2,0),

of area(u-1)/2. E_y contributes the transpose with the same area.
These two additions, A and the deletion inside T have pairwise disjoint
interiors. Summing gives

    1/4+(4v-11/2)+(u-1)=4v+u-25/4.

All original receiver rectangles are open, as required by the strict halo.
Changing finitely many rectangle boundary segments in area calculations
has no Lebesgue effect; no threshold equality is promoted to a witness.

## A small rational strict certificate

The rational comparisons

    2-(707/500)^2=151/250000>0,
    5/4-(559/500)^2=19/250000>0

give v>707/500 and u>559/500. Hence

    |D intersect O|>131/250=.524>1/2,

with strict excess greater than3/125. This certificate does not depend
on a displayed floating decimal.

## Actual directed computation and scope

The saved exact input is `inputs/r075/three-atom-collapse-d2.json`.
Python3.13.18 actually enumerated all7 full and3 residual source subsets.
For each root radius, integer isqrt gives lower/upper bounds in units
2^-72. The inner-full minus outer-residual area is a lower bound; outer-full
minus inner-residual is an upper bound. Rectangles are clipped only at the
orthant's two upper boundaries0. There is no artificial lower receiver cutoff.

Two separate exact integer area paths agree on both endpoint numerators:
vertical x-events with active y-interval merging, and a compressed2D
Cartesian arrangement using exact midpoint coverage predicates. The saved
interval is approximately[.524888238242275,.524888238242275]; its rational
endpoints and width are retained in `d2-counterexample.json`. A separately
formed root interval for the hand formula intersects this enclosure.
These are two algorithms in this production code; independent reviewer
reconstruction is owned and recorded separately by audit/core.

The first directed candidate d=2 already gives a strict violation.
No random batch was launched, and d=4,8,16 were not evaluated. This is
one input, not an optimization or a claim of a sharp uncapped constant.
The mathematical bound with only one outside atom remains unaffected.

The process took.017860 seconds internally, exited0, and PID5411 was
checked absent. Actual source/input/result hashes and the execution
command are retained. No additional agent, MATLAB, commit, or push was used.
All R74 main/preliminary/supplement and R73 frozen files remain unchanged.

From the topic root, after creating the output directory, reproduce with:

    python computation/round075/certify_three_atom_collapse.py --input computation/inputs/r075/three-atom-collapse-d2.json --out computation/results/round075-collapse/d2-counterexample.json --bits 72
