# R71: independent all-radius proof for five genuine deletion cells

2026-11-16. Threshold 1, Lebesgue area in R², centered axis-parallel squares
with one common radius in both coordinates. This appendix certifies a small
positive area contained in actual deletion for the specified 33-atom input.
It neither computes the full halo nor bounds deletion area from above.

## 1. A robust-cell lemma from complete closed prefixes

Let nu be a finite positive atomic residual measure. Fix a center x and list
the distinct positive L-infinity distances of its residual atoms as

    0 < d_1 < ... < d_m,
    M_j = nu({p : ||p-x||_infinity <= d_j}).

Every atom tied at d_j must be included in M_j. Suppose epsilon>0 satisfies

    epsilon < d_1,
    M_j <= 4(d_j-epsilon)^2 for every j=1,...,m.          (1)

Then **every** y in the closed square C=x+[-epsilon,epsilon]² is outside the
residual halo, with all residual radii retained:

    nu(y+[-r,r]²) <= 4r² for all r>=0.                   (2)

To prove this, fix y and r, and let T be the residual atoms in the closed
square y+[-r,r]². If T is empty there is nothing to prove. Otherwise let d_j
be the largest distance to x among atoms of T. The triangle inequality gives
d_j<=r+epsilon. The mass of T is at most the **complete** closed prefix M_j,
so

    nu(T) <= M_j <= 4(d_j-epsilon)² <= 4r².

The last comparison uses d_j>epsilon, hence d_j-epsilon>=0. This argument
uses one square radius r in both coordinates throughout; no independent
coordinate radii or finite radius menu is substituted. It includes r=0,
all inter-event radii and the whole unbounded tail. It works with non-strict
inequalities in (1), although the first five cells below have strict margins.
Checking just one residual radius would not establish (2).

Now let mu be the full measure and let U be its atoms in the **closed** square
x+[-r0,r0]², with mass m_U. If a fixed number R satisfies

    R > r0+epsilon,           m_U > 4R²,                 (3)

then for every y in C, every p in U satisfies

    ||p-y||_infinity <= r0+epsilon < R.

Thus the same open square radius R is a strict full-halo witness at every y
in C. Combining (2) and (3),

    C subset H(mu) \ H(nu) = actual deletion.            (4)

For several centers all measures are the same mu and nu; a single epsilon
and R may be used when (1) and (3) hold at each center. Their certificates
therefore represent simultaneously realized outside sources, rather than
separate source choices at each receiver. This lemma needs no maximal-block
assumption; admissibility of the specific chosen S is a separate audit.

## 2. An independent whole-cell residual check

For an arbitrary closed receiver rectangle C=[u,v]x[s,t] and nonempty
residual subset T with source bbox [a,b]x[c,d], put

    delta_C(T)=max{(b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,0}.

This is exactly min_(y in C) max_(p in T)||p-y||_infinity. Indeed, in the
first coordinate a radius r permits a center in C if and only if
[b-r,a+r] intersects [u,v]; this gives the first coordinate maximum of the
three displayed constraints. Do the same independently in the second
coordinate and then take their maximum. The minimum is attained on compact C.
Consequently

    C intersect H(nu) = empty
    iff nu(T) <= 4 delta_C(T)² for every nonempty T.      (5)

If the inequality fails, a minimizing center and a slightly larger radius
give a strict halo witness, including when the minimizing radius is zero.
Conversely a halo witness supplies a violating captured subset. Formula (5)
is the core proof's Section 9 criterion, independently checked here. For this
8-atom residual there are only 255 nonempty subsets, so the actual independent
run checked all of them for each of the six proposed cells. This second route
does not use the complete-prefix reduction or the LP matrix. In larger inputs
one may compress to complete source bounding-box closures; no such compression
was needed in this small check.

## 3. Actual input, exact parameters and positive cells

The actual input is the compute worker's
`computation/inputs/r071/positive33-actual-input.json`: 25 internal atoms at
(6i/25,6j/25), 0<=i,j<=4, each of mass 1/25, and eight positive outside atoms.
This independent appendix reads their exact positions and masses directly;
it does not import the LP or the production geometry program. The input SHA-256
is recorded in the result receipt. We use

    r0=1/4,  epsilon=101/40000,  R=253/1000.

In particular R-(r0+epsilon)=19/40000>0. The five certified centers are

    x0=(-1/10,6/25),   x1=(-1/10,12/25),
    x2=(-1/10,24/25),  x3=(0,53/50),  x4=(12/25,53/50).

Define the five **closed** squares C_k=x_k+[-epsilon,epsilon]². These are exact
containment claims, including the cell edges and corners. The independently
reconstructed complete residual prefix counts and their smallest displaced
slacks are:

| Center index | Complete prefix count | min_j [4(d_j-epsilon)²-M_j] |
| --- | ---: | --- |
| 0 | 5 | 30314201/400000000 |
| 1 | 6 | 38030601/400000000 |
| 2 | 7 | 18191001/400000000 |
| 3 | 7 | 18191001/400000000 |
| 4 | 6 | 29519001/400000000 |

All are strictly positive and all d_j exceed epsilon. At each center the
original closed radius 1/4 captures full mass exactly 2601/10000. Thus the
uniform open radius R has strict mass margin

    2601/10000 - 4(253/1000)² = 127/31250 > 0.

The actual worst distance of the captured full label over the whole cell is
9701/40000 for centers 0,2,3,4 and 9901/40000 for center 1. These stronger
direct geometric checks also lie strictly below R. Every full witness uses
the **same** full input and the **same common open square radius** R. By (4),

    C0 union C1 union C2 union C3 union C4 subset H(mu)\H(nu).

All ten pairs of these closed squares are disjoint (in fact separated by a
positive coordinate gap). The independent code also computed their union area
by an exact x-slab/y-interval union routine, rather than merely multiplying
the number of receiver points. Both calculations give

    area(union_(k=0)^4 C_k)
       =5(2epsilon)²=10201/80000000=0.0001275125.         (6)

This is a rigorous area lower certificate for actual deletion for this one
input. It is not a complete deletion area, and it is not automatically a lower
certificate for the smaller mixed region H(mu)\(H(S) union H(nu)): the cells
were not required to avoid the internal halo. No general C2 bound follows.

## 4. Why the sixth point is omitted

The sixth proposed center is x5=(18/25,53/50). Its first closed residual
prefix is a single outside atom at (18/25,6/5), of mass 49/625, at distance
7/50. There is exact equality

    49/625 = 4(7/50)².

For the proposed positive epsilon the displaced-prefix slack is
-1120999/400000000. The all-subset whole-cell check finds the same negative
slack. More directly, move the receiver a small positive distance upward toward
that atom; its distance decreases while its mass remains fixed, producing
a strict residual-halo witness. Therefore no positive symmetric neighborhood
of x5 can be residual-halo-free. Its full-halo positive margin cannot repair
this residual obstruction. No area was counted for x5.

## 5. Actual independent execution and reproducibility

The new independent files are

* `computation/round071/check_robust_deletion_cells.py`;
* `computation/inputs/r071-robust-cell-independent-config.json`;
* `computation/results/r071-robust-cell-independent.json`.

The standard-library Fraction program actually ran under Python 3.15.0 on
macOS arm64, exit 0, about 0.0372 seconds internally. It reconstructed all six
complete prefix tables (37 total events), all captured full labels, all five
box bounds and all ten separation tests. It independently checked 1530
residual subsets (255 at each of six cells), verified (5), and computed the
exact union area. The final comparison against the production certificate
confirmed every radius, mass and central slack in its six complete prefix
tables, the five bounds, common open radius and final area. The result receipt
pins the actual input, configuration, compared production certificate and
independent script by SHA-256; floating point appears only in the optional
display of the exact area. No optimization or random search was run here.

From the topic directory, after integrating the compute worker's input and
production certificate:

    python3 -B computation/round071/check_robust_deletion_cells.py \
      --input computation/inputs/r071/positive33-actual-input.json \
      --config computation/inputs/r071-robust-cell-independent-config.json \
      --claim computation/results/round071-joint/primal-geometry-certificate.json \
      --output /Volumes/data/CodexTemp/T2/r071-robust-cell-recheck.json

The `--claim` comparison is optional. Without it, the same mathematical checks
run solely from the actual atom input and declared cell configuration. This
independent worker did not rerun the production maximum-block/terminal tests,
the LP, or full halo area; their evidence belongs to their respective audits.
No commit or push is made from this delegated worktree.
