# Eight-point perimeter obstruction to insertion-proxy bounds of 2

Date: 2026-11-15 (Asia/Shanghai). New research result in `pcw-sol-cover-counterexample`, checked by the covering core. No changes to frozen R34 files. No commit or push; root integrates.

## Object and exact finite result

For an ordering π of finite positive atoms, insert p_j with mass m_j after an old subset P. For each nonempty S⊆P let d(S,p_j) be the maximum coordinate span of the bounding box of S∪{p_j}. Define

    M_eff(P,p_j)=max{m(S): d(S,p_j)²<m(S)+m_j},

with value 0 when no old subset is eligible. Let Q_j be the open square centered at p_j with area M_eff(P,p_j)+m_j, and U_π=∪_j Q_j. All squares are centered, axis parallel and have a common radius in their two coordinates. The strict eligibility inequality is part of the definition.

This is an upper-bound proxy: the point-insertion shell inclusion, applied at each step, gives H(μ_final)⊂U_π. It is not an exact expression for the true halo.

Take a=49/100 and the eight perimeter points {−a,0,a}²\{(0,0)}, each of mass 1/8. Their total mass is 1. The exact minimum proxy area over all 8! orderings is

    min_π |U_π| = 18121/10000 + 3sqrt(6)/25 − sqrt(2)/200
                ≈ 2.098967701322116 > 2.                 (1)

Thus no uniform inequality min_π|U_π|≤2M can hold for this proxy. This does not give a C₂ lower bound, and does not exclude a larger uniform proxy upper bound, including one below the inherited general upper bound 3.590865312933.

## Exact eligibility reduction

All coordinate spans are 0, a, or 2a, and

    a²=2401/10000<1/4,
    (2a)²=2401/2500>7/8,   (2a)²<1.

Every box of maximum span at most a contains at most three of the eight perimeter points. Therefore before the last insertion the proxy square area is exactly one of 1/8, 1/4, 3/8: a singleton, an eligible neighboring pair, or an eligible three-point corner block. Every subset with four or more points has span 2a and fails eligibility before all eight points are present. The last insertion has area 1.

The production calculation uses exact Fraction comparisons and a subset dynamic program for M_eff. The independent verifier uses the preceding local classification directly: one neighbor supplies a pair, two mutually neighboring old points supply a triad, and only the final insertion supplies the global square. These two eligibility implementations agree on every ordering record.

## All-order finite certificate

Index vertices in increasing x and then y order:

| j | p_j | Square area times 8 in the minimizing pattern |
| --- | --- | --- |
| 0 | (−a,−a) | 1 |
| 1 | (−a,0) | 8 |
| 2 | (−a,a) | 1 |
| 3 | (0,−a) | 2 |
| 4 | (0,a) | 2 |
| 5 | (a,−a) | 1 |
| 6 | (a,0) | 3 |
| 7 | (a,a) | 1 |

D4 symmetry reduces the last vertex to a corner (index 0) or an edge midpoint (index 1). Every possible order of the remaining seven vertices is retained: 5040+5040=10080 representative orderings. Their area patterns reduce to 348 distinct signatures. The CSV ledger records the signature of every ordering, and the JSON stores every signature's eight area types, representative order, multiplicity and exact area bounds.

For each area type T=num/8, half-side r satisfies r²=num/32. With b=32, set

    j_lo=isqrt(num·2^(2b−5)),
    j_hi=j_lo if j_lo²=num·2^(2b−5), else j_lo+1.

Then j_lo/2^b≤r≤j_hi/2^b. With common coordinate denominator D=100·2^b, all centers and these inner/outer square endpoints are integers. The exact union area of the inner squares is a rigorous lower bound; that of the outer squares is a rigorous upper bound. Production uses an integer x-sweep and interval union. A separate verifier recomputes every area by integer rectangular-cell coverage and imports no production area helper.

The resulting global minimum interval is

    378116406113191688339/180143985094819840000
      ≤ min_π|U_π| ≤
    12099725004195137798823/5764607523034234880000.

In particular 2.0989677002<min_π|U_π|<2.0989677018. Every signature's inner area is strictly greater than 2D². The only signature whose interval can attain the minimum is number 237, with the area multipliers in the table. Its representative insertion order is

    (0,2,5,3,7,4,6,1).

All other signatures' lower bounds exceed signature 237's upper bound; the next lower bound is approximately 2.15760323434. Thus the finite interval certificate also identifies the exact minimizing area pattern, rather than merely finding a good order.

## Closed-form area of the minimizing pattern

Let r=sqrt(2)/8 be the singleton half-side and u=sqrt(6)/8 the triad half-side. The global square is centered at (−a,0), the triad at (a,0), pair squares at (0,±a), and singleton squares at all four corners. The strict inequalities

    a>r+u,   a−r>1/4,   a−u>1/2−a,   u>6/25

are certified using 64-bit dyadic root bounds in the independent verifier. The only positive-area pairwise intersections are: the global square with the two pair squares; the global square with its two same-side corner singletons; and the triad with the two pair squares. These six intersections are mutually disjoint, so there are no positive-area triple intersections.

The total square area is 19/8. The preceding intersection areas sum to

    2(13/50)² + 4r(r+1/100) + 2(u−6/25)².

Subtracting and using r²=1/32, u²=3/32 gives exactly (1). The verifier checks the rational scalar simplification, overlap inequalities and 64-bit interval for the formula against the all-order interval certificate.

## Actual runs and reproducibility

Python 3.15.0, standard library only, one process per run, no pseudorandom inputs and no floating tolerance in the certificate. Production: 1024 exact subset eligibility checks, 10080 orderings, 348 signatures; elapsed 0.063746 seconds, exit code 0. The separate verifier checks every ordering and all 696 inner/outer areas; its actual final run is recorded separately. Floating values are included only to display approximate areas.

Necessary files are `computation/eight_point_order_proxy.py`, `computation/verify_eight_point_order_proxy.py`, and the three result files in `computation/results/eight-point-order-proxy/`: `exact-eight-point-proxy.json`, `all-order-signatures.csv`, `independent-cell-verification.json`. A small manifest records hashes and scope. Run from the topic directory:

```sh
python3 -B computation/eight_point_order_proxy.py
python3 -B computation/verify_eight_point_order_proxy.py
```

No larger point count or additional random batch was run. The general behavior of this proxy above 2.0989677013, and any useful uniform bound below the inherited 3.590865312933, remain open.
