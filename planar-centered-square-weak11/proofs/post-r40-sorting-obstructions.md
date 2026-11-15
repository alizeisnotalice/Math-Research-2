# After R40: projection-order obstruction and an interleaving caution

Date: 2026-11-15 (Asia/Shanghai). Independent adversarial verification requested by the covering core. New files only; R34 and R40 are unchanged. No commit or push. This note does not refute the conjectural free-order bound min_π|U_π|≤3M.

## Every monotone projection order can be bad while free ordering is almost optimal

Use the same effective-mass insertion proxy as in `point-insertion-shell-budget.md` and the R40 proof: a new point of mass m receives a centered square of area M_eff+m, where M_eff is the largest eligible old mass and eligibility means strict bounding-box maximum span squared less than combined mass.

Let one heavy point of mass 1 be at the origin. Put eight light points, each of mass ε=1/10000, at

    {−a,0,a}²\{(0,0)},   a=99/100.

Total mass is M=1+8ε=1251/1250. Consider increasing order of any nonzero linear functional ℓ(x,y), with arbitrary tie breaking.

The eight light points form four antipodal pairs in four distinct unoriented directions: horizontal, vertical, and the two diagonals. A nonzero ℓ has a one-dimensional kernel, which can contain at most one of these pairs. Every other pair contributes exactly one strictly positive projection. Thus at least three light points have ℓ(p)>0 and must follow the heavy point, regardless of ties.

For every such light point q, the old singleton consisting of the heavy point is eligible:

    ||q||∞²=a²=9801/10000<1+ε.

Its proxy therefore has area at least 1, as does the heavy point's own proxy. The union contains four unit squares centered at the heavy point and three distinct later light points. Every pair of these centers has maximum-coordinate distance at least a; therefore the intersection area of any two unit squares is at most 1−a=1/100. The elementary first Bonferroni inequality gives

    |U_π|≥4−6/100=197/50,
    |U_π|/M≥4925/1251≈3.93685051958.                    (1)

This conclusion holds simultaneously for every nonzero ℓ, including a direction selected after inspecting the masses. Thus changing the direction alone cannot repair a sorting rule constrained to a single monotone projection.

In sharp contrast, the free-order optimum for this exact input is

    min_π |U_π|=1+11ε=10011/10000,
    min_π |U_π|/M=3337/3336≈1.00029976019.               (2)

To prove this, insert all eight light points before the heavy point. Any group of at least two light points has span at least a, whereas their total mass is at most 8ε<a². Hence all light proxy squares are singletons, each of area ε, and are pairwise disjoint. At the final heavy insertion a box of span at most a contains at most three light points, and a three-light corner block exists. Any four light points have span 2a, with

    (2a)²>1+8ε.

Thus the heavy effective old mass is exactly 3ε, and its proxy area is 1+3ε. This proxy is disjoint from every tiny light square, since

    sqrt(1+3ε)+sqrt(ε)<2a.

So the union area is exactly 1+11ε. If the heavy point is not last, some later light point and the heavy point supply two unit squares at maximum-coordinate center distance a. Their union has area at least 2−1/100=199/100, strictly larger than 10011/10000. This proves (2) for all orders without a permutation search.

A family version makes the separation stronger. Replace a by 1−δ and let δ,ε decrease to zero, preserving 8ε<a², 1+8ε<4a² and sqrt(1+3ε)+sqrt(ε)<2a. Every projection order then has ratio at least (4−6δ)/(1+8ε), tending to 4, while the free-order optimum ratio tends to 1. This is an obstruction to this sorting class, not a lower bound on the free-order proxy and not a C₂ lower bound.

## Increasing individual mass also fails without cluster aggregation

Replace the heavy source by 100 distinct center points

    q_k=(k/100000,0),   k=0,…,99,   each of mass 1/100,

and put eight perimeter points at {−a,0,a}²\{0}, a=99/100, each of mass 1/50. Total mass is 29/25. Any order of nondecreasing individual masses must insert every center point first.

Regardless of ties within the center group, its final insertion has proxy area 1: the full center diameter is at most 99/100000, so all center mass is eligible. Every later perimeter point can use the entire center mass, because its common-capture span is at most 99099/100000 and this span squared is less than 1+1/50. Thus the proxy union contains eight unit squares centered at the perimeter points, plus a unit square centered at some q_k.

The eight perimeter unit squares leave only the central hole (−49/100,49/100)² inside the large square [−149/100,149/100]². Since the center displacement is less than 1/100, the center unit square contains this entire hole. Therefore every nondecreasing-individual-mass order has

    |U_π|≥(149/50)²=22201/2500,
    |U_π|/M≥22201/2900≈7.65551724138.                  (3)

There is again a simple good free order. Insert the eight perimeter sources first; their proxies are singleton squares, total area 4/25. Insert the center group afterward. Any eligible set involving a new center point can contain at most three perimeter points: four perimeter points have span 2a, whose square exceeds the total input mass 29/25. Thus every center proxy has area at most 1+3/50=53/50. These proxy centers have x-span below 1/1000, so their entire union lies in a rectangle of area below 53/50+(1/1000)sqrt(53/50). Since sqrt(53/50)<103/100, this order gives

    |U_π|<122103/100000,
    |U_π|/M<122103/116000≈1.05261206897.

Thus neither a single geometric projection nor increasing individual masses supplies a useful universal bound for this proxy. This does not refute increasing aggregate cluster mass or another mass-aware ordering.

## A minimal corner-cluster interleaving caution

Fix a candidate final point p=(0,0) of mass 1. Its old sources q=(1,0), mass 1, and r=(4,0), mass 14, lie in one closed coordinate quadrant relative to p. At p the only maximal eligible old cluster is C={q}: q is eligible because 1<2; r is not because 16≥15; the union q+r is not because 16=16, and strict eligibility is required.

Yet in the order r→q→p, q captures the effective old mass 14 because the q,r span squared is 9<15. Its proxy has area 15. The final p proxy has area 2 and is contained in q's proxy: sqrt(15)>2+sqrt(2), as follows from 9>4sqrt(2). Thus the proxy belonging to C's own point has area 13 outside the final p square, much larger than 2m(C)=2.

This refutes any claim that C's proxy cost outside the final square is bounded solely by 2m(C) under arbitrary interleaving with external sources. It does not refute the core's existence-of-an-order total budget: one may need C-before-external ordering or joint payment for outside mass that enlarges C's proxy.

## A valid restricted separator interface

The structure core's eligible-hypergraph reduction survives this audit. After removing a proposed last point p, form a hypergraph on old atoms with edges B satisfying bbox-span(B)²<m(B). Every eligible group used by a prefix proxy is a full old-input hyperedge and therefore lies within one hypergraph connected component. Old proxies from different components can be constructed independently; interleaving other components cannot enlarge their effective masses.

For any component of at most three atoms, mass-increasing insertion order has total proxy square area at most twice its mass: for a≤b≤c the sum bound is a+(a+b)+(a+b+c)=3a+2b+c≤2(a+b+c). Hence if some deleted p leaves only such components, inserting p last gives

    |U_π|≤M+2(M−m_p)=3M−2m_p.

For components of at most four atoms the same calculation gives the factor 5/2 and bound (7/2)M−(5/2)m_p. These are conditional interfaces; there is no proof here that every input admits such a separator.

## A minimal obstruction to a uniformly small component separator

Let a=49/100 and put equal mass 1/6 at the six points, in cyclic order,

    (a,0), (a/2,a), (−a/2,a), (−a,0), (−a/2,−a), (a/2,−a).

For the old-input eligibility hypergraph used in the preceding separator interface, its eligible-pair graph is exactly the six-cycle. Adjacent pairs have maximum-coordinate span a, with a²<1/3. Every nonadjacent pair has span at least 3a/2, whose square exceeds 1/3. Hence deleting any point leaves a connected five-point path, and therefore one hypergraph component containing every remaining old point.

Moreover every such five-point component has bounding-box maximum span 2a, and

    (2a)²=2401/2500>5/6.

Thus that connected component is not itself an eligible set. The full six-point input is eligible, since 2401/2500<1. This explicitly separates connected-component closure from saturation and shows that a separator with components of size at most four need not exist. Six is the minimal total source count for this failure by cardinality. It does not contradict the conditional separator theorem or prove a free-order lower bound above 3.

## A valid rank interface covers the six-cycle despite its giant components

The covering core's stronger all-source rank observation is valid. If every eligible hyperedge of the final input has cardinality at most k, mass-nondecreasing insertion gives each proxy area at most k times the mass of its newly inserted atom: every old atom in a winning eligible set has mass at most that new mass. Summing individual square areas gives |U_π|≤kM, regardless of component sizes.

If only the full input may have cardinality above two among eligible edges, insert a maximum-mass atom last. Every earlier proxy has area at most twice its new mass, and the last proxy area is at most M, yielding

    |U_π|≤3M−2m_max.

The six-cycle input above satisfies this latter condition. Every three-point set has span at least 3a/2, whose square is greater than 1/2, and every four- or five-point set has span 2a, whose square exceeds 5/6. Thus all proper eligible edges are only singletons and pairs; the full six-point set is eligible. It follows that this input has an order with proxy area at most 3−2/6=8/3. Its giant components therefore obstruct only the small-component separator criterion, not this stronger mass-rank criterion or the general target of 3.

## Actual evidence and remaining target

`computation/post_r40_sorting_obstruction_check.py` uses Fraction arithmetic only. It checks the exact numerical constants, all 16 projection sign sectors/boundary rays, all 70 four-light-point bounding boxes, the 100-center mass-sorting inequalities, the six-cycle graph, all six deletion components and all 41 proper subsets of size at least three, scalar gap conditions, and the three-point interleaving identities. The antipodal-pair argument proves the all-direction statement; the finite sector check is supporting arithmetic, not a replacement for that proof. No random inputs, free-order search, or additional point-family search was run.

The R40 proxy obstruction above 2 remains frozen. A general free-order bound of 3, an obstruction above 3, and a useful mass-aware global ordering or saturated-component allocation remain open. Mechanical nesting of the R40 eight-point pattern has not supplied a counterexample: preserving each child pattern's relative diameter consumes order sqrt(child mass) spatial span, while choosing an inner child point last can make the global square cover much of the proposed added proxy area. A valid recursive obstruction must control that free last-point choice rather than prescribe a bad order.
