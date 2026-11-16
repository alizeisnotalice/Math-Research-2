# R72: exact receiver-cell enclosures for a Borel shared-source family

2026-11-16. This is a direct computation for atomic sigma plus an absolutely
continuous outside measure. It is not a finite-atomic approximation, subtraction
transfer, or a counterexample in a finite-source selection domain.

## Family and the separate global premise

For integer n>=2 put h=(n+1)/n², w=1-1/n², t=n/(n+1), and

    sigma_n = n^-2 sum_(i,j=0)^(n-1) delta_(ih,jh),
    J=[-1/n²,1]²,  B=[-1,2]²,
    nu_n=t 1_(B\J) dx,       mu_n=sigma_n+nu_n.

The global premise supplied by the root/Astra proof and independently reviewed
there is that every axis-parallel square Q satisfies

    mu_n(Q) <= max{1, area(Q)}.                           (1)

This computational appendix verifies receiver cells conditional on (1); it
does not certify (1) merely by checking finitely many squares. The separate
proof is Section 3 of
`proofs/r072-uniform-grid-borel-family-and-atomic-discontinuity.md`; its Section 5
also proves the same critical-event reduction and confirms that replacing B
by [-1,w+1]² gives the identical halo. The integrated round report binds the
final source version of that proof. Since nu has density
at most t<1, H(nu) is empty, directly at all radii. By (1), a strict full-halo
witness has common radius r<1/2. With no captured sigma atom, density below 1
cannot produce a strict witness. Hence the full halo is contained in

    C=[-1/2,w+1/2]².

For x in C and r<1/2, Q(x,r) is contained in B. Thus the closed-square mass
excess is positive exactly when

    sigma_n(Qclosed(x,r))
       > t area(Q(x,r) intersect J) + (1-t)4r².           (2)

The density has zero mass on square boundaries; sigma is counted on the closed
square. Strict (2) at a critical radius also yields a slightly larger open
square witness, because the right side is continuous in r and the already
captured atoms remain captured. This includes critical radius zero. For a
fixed point, sigma counts change only at the union of the n distances on each
axis. Between events, both right-side costs are nondecreasing, so checking
all closed critical events below 1/2 is complete. Ties must be taken together.

## Whole-cell classification

For a closed receiver cell D=[a,b]x[c,d] and one axis grid point p, define

    d_guaranteed(p,[a,b])=max{|p-a|,|p-b|},
    d_possible(p,[a,b])=max{a-p,p-b,0}.

At radius r, let Nx_g(r), Ny_g(r) count grid coordinates with guaranteed
distance at most r; define Nx_p, Ny_p using possible distance. Then at every
x in D,

    Nx_g(r) Ny_g(r) / n² <= sigma_n(Qclosed(x,r))
       <= Nx_p(r) Ny_p(r) / n².                          (3)

The possible atoms need not be simultaneously realizable at a common center;
the right side is deliberately an upper bound, so that coupling loss is safe
for exclusion. The guaranteed set, in contrast, is captured at every center.
Both products use the same square radius r in both axes.

For the one-dimensional J interval [j0,j1], the overlap length with
[x-r,x+r] is symmetric about (j0+j1)/2 and nonincreasing with distance from
that midpoint. Its maximum for x in [a,b] is attained by projecting the
J midpoint into [a,b], and its minimum is attained at the farther endpoint.
The two coordinates vary independently, so the exact minimum and maximum
overlap areas A_min(D,r), A_max(D,r) are the products of these axis extrema.
These functions are nondecreasing in r; in particular the minimum over D is
nondecreasing because every individual overlap area is nondecreasing.

An INNER certificate chooses any fixed r<1/2 for which

    Nx_g(r) Ny_g(r) / n² > t A_max(D,r)+(1-t)4r².         (4)

It proves the **entire closed cell** lies in the full halo. The algorithm tests
all guaranteed-distance events, at most 2n. A strict uniform margin in (4)
allows a slight radius expansion if open square witnesses are required.

An OUTER exclusion checks

    Nx_p(r) Ny_p(r) / n² <= t A_min(D,r)+(1-t)4r²         (5)

at every possible-distance event below 1/2, again at most 2n and retaining zero
and all ties. Between these events the left side is constant and the right
side nondecreasing, so the complete radius interval is covered. Before the
first event the left side is zero; after the last event the same monotonicity
covers the tail up to 1/2. Radii at least 1/2 are excluded by the separate
global premise (1). Therefore (5) proves every point in the closed cell is
outside the full halo at every radius. Equality is safe for exclusion.

Cells passing neither test are UNRESOLVED and charged only to the upper area.
The code never discards them because a sampled center is negative.

## Exact integer implementation and actual runs

`computation/round072/borel_grid_adaptive.py` uses an integer coordinate scale
D0=2n²*2^16, accommodating all grid coordinates, J midpoint, domain endpoints
and subdivisions through depth 16. If r and overlap area are measured in
integer units, (4) and (5) reduce to exact integer comparisons of

    count*(n+1)*D0²    versus    n³*overlap_units+4n²*r_units².

There are no floating-point classification decisions. A largest-area-first
quadtree partitions C into classified and unresolved leaves. Replacing a
parent by its four children preserves area; the final certificate saves every
leaf, its depth, closed integer box, state, and positive inner witness or
complete-event exclusion summary. Pending cells at a budget stop remain
explicitly unresolved. Shared leaf boundaries have area zero, while the
classification itself remains valid on the closed cells.

The actual runs used Python 3.15.0 standard library, two independent processes,
target unresolved area 1/50, depth cap 16 and node cap 200000:

| n | Certified area lower | Certified area upper | Processed cells | Target reached |
| --- | --- | --- | ---: | --- |
| 5 | 27886815467/10485760000 | 28096528411/10485760000 | 135462 | yes |
| 17 | 233950207087/87578116096 | 14988796509/5473632256 | 200000 | no |

The decimal intervals are about [2.659494,2.679494] and
[2.671332,2.738364], respectively. Both upper endpoints are strictly below
7/2 by exact rational comparison. The n=17 unresolved area is about 0.06703,
larger than the requested 0.02; this budget stop is retained explicitly. Since
the certified upper endpoint is already far below 3.5, no further refinement
was performed. These runs took about 1.69 and 4.33 seconds internally and
returned exit 0. The saved compressed leaf certificates are about 2.03 and
3.33 MB; both are complete finite certificates, not sampled receiver lists.

For either n, Delta=H(mu_n) because the residual halo is empty. The result
establishes only the displayed actual-family area interval conditional on the
separate cap proof. It gives no uniform result over all n, no finite-atomic
deletion transfer, and no general sharp weak constant.

## Files and reproduction

The result files are `computation/results/round072-borel/n5.json` and
`n17.json`, each pinning its compressed complete leaf certificate, script,
parameters, Python/platform and process id. To reproduce in a fresh output
directory on data:

    python3 -B computation/round072/borel_grid_adaptive.py --n 5 \
      --target-error 1/50 --max-depth 16 --max-nodes 200000 \
      --out-dir /Volumes/data/CodexTemp/T2/r072-borel-recheck

Replace 5 by 17 for the second case. A reader can independently reclassify
each saved leaf from its n, D0 and endpoints without reconstructing any search
history. Full partition verification additionally checks that these dyadic
leaves form the saved domain and have the stated total classified/unresolved
areas.

The independent review actually checked the entire n=5 partition (198139
leaves) and reclassified all 30669 exclusion leaves and all 38747 inner leaves
using literal two-dimensional atom distances instead of the production axis
products. Its 1735400 atom tests and 54705 actual two-dimensional possible
events reproduced both exact area endpoints. This was saved-certificate
verification, not another adaptive search; it took about 1.35 seconds, exit 0.
For n=17 it checked the full 306433-leaf partition and complete area ledger,
then independently reclassified 16 deterministic quantile leaves of each
classified state, 32 classified leaves in total. Those checks passed, taking
about 1.43 seconds, exit 0. It did **not** independently reclassify every n=17
leaf or rerun that adaptive search. The source-code and all-radius lemma were
also independently reviewed. The separate audit result/receipts pin those
different verification scopes; neither the n=17 budget target nor a full
n=17 classification replay is claimed.
