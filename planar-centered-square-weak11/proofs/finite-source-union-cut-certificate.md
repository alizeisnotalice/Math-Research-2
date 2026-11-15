# R61: a positive, convergent certificate for finite source-union allocation

2026-11-15 (Asia/Shanghai). This theorem gives a finite, checkable, nonproportional transport certificate for each prescribed finite atomic input at threshold 1. It preserves all common-radius witnesses and all overlaps. It does not give a new uniform C2 bound. Unlike the disproved R56 proportional rule, the source allocation is free within its actual eligible source set.

## 1. Complete geometry and the actual eligible sources

Let mu=sum_{j=0}^{N-1} m_j delta_{p_j}, with distinct p_j in R² and m_j>0. Merge coincident points first. All squares are closed, centered at their receiver, axis-parallel and use one common radius. Let H={x: exists r>0, mu(Q(x,r))>4r²}.

For each tight coordinate bbox B of a nonempty subset of the source positions, take **all** sources lying in B, denoted I_B, their total mass M_B, and s_B=sqrt(M_B)/2. Deduplicate equal source closures. With B=[a,b]x[c,d], set

    E_B=(b-s_B,a+s_B)x(d-s_B,c+s_B).

Discard empty E_B, retaining degenerate source boxes. The complete identity H=union_B E_B holds. Indeed every strict square's full captured set has such a tight bbox and the same captured mass; conversely x in E_B admits a radius between its maximum distance to B and s_B, giving a strict witness which may also capture extra sources.

Define the actual eligible-source set

    N(x)={j: some actual strict Q(x,r) contains p_j}.

Then

    N(x)=union_{B:x in E_B} I_B.                          (1)

The reverse implication in (1) uses the same strict witness as above. The forward implication uses the tight bbox of all sources in that actual square. Thus the rule does not discard additional captured sources, replace a true bbox mass by a single atom, or combine unequal coordinate radii.

All quantities are finite. In particular a strict witness has r<sqrt(mu(total))/2. The finite rectangle-endpoint arrangement partitions H, modulo a null set, into finitely many positive-area cells on which N(x) is constant and nonempty.

## 2. The exact finite Hall cost

For S subset {0,...,N-1}, define m(S)=sum_{j in S}m_j and

    D(S)=|{x in H:N(x) subset S}|,    D(empty)=0.

A free source allocation may split a receiver's demand among its eligible sources. Its cost is the maximum source load divided by m_j. Every used source is realized by some actual strict witness at that receiver; different portions of the demand may choose different witnesses.

The minimum cost is exactly

    K_exact(mu)=max_{nonempty S} D(S)/m(S).               (2)

Necessity follows by adding the demand whose only eligible sources lie in S. For sufficiency, use the finite true arrangement cells as demand nodes and source atoms as supply nodes. Every subset of demand cells has source neighborhood S, and its total area is at most D(S). Hence all finite Hall inequalities hold with capacities K_exact m_j. Finite capacitated-flow duality, or the usual maximum-flow/minimum-cut proof on the compact finite flow polytope, supplies an allocation. This argument concerns a finite graph with real areas, not an unproved infinite transport theorem.

The allocation is legal in the original geometry. If j in N(x), choose a rectangle B with j in I_B and x in E_B, resolving ties in a fixed finite order. The radius

    r(x)=(max_{p in B}||x-p||_infinity+s_B)/2

is positive and strictly below s_B, includes p_j, and is a measurable actual strict witness. If the maximum distance is zero the displayed radius remains positive. The cell/source flow can be spread uniformly on its cell, followed by this finite measurable witness rule. No strict witness radius supremum needs to be attained.

Formula (2) alone is an exact reformulation. The following rational construction is the positive computable certificate and supplies an explicit precision error.

## 3. Rational source-union upper cuts

Choose rational inner/outer radii s_B^-<s_B<s_B^+ (allow s_B^-=0 with empty inner rectangle). Build E_B^- subset E_B subset E_B^+ using the same bbox, and let H^+=union E_B^+. On H^+, label

    T^-(x)=union_{B:x in E_B^-} I_B.

This label can be empty on the uncertain fringe. Use the rational endpoint arrangement, and let

    W_T=|{x in H^+:T^-(x)=T}|.

All W_T are exact nonnegative rational rectangle-union areas. Source overlaps are combined by set union before any capacity is charged. Put

    Z(S)=sum_{T subset S} W_T.                            (3)

For every nonempty S,

    D(S)<=Z(S).                                          (4)

Indeed T^-(x) subset N(x), and H subset H^+. Thus every receiver in D(S) is counted in (3). It is essential that the uncertain W_empty is included in every nonempty upper cut. One must **not** demand the artificial inequality Z(empty)<=0: the true empty-source cut is D(empty)=0. One must also not interpret outer-only rectangle incidences as actual strict witnesses.

Consequently the rational number

    K_raw=max_{nonempty S} Z(S)/m(S)                     (5)

certifies a legal allocation for the full true H by (2). The complete family of 2^N cuts is evaluated by the subset zeta transform in O(N 2^N) arithmetic operations after the geometric cell weights have been computed. The geometric preprocessing has its own finite cost and is not included in this bound.

There is also an explicitly stored rational flow. For T nonempty, the W_T cells are inside H^-:=union E_B^- and all j in T are legal. Construct the finite graph with demand W_T and capacity K_raw m_j. Its Hall cuts hold because sum_{nonempty T subset S}W_T<=Z(S). The implementation scales all data to integers and records an exact Dinic flow. This is an actual nonproportional strict-witness allocation on H^-; it is not falsely presented as an explicit allocation of the outer fringe. Full-H existence follows separately from the certified cuts (4)–(5).

## 4. A quantitative finite-precision error

Define the computable nonnegative area

    Gamma=sum_B (|E_B^+|-|E_B^-|).

Outside union_B(E_B^+\E_B^-), the true, inner and outer rectangle memberships agree. Therefore for every nonempty S,

    0<=Z(S)-D(S)<=Gamma.                                 (6)

More explicitly, if a point counted by Z(S) is not in the true halo, some outer rectangle but no corresponding true rectangle accounts for it. If it is in H but not D(S), a true eligible rectangle contains a source outside S, whereas that rectangle's inner approximation cannot contain the point. Either way the point belongs to one of the displayed uncertainty bands; the sum of their areas bounds the union.

Writing m_min=min_j m_j, (2), (5) and (6) give

    K_exact <= K_raw <= K_exact+Gamma/m_min.              (7)

For fixed finite input, shrinking all endpoint gaps to zero makes Gamma tend to zero. Thus these positive rational certificates converge to the optimal free-allocation cost. This is not a uniform error estimate in N, geometry, or masses, and it gives no uniform C2 constant.

## 5. R53 supplies a genuine geometric source-union envelope

Let d(S) be the largest coordinate span of the source set S. If x is counted by D(S), a strict witness at x captures only sources in S, since every source of every actual witness is in N(x). Hence

    {x in H:N(x) subset S} subset H(mu restricted to S).

The independently proved R53 theorem gives, with

    B=1994727599733432151936139/604462909807314587353088 <33/10,

    D(S)<=B max(m(S),d(S)^2).                            (8)

The theorem is for the actual centered common-radius all-scale halo of any finite positive Borel measure supported in a square of side sqrt(M). The arbitrary-bbox form follows by adding mass at an existing point until the total is max(M,d²); halo monotonicity supplies the inequality. No source is removed from its overlap neighborhood in this argument.

Combining the independent bounds (4) and (8) yields the positive source-union certificate

    U(S)=min{Z(S),B max(m(S),d(S)^2)},
    K_geom=max_{nonempty S} U(S)/m(S),
    K_exact<=K_geom<=K_raw<=K_exact+Gamma/m_min.           (9)

In particular all cuts with d(S)^2<=m(S) are already paid at capacity B. To certify capacity B for an input, only its noncompact source cuts require the additional accurate Z(S) check. This is a concrete geometric pruning rule and a quantitative source-union compensation interface. It does not prove that every noncompact cut is paid at B, or provide a general mechanism to establish that assertion without further geometry.

R53 source pins: author `proofs/r053-nine-cell-mass-simplex-certificate-33over10.md`, SHA256 `f7e9223403939788d940ebc1cfd47d69ad772eb44da73e42cdeb08d6f80fd784`; independent Borel-measure proof and mass-inflation extension, `reviews/r053-covering-core-independent-review.md`, SHA256 `f466f85312c98c116f8276a10b14fbde3f8d9870be1ac5f970c4e196d76ece66`. These are inherited mathematical results, read and used but not rerun in R61.

## 6. Exact implementation and actual scope

The saved inputs use positions in (1/4)Z² and masses in (1/16)N. For D0=2^bits, the radius floor for a bbox of integer mass q is computed by

    k=floor(D0 sqrt(q/16)/2)=isqrt(q (D0/8)^2).

The implementation checks the exact squared inequalities, uses (k-1)/D0 and (k+1)/D0 as the inner/outer radii, and discards a true empty bbox only after the exact test span_quarters²>=q. All sweep areas, subset sums, ratios, flows, and post-run Gamma values use integer/Fraction arithmetic. Finite endpoint arrangements ignore only zero-area boundaries.

The 12 completed batches used Python 3.13.18, each one process. They evaluated 768 cases (725 distinct inputs up to source reordering and 43 repeated controls), with 8–12 sources. The exact ledger has 1,217,024 cuts: 397,434 compact and 819,590 noncompact. The most restrictive computed cut is proper in 430 cases and noncompact in 165. All saved input-specific bounds are below B. R53 did not numerically lower Z(S) on this particular dataset: the recorded number of improved cuts is zero. Its role here is the proved general envelope and compact-cut pruning, not an invented measured improvement.

The highest bound in this limited test family is 611497158639875/525566558076928 (batch 7, case 46). Its critical mask has 10 of 11 sources and span²/mass=256/239>1. This is an upper certificate for that prescribed finite input at threshold 1, not an extremal search result and not a bound on C2. Execution receipts are separate artifacts. Independent audit remains pending; no completed independent audit result is claimed.

The maximal R54/R59 boundary coefficients are family-specific necessary conditions and have not been used as arbitrary source-cut bounds. The R56 L1 baseline can be chosen freely by future allocations; this atomic implementation does not impose a maximal identity baseline or claim to discretize that saturated L1 example without error.

Open: prove a useful uniform compensation estimate for the noncompact cuts in (9), or exploit additional spatial structure to certify them without exponential enumeration. The planned 16–20-source structural frontier was not implemented or run after the user's new team limit. No new general C2 bound, commit, or push is claimed by this delegated branch.
