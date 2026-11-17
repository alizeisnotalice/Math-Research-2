# R22: local rigidity of all legal blocks below twice the minimum weight

Scope: X=F2^(2k),q=2^k,k>=2; f,f+1_A bent, q<=|A|<2q. This is a local necessary classification, not an A28 solution. KTA forms are from SST Theorem1.8 and Lemma1.10 (PDF pages2–3), read 2026-11-17.

## Excluding type I mu>=3

The type-I support is B times F2^r in an affine container of dimension k+2, r=k+2−2mu>=0. B={Q=1}, Q the hyperbolic quadratic on 2mu coordinates, |B|=2^(2mu−1)−2^(mu−1). Fourier-normalize each free-frequency slice by 2^r; coefficients a are real, supported B, |a|<=1, and Wa is 0 or +/-T, T=q/2^r=2^(2mu−2).

If m is the number of nonzero Walsh coefficients, Parseval gives
m <= 2^(2mu)|B|/T^2 =8−2^(3−mu)<8.
Thus m<=7. But the characters indexed by any m<=7 frequencies are linearly independent on B^c: their Gram diagonal is D=2^(2mu−1)+2^(mu−1), and all off-diagonal entries have absolute value 2^(mu−1), by the explicit quadratic Gauss sum. The smallest eigenvalue is at least 2^(mu−1)(2^mu+2−m)>0 for mu>=3. Hence no nonzero combination of these characters can vanish on B^c. Every free slice must be zero, contradicting nonzero signs on A. Type I mu>=3 is impossible. Type I mu1 is minimal, mu2 was R21.

## Type II, all mu>=3

After translation, the core consists of (u,0),u!=0 and (0,v),v!=0 in F2^mu x F2^mu. Put N=2^mu. Free dimension r=k−mu. A normalized real slice has arm coefficients a(u),b(v), each in [-1,1], and zero at 0. Its Walsh is Ahat(xi)+Bhat(eta), which belongs to {0,+N,-N} for every independent xi,eta.

If either arm transform is constant, its mean is its value at physical0=0, so it is zero. The other arm cannot have a nonzero coefficient of magnitude N, because Parseval gives total squared Fourier mass <=N(N−1)<N^2. Thus the slice is zero.

Otherwise both arm transforms are nonconstant. Differences within either value set are multiples of N. A range of 2N would force the other value set to be a singleton, impossible. Thus both have exactly two values at distance N. Write Ahat=N1_E−r, with r=|E| between1 andN−1; this follows from its zero mean. Parseval gives sum_u a(u)^2=r(N−r)<=N−1, forcing r=1 orN−1. Thus a(u)=+chi_alpha(u) or -chi_alpha(u) for every u!=0. For Bhat the analogous multiplicity s satisfies r+s=N, since its smallest sum must be -N. Consequently b(v) has the opposite sign and is another character. Each nonzero slice is fully saturated on both arms.

Parseval along free coordinates forces exactly one common nonzero free-frequency slice. Hence every admissible phase on A is a fixed opposite-arm reference plus an affine function of its (k+mu)-dimensional container. This extension is unique, since the punctured arms span the container affinely (mu>=2), and restrictions of all affine functions have dimension k+mu+1. Relative directions of any common bent family are affine on A.

The intrinsic core spectrum is itself the cross {(alpha,eta):eta!=beta} union {(xi,beta):xi!=alpha}, up to translations. Ambient partial spectrum is the symmetric difference of two affine k-flats with intersection dimension k−mu. Splitting each punctured arm into mu affine pieces gives 2mu pieces, and local kernel dimension1. A complete independent block partition therefore has a crude contribution 2mu+1 per type-II block, combined with R2. This alone has an unbounded mu factor and does NOT yield O(q).

## Computation

Full mu3 core sign scan in type2_core.py: 16384 patterns; exactly128 admitted, exactly one affine restriction orbit. Finite check only. Two independent audits agree with the type-I exclusion; the full type-II argument has an independent audit. The dyadic complete-partition bound below has a separate independent audit. A common-rank argument avoiding the 2mu cost remains unresolved. Both local and MacBook CPU/Metal finite checks passed; their full legal-mask sets match Python.

## Explicit complete-partition bound (limited scope)

Let m1 minimum blocks, m3 three-half blocks, and t_mu second-type blocks (mu>=3) form a complete independent partition. The dyadic partition of F2^mu minus a point has one affine piece of each dimension 0,...,mu−1. In particular each arm includes a point piece. After lifting with annihilator and parameter variables these are affine cylinders. A direction in the annihilator of projected normals must leave every such piece invariant; the point piece forces the full intrinsic translation to vanish. Thus only one constant per block survives, with no hidden mu-dimensional local kernel. Applying R2 to the total number of cylinders gives

    dim V <= 2 m1 + 4 m3 + sum_mu (2mu+1)t_mu − 1.

The mass identity is m1+(3/2)m3+sum_mu (2−2^(1−mu))t_mu=q. For only mu3 in addition to minimum/three-half blocks, the maximal cost per unit mass is 7/(7/4)=4, so dim V<=4q−1. For unbounded mu this estimate has a logarithmic loss and is not the target O(q). Completeness and independent indicators in V are essential hypotheses.

## Source scope and finite checks

The classification input is Schillewaert–Storme–Thas, *Minimal codewords in Reed–Muller codes*, author manuscript, Theorem 1.8 and Lemma 1.10, PDF pages2–3: <https://cage.ugent.be/geometry/Files/342/MinimalRevision.pdf>. The theorem and weight formulas were directly reread this round; the original KTA proof was not independently reconstructed. The downloaded-source hash is recorded in references/r021-second-weight-source.md.

The exact local audit gauss_and_cross_audit.py verifies the quadratic complement transforms and Gram margins for mu3..6, as well as opposite-arm reference spectra. cross_partition_audit.py checks all intrinsic translations for mu3..6 (64,256,1024,4096), reconstructs the dyadic pieces and computes projected-normal ranks by F2 elimination. These are finite checks, separate from the general argument.

## Every second-type parameter has actual bent-pair witnesses

For mu>=3 and r>=0 let k=mu+r, x,y in F2^mu and z,w in F2^r. Set

    f=x·y+z·w+1_{x=0,z=0},
    g=x·y+z·w+1_{y=0,z=0}.

The first is MM under (x,z)|(y,w), the second under (y,z)|(x,w), with identity permutation. Thus both are bent. Their difference support has z=0 and exactly one of x,y zero; w is free. Its size is 2(2^mu−1)2^r=(2−2^(1−mu))q. Hence large-mu second-type legal blocks cannot simply be ruled out. This is not a high-dimensional common-code construction or an A28 counterexample.
