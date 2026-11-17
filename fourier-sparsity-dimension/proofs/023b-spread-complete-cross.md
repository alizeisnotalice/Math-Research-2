# Complete subdouble partitions with unbounded cross parameter

For k>=2,q=2^k take a k-spread of X=F2^(2k): q+1 linear k-flats L_j meeting pairwise only at0 and partitioning X minus0. Desarguesian spreads exist by identifying X with F_q^2 and taking the q+1 one-dimensional F_q subspaces. Their orthogonals also form a spread for any fixed nondegenerate F2 pairing, since L_i+L_j=X.

Choose eps_j in {+1,-1}, q/2 negative and q/2+1 positive, so sum eps=1. Define F(0)=1 and F(x)=eps_j for nonzero x in L_j. For nonzero frequency xi, exactly one L_j has xi in L_j^perp. Thus

 WF(xi)=1-sum_j eps_j+q eps_j=q eps_j,
 WF(0)=1+(q−1)sum_j eps_j=q.

Hence F=(-1)^f is bent. Pair each negative L with a different positive L, leaving one positive L_*. Each A_i=L_i△L_i' has weight2q−2 and parameter mu=k; the final block is L_* of weightq. All q/2+1 blocks are disjoint and cover X. Flipping any paired block exchanges its two eps signs and preserves their sum0. Flipping L_* changes both F(0) and the remaining eps from+1 to−1, leaving F(0)=sum eps. The same Walsh calculation gives bentness for every independent block flip.

Let W span all block indicators and V=W+RM(1,2k). All f+V are bent. W has dimension m=q/2+1, and its intersection with affine functions consists exactly of constants. Indeed after subtracting a constant, any linear ell in W is constant on each punctured L_j. For k>=2 choose nonzero u,v,u+v in a given L_j: linearity forces that constant to be0, hence ell vanishes on every L_j and on X. Consequently

 dim V=m+(2k+1)−1=q/2+2k+1.

This constructs genuine complete partitions for arbitrarily large mu=k. It is not a high-rank A28 counterexample; dimension is O(q). The supports of all block partial spectra are the analogous paired cross blocks of the orthogonal spread. Block flips only change their signs, while adding a linear phase translates all spectra by the same character. Thus selected-patch projected normal rank is exactly2k, and the remaining parameter kernel is the m-dimensional constant-per-block space.

Independent construction and finite validation passed. The main agent reran the archived final script: k3 all2048 phases; k4 all512 block-code phases plus5120 affine-basis translates (5632 checks including duplicates); k5 4097 deterministic samples, seed20261117. Generator ranks11,17,27 agree with the formula. These finite checks do not exhaust the k4/k5 full codes; the sign-sum identity proves the general claim.

The two full punctured spectral arms are not affine flats. Filling their apices would make them intersect. R23 instead selects one affine hyperplane avoiding the apex within each arm; those actual selected patches are disjoint and their direction intersection is zero. No affine completion is used in the rank argument.
