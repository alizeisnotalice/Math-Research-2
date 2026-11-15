# Dedekind cycle-type source audit — 2026-11-15

## Verified author-hosted source

Keith Conrad, University of Connecticut, [*Galois Theory at Work*](https://kconrad.math.uconn.edu/blurbs/galoistheory/galoisapp.pdf), PDF p. 16, Theorem 4.29.
The retrieved author-hosted PDF has 35 pages and PDF creation/modification
metadata `2013-01-16`; its SHA-256 in this temporary audit directory is
`d394dc4230ebb8af5b113e519b41b74db066528e880e305eca92a740823c1fa7`.

The theorem assumes that `f in Z[X]` is monic and irreducible over `Q`, of
degree `n`, and that the prime `p` does not divide `disc(f)`.  If the distinct
monic irreducible factors of `f mod p` have degrees `d_1,...,d_k`, it concludes
that the Galois group, acting on the roots, contains a permutation of cycle
type `(d_1,...,d_k)`.  Page 18, Corollary 4.34 explicitly identifies the type
`(1,...,1,2)` as a transposition.

## Proof dependency checked

Conrad's one-page author note [*Galois Groups over Q and Factoring Mod p*](https://kconrad.math.uconn.edu/blurbs/gradnumthy/galois-Q-factor-mod-p.pdf), Theorem 1 on p. 1, gives the Frobenius proof: take a prime ideal above `p` in the splitting field; the discriminant condition gives unramifiedness and distinct reduced roots; Frobenius acts by the `p`th-power map; roots of an irreducible finite-field factor of degree `d` form a Frobenius orbit of length `d`.

## Scope for the current use

This cited formulation requires `Q`-irreducibility.  It should therefore be
applied only after an independently verified irreducible reduction modulo a
prime and Gauss's lemma have shown that the concrete `C` is irreducible over
`Q`.  For monic integral `C`, a squarefree reduction modulo `p` is equivalent
to `p` not dividing its discriminant, so it supplies the theorem's good-prime
condition.  This note does not use the source to assert the more general
reducible-squarefree version.
