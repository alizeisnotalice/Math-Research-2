# Rank-biased fixed-prefix move

The modified move selects ranks 0 and 1 more often while preserving the
uniform law on any fixed valid ordered prefix space.  Exact Fraction
enumerations verify detailed balance for all P2 prefix ranks and P3 ranks
through 5, both for this choice and for generic fixed rank weights.  A small
GPU trace at n=3 is replayed independently at the proposal level; its
verification file records accepted and rejected root-moving categories.

For a legal root--nonroot active swap, the new proposal probability is the
old one times `r/8+1/2` (strictly larger only when `r>4`); for a root--outside
replacement it is multiplied by `r/4+1/2`.  This is a local transition comparison only.  It proves neither
faster mixing nor lower SMC variance.  This output uses a new implementation
identifier and cannot be accepted by the frozen r023 aggregator.
