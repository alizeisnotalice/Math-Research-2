# Recovered local n=24 terminal batch

This directory records a completed local `n=24,P=131072,R=64,L=16,strat`
batch.  Its GPU computation completed normally, but the original driver
postcheck rejected fifteen printed values of `F-2G-4D` that were nonzero by
only approximately 1e-216 from final floating arithmetic.  The per-group
integer identity is exactly zero in every group.  The frozen aggregator uses
a documented scale-aware check for that printed field and independently
checks the integer identities, metadata, seed range, terminal invariant, and
distance histograms; it accepted all 64 records.

The finite complete-law diagnostic has lift ESS about 1.45 of 64 and maximum
replicate share about 0.808.  It is kept to expose this concentration, not
as a precise probability estimate or an asymptotic result.

The separately archived normalizer--distance consistency audit also shows
that this batch's plug-in \(A_2/B\) and complete-law distance histogram
violate an exact first-merge necessary inequality.  Thus even before assigning
a mechanism to the discrepancy, its reported Q and tail values cannot be
used as exact original-law values.
