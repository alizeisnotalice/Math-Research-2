# Bounded continuation after the beam-depth correction

The old corrected 16-seed receipt (`corrected-results.json`) was produced by
source SHA-256 `9e26ee2b...`; that exact source is retained as
`mm_affine_flat_gpu_search_prebeam_9e26ee2b.py` (file SHA-256
`9e26ee2b...`, byte-for-byte equal to the source hash embedded in that
receipt).  It is not rewritten by this continuation.

The continuation source SHA-256 is `598cc6a...`.  Unlike the earlier beam, it
does not start from the empty code for only `9-r` layers after random greedy.
For each of four continuation runs it starts from the first six independent
directions of a verified random-greedy basis, expands with beam width 32 toward
quotient dimension 9, and records the actual depth.

Parameters: four runs, 96 random orders per run, beam width 32, prefix
dimension 6, target dimension 9.  All four beams reached dimension 8 after
two added layers and had no legal ninth extension in the bounded candidate
search.  The respective arbitrary-XOR exact-oracle cache sizes were 60,846,
63,589, 61,501, and 47,689.  This is a bounded continuation from selected
prefixes, not an exhaustive exclusion of ninth directions or a general upper
bound.  The command retained the default RNG root `20261117`, so these are not
four independently newly sampled RNG roots.  Run 0 has the same MM seed as
run 0 in the old 16-run receipt; after that, the corrected beam consumes a
different random stream and runs 1--3 use different recorded MM seeds.  They
must be described as four continuation repeats, not four independent seed
families.

`beam-continuation-4seed-verification.json` independently rechecks the four
stored bases: independence modulo `C` and all 1,024 W-slices with Python
integer FWHT and zero tolerance.
