# Exact n=7 central-hole CRT certificate

This addendum completes the four modular residues for the single central
selected hole `v=24` of the 7-by-7 grid.  The local `p0` checkpoint is paired
with three finished remote checkpoints (`p1`--`p3`); the latter's exit and
resource receipts are retained.  The reusable modular source is the sibling
`../punctured_layer_mod.cpp` (SHA-256
`5f93fa39695e0f4badabcbabaad2befab920fa1590ea2337fe47a09201587db9`).

`verify_n7_crt.py` independently rebuilds the 49-node Pratt certificate, an
integer independent-set layer histogram after deleting the central vertex and
its four neighbours, both CRT formulae, and the strict uniqueness comparison
with `48!`.  Its output certifies only the two selected-hole integers
`A2v` and `Bv`; it is neither a full `A2(7)` computation nor an asymptotic
result.  In particular, `5 Bv - A2v > 0` is a fact about this one hole.

The separate `../n7-partial/` record remains intentionally present: its
local second-prime attempt had an unknown nonzero termination.  It did not
produce a checkpoint and is not used as a cause diagnosis; the completed
remote `p1` residue is the one used here.
