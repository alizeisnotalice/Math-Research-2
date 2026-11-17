# Corrected modular short-tail implementation

Only the corrected source, its final aggregate finite cross-checks, boundary
test, and n=7 opt-in guard are retained.  The earlier `candidate-small/`
outputs used a directed-edge union bug and are deliberately excluded.  No
binary is archived.  The formula itself is justified independently in the
sibling `r029-short-tail-ie` certificate; this implementation alone is not a
proof of it or of an all-n tail bound.
