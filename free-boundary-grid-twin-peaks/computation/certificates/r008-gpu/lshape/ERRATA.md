# Erratum: boundary multiplier

The original `smc.jsonl` recorded the correct raw `p` and `logZ`, but its derived `scaled` field used a hard-coded multiplier 6.  The specified L shape has initial boundary 7.  Therefore every original L `scaled` entry must be multiplied by `7/6`; the raw run is retained unchanged.  The corrected source computes the boundary directly from the start mask and asserts `B==7` for this run.
