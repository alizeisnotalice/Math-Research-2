# R78 actual prefix computation, one fixed input

The original R77 corrected S20 input is used without moving any atom:
400 old points `(21i/400,21j/400)`, `0<=i,j<20`, each of mass 1/400;
q1=(-71/200,0) of mass 1/4 and q2=(0,-171/200) of mass 1/2.
The minimal input is a byte-identical copy of the frozen R77 input, SHA256
`d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a`.
Its inherited `scope` string describes the R77 qualification run; this new
component additionally computes the actual halo difference on Gamma only.

The region is the closed square Gamma=[-1/2,-1/400]^2. Halo witnesses use
strict open squares and mass > 4r^2. Gamma boundary inclusion has no effect
on area. The frozen R77 all-square caps exclude all witnesses with r>=1/2.
For x in Gamma, an old atom in an open witness must have each coordinate
<199/400, so only i,j<=9 can occur. The nonempty old capture is a Cartesian
origin-anchored prefix with 1<=k,l<=10. The menu therefore consists of
100 old prefixes times all four outside subsets, plus the three pure-outside
labels. A label may also capture additional atoms: it remains a sufficient
witness, and all actual capture labels occur in the menu. The separate
theory component supplies the full completeness proof.

Every label retains its actual bbox. With mass M and bbox
[xmin,xmax] x [ymin,ymax], its open center rectangle is
(xmax-sqrt(M)/2,xmin+sqrt(M)/2) x
(ymax-sqrt(M)/2,ymin+sqrt(M)/2), empty if M<=span^2.
The full menu includes pure old and pure outside labels. The residual is
the same complete set of three pure-outside labels. No internal mass is
collapsed and no mixed-only surrogate replaces the full halo.

All roots are enclosed by exact rational intervals of width at most
1/(400*2^72). Integer area lower bounds use full inner rectangles minus
residual outer rectangles; upper bounds reverse both directions. Each
endpoint is evaluated twice: by vertical x-event/direct y-interval merging
and by compressed Cartesian midpoint-cell coverage. The exact integer
answers agree. Both union routines are imported read-only from the frozen
R75 core, whose hash is checked at import.

The full exact fractions are in `actual-prefix-area.json`. Decimal displays:

| Quantity | Area |
|---|---:|
| Actual complete Gamma deletion | .007652366201647474 |
| Pure-old labels after complete residual | .00050625 (exact 81/160000) |
| Old plus exactly one outside, after complete residual | .004241193795207637 |
| Old plus both outside, after complete residual | .0029049224064398374 |

The certified full upper endpoint is below 1/125, hence strictly below both
half the old P mass, 1/8, and the main corner threshold 80003/640000.
The group areas are individual unions, not additive credits asserted by
the computation.

Of the 100 labels for each outside mask 0,1,2,3, respectively 1,22,0,9
have nonempty open centers in Gamma. The nine joint-active labels are
(k,l)=(1,1),...,(9,1). **Each of the nine has strictly positive area outside
the complete residual halo.** Their individual areas are recorded with
directed endpoints and are not summed because they overlap. Their union
after residual subtraction is the last row of the table. None is classified
as wholly residual-covered or unresolved.

This is one exact actual-input regional computation. It explains the sharp
loss between the frozen R77 collapsed proxy (.1291124425...) and the actual
input (.0076523662...). It does not improve a general C2 bound or establish
a theorem for other inputs, and no parameter search or full-plane 400-source
halo was run.

## Reproduction

From the topic directory, using the existing Python 3.13.18:

```
/Users/zzh/.venvs/main/bin/python3 computation/round078/certify_actual_prefix.py --input computation/inputs/r078/n20-corrected-actual-input.json --out NEW-OUTPUT.json --bits 72
```

Use a new output path to preserve the frozen production result. The actual
run used the `python-math` run helper; `process-receipt.json` preserves its
exact command, Python version, source hash, PID 94876, exit 0 and timing.
The manifest records observed process absence. Random seed and numerical
tolerance are inapplicable: this is exact deterministic integer/Fraction
arithmetic. No MATLAB, new agent, commit or push was used.

Frozen dependencies:

- R75 exact union core `computation/round075/certify_three_atom_collapse.py`,
  SHA256 `a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`.
- R77 qualification `computation/results/round077-mixed-proxy/n20-qualification.json`,
  SHA256 `c6e1671967a52587f3e91da303420bda1654bc82d07b64fad080c62c637e2180`.
  This supplies the independently audited complete caps, unique maximum,
  connectivity and nonterminal status; those global checks were not repeated.
