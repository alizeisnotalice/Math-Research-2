# R77 mixed-only collapse diagnostics

This component certifies only three specified three-source proxies. All use
strict open square receivers at threshold 1. Full labels are exactly
`{old,q1}`, `{old,q2}`, `{old,q1,q2}` (masks 3,5,7); the old singleton is
excluded. The complete residual halo uses masks 2,4,6. Empty receivers are
discarded by the exact test mass > squared bbox span before clipping.

| Specified input | Region | Certified area, decimal display | Decision |
|---|---|---:|---|
| old mass 1; q1=(-28/25,0), mass 4; q2=(0,-161/50), mass 8 | entire SW orthant | .522033990724985 | >1/2 |
| collapsed old mass 1/4; q1=(-31/100,0), mass 1/8; q2=(0,-83/100), mass 1/2 | Gamma=[-1/2,-1/400]^2 | .11072244672167358 | <1/8 |
| collapsed old mass 1/4; q1=(-71/200,0), mass 1/4; q2=(0,-171/200), mass 1/2 | same Gamma | .129112442502507 | >1/8 |

The full rational endpoint fractions in `directed-areas.json`, rather than
these rounded decimals, establish the comparisons. Radius intervals have
width at most 2^-72 divided by the rational coordinate scale. For each endpoint,
a vertical x-event/y-interval union difference and a Cartesian compressed
midpoint-cell algorithm give identical exact integer areas. Both functions
are imported read-only from the already frozen R75 file, pinned below.

The first case gives an uncapped mixed-only counterexample. Its exact formula
is 78/25 - 53sqrt(5)/50 - 161sqrt(13)/100 + 161sqrt(3)/50. The short root
bounds in `fixed-formula-checks.json` give the lower bound 652519/1250000 >
261/500 > 1/2. The original old block is not maximum: the whole three-source
bbox has mass 13 and squared span 25921/2500 < 13.

For the original n20 proposal the residual pair is empty: 5/8 < (83/100)^2.
Subtracting a negative rectangle area caused the earlier hand estimate to
look too large. The failed proposal is retained separately; no connector or
additional parameter search was run for it.

For the corrected proposal the exact proxy formula is
43sqrt(3)/200 - 101sqrt(2)/1600 - 24641/160000. The elementary bounds
sqrt(3)>433/250 and sqrt(2)<283/200 give 206483/1600000 > 1/8.
The proxy total mass is 1, so every strict witness has radius <1/2.

The original, uncollapsed S20 input is separately supplied with all 402
positions/masses. Here h=21/400, every old atom has mass 1/400. The 440
internal-prefix cap constraints and three pure-outside constraints hold for
both n20 proposals. For the corrected actual input a second computation
enumerates all 45,302 distinct nonempty coordinate-rectangle closures,
recovers each tight bbox, and verifies every cap mass <= max(1,span^2).
There are 7,386 eligible closures; S20 is the unique eligible block of maximum
mass 1. The eligible S20 block connects all old sources, the pair origin-q1
connects q1, and the eligible q1-q2 pair connects q2. The nonterminal gap
span(all)^2 - (35/33) mass(all) is 8319673/5280000 >0.

The 100 actual old-prefix+both-outside labels within the SW source window
have 9 strict Gamma witnesses and 91 null witnesses. These are full-label
activity tests, with no residual subtraction; they do not give actual
S20 corner deletion area. No 400-source halo was computed.

The collapse also destroys a unit cap. The closed unit square
[-71/200,129/200] x [-171/200,29/200] captures original old mass 39/400,
of which 30/400 lies in P. The original full mass is 339/400; after all
P old mass 100/400 is moved to 0, the full mass becomes 409/400.
Thus the corrected example disproves the mixed collapse proxy even when
its pre-collapse input is legal, unique-maximum, connected and nonterminal.
It does not assert an actual S20 corner-area violation, a 3.5 violation,
or a new bound on C2.

## Reproduction and dependencies

Run from the topic directory with the existing Python 3.13.18 interpreter:

```
/Users/zzh/.venvs/main/bin/python3 computation/round077/certify_mixed_proxy.py --input computation/inputs/r077/mixed-proxy-cases.json --out computation/results/round077-mixed-proxy/directed-areas.json
/Users/zzh/.venvs/main/bin/python3 computation/round077/check_n20_qualification.py --input computation/inputs/r077/n20-corrected-actual-input.json --out computation/results/round077-mixed-proxy/n20-qualification.json
/Users/zzh/.venvs/main/bin/python3 computation/round077/check_fixed_formulas.py
```

For archival replays, select new output names or a copied topic directory;
the frozen production outputs should remain unchanged. The actual three
production runs used `python-math/scripts/run_math.py`; their exact arguments,
source hashes, interpreter versions, PIDs and successful exit codes are in
the three process receipts. The first two processes were independently
observed absent by the final deterministic check. The component manifest
records final absence of all three. There was no random batch, new agent,
MATLAB invocation, commit or push.

Frozen dependency: `computation/round075/certify_three_atom_collapse.py`, SHA256
`a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`.
The analytic derivations are in the separately owned R77 theory proof;
formula interval overlap is a consistency check, while the two exact union
algorithms directly establish the geometric area inequalities.
