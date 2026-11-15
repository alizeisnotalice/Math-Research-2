# R69: exact heavy-core certificates on three existing fixtures

2026-11-15. Apply the frozen general theorem
`r069-heavy-core-mass-sensitive-payment.md` to three already fixed inputs.
There is no new source-family search or full-halo/erosion-area computation.
All have L=s=1, a unique maximum S, strict positive masses, one connected
nonterminal component, fat bbox widths and failed support-tree criterion.
The script checks these assertions by all63,511,127 nonempty subsets.

The R69 minimum-input file copies the mathematical positions/masses from
its individually hashed source inputs, normalizing their different JSON
schemas. Its new hash does not claim byte identity with those source files.
In particular, the seven-source input is the core's frozen net-mixed
mechanism, not another independently generated input.

| Existing input | Core H | M | E | r | beta | Gamma | Gamma+9(E+beta) |
|---|---|---:|---:|---:|---:|---:|---:|
| R68 final six-source surplus | its two main sources |499/500|1/500|43/50|0|8377/2500|4211/1250|
| R68 eight-point heavy-corner tree obstruction | heavy corner |393/400|7/400|17/20|0|4309/1600|4561/1600|
| R69 seven-source net-mixed mechanism | two main sources |499/500|1/500|43/50|0|7927/2500|1993/625|

Every final number in the last column is strictly below7/2. These are
proved upper bounds for W and Wprime through the new class, not the actual
areas of either state. The support-tree minimum costs are respectively
27/10,21/8,27/10, all strictly above5/2. Thus the mass-sensitive class pays
cases not paid by the R68 tree criterion.

## Exact clipping and near-source checks

For the six-source input, C=[-1/5,17/10]^2 and K=[17/50,29/25]^2.
Their intersection area is1681/2500 and the core-support intersection
area is1033/2500. Hence Gamma=8377/2500. Both outside sources have squared
distance to V greater than E, so beta=0.

For the heavy-corner input, C=[-1/2,11/8]^2 while
K=[1/40,69/40]^2 extends outside C. The actual intersection has area729/400,
while Z_H intersection K has area1. Therefore

    Gamma=225/64-729/400+1=4309/1600.

The nearest outside-to-V distance is31/60, whose square exceeds E=7/400.
This case also satisfies the theorem's uniform arbitrary-rank corner
condition E+beta<=1/32, independently giving constant2767/800.

For the seven-source input, C=[-1/2,7/5]^2 and
K=[-9/25,43/50]x[-43/50,43/50]. K again extends outside C. Its clipped
area is1037/625, whereas Z_H intersection K has area61/50. Thus

    Gamma=361/100-1037/625+61/50=7927/2500.

The two outside squared distances to V are81/100 and4/25, both aboveE.
This explicitly tests why subtracting |K| instead of |C intersection K|
would be invalid. The script independently integrates the full predicate
C AND(not K OR core-support) over rational coordinate cells and obtains
the same Gamma as the clipped formula.

## All-label branch tests without a W-area run

The exact all-subset eligible labels meeting S split into core-touching,
low-mass missing-core and high-mass missing-core counts:

    six-source:12,2,0;
    eight-point:160,7,0;
    seven-source:17,3,1.

These are all eligible subsets, not a claim about the number of bbox-closed
production labels. Every low label is included in V union B. For every
high label the script checks M+m>L, the H/T span inequalities, union
ineligibility and exact emptiness of E_T(m) intersection CLOSED K.
The latter uses rational comparisons sqrt(m)>q via q<=0 or m>q^2,
so no floating square-root boundary is involved. The seven-source high
label makes this test nonvacuous.

The low rank refinement from the general theorem is also recorded. It
can reduce the certificate using min{9(E+beta),rank(V)E}; it does not change
any of the original three sources or introduce another geometry run.

Actual final execution used Python3.13.18 and SymPy1.14.0,0.017405seconds
internally,exit0. Launched PID30634 was verified absent afterward. The
three symbolic identities concern only the rational corner constants and
3-times-area factor; the previously established scalar minimum was not
rerun in CAS. All parameters, exact areas, source pins, branch labels and
near profiles are saved in the result/receipt. No MATLAB, LP, full halo,
W/Wprime area, general3.5 selection claim, commit or push is included.
