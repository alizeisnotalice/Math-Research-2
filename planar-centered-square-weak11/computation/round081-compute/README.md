# R81 shared-side alternating family: cap interval and two actual controls

This component evaluates one fixed twelve-position family at its two
approved feasible endpoints. It does not sample the parameter interval,
run uniform parameter B&B, or claim a general four-corner theorem.
The already frozen R79 and R80 components remain unchanged.

## Why the shared side is the state

R63 retains the actual full-minus-residual deletion and its exact deficit.
R68/69 already rule out reachable-state and universal net-mixed shortcuts;
R70 covers balanced corner supports only in a smaller width interval.
R72/73 exclude transferring diffuse deletion through ordinary atomization
and extracting a fixed area payment solely from bounded unit-cap slack.
R74 proves corner localization but refutes independent pasting and pointwise
D4 preservation. R80 further refutes individual corner half-mass payment.
These inherited results were selectively read, not rerun; source paths,
hashes and reading scopes are in `source-reading-ledger.json`.

For a normalized general fat bbox [0,w]x[0,h], the four corner squares have
total area 4(w-1/2)(h-1/2). The complementary fifth region has area
3(w+h-wh). Paying that fifth region by its full area leaves aggregate quota
K=1/2+3(1-w)(1-h). The sufficient target is sum of actual corner deletion
<=K, not a separately imposed quota at every corner.

Our control has w=h=999/1000 and four old corner atoms of mass 1/4.
Start with near=(-71/200,0) and far=(0,-171/200), then rotate both through
the four corners by R(x,y)=(w-y,x). There are eight outside positions.
Index a corner by i=0,1,2,3. Its near and far masses are

```
i even: a_i=1/4,       b_i=t
i odd:  a_i=1/2-t,     b_i=1/4.
```

The side unit window captures its two old endpoints (mass 1/2), near_i
and far_(i-1). Thus a_i+b_(i-1)<=1/2 is necessary, and is saturated by this
family for every t. The source windows of different corners are disjoint,
but these physical square caps couple their masses. Every endpoint output
records the actual four closed unit windows and all captured source indices.
The total outside mass is 2 and total mass is 3, independently of t.

## Exact cap feasibility and actual qualification

The fixed positions have six coordinate values in each axis. Enumerating
all 21x21 closed axis-interval products gives 441 candidate rectangles,
which reduce to 157 distinct nonempty full bbox captures. The independently
rebuilt residual has 104 distinct nonempty captures from its own 441
candidates. No 2^12 source-subset enumeration is used or needed.

Each captured mass is affine in t. Intersecting the exact inequalities
mass(t)<=max(1,bbox-span^2) gives the raw cap interval [0,1/2]. Strictly
positive outside masses require 0<t<1/2; both requested endpoints
t=1/4 and t=499/1000 are therefore legal, and the requested interval is
unchanged by feasibility clipping. No replacement or extra parameter point
was needed. The cap interval covers the continuum by linear inequalities;
the area calculation covers only its two requested endpoints.

At both endpoints there are 25 eligible full bbox captures. Maximum eligible
mass is 1, with all five maximum masks
15,99,390,1548,2073. Mask 15 is the old block; the others are the four side
blocks. These are all maximum eligible sets because a maximum set must be
bbox-closed when every mass is positive. Eligible hyperedges make the input
connected, and its nonterminal gap is exactly 45725491/11000000 >0.
Ties are recorded rather than discarded. The R72 common-contraction theorem
explains why a strict universal-budget counterexample could be made unique;
no such extra contracted input was computed here.

## Complete actual areas and five-region credit

Every full and residual bbox receiver uses its actual captured mass and
actual bbox. Exact integer-isqrt root brackets have width at most
1/(1000*2^72). Full labels containing no old source cancel against the
identical complete residual halo before bounding actual deletion. All
old-involving outer receivers are contained in C=[-1/2,w+1/2]^2.

The five regions are SW,SE,NE,NW corner squares of side 499/1000 and
N=C minus those four squares. Strict open receiver boundaries have zero
area. For every region, the code saves

```
I = area of the internal halo,
T = area of internal halo intersect complete residual,
G = area of full halo outside (internal halo union complete residual),
D = I-T+G.
```

Here the internal halo has exact rational endpoints. Its retained part
plus G equals D at each directed integer endpoint, not merely after decimal
rounding. A separate postprocess reclips the two saved geometries and also
records full halo, complete residual, uncovered-full and true deficit in
all five regions. It introduces no input and does not recompute caps.
Both passes use the frozen R75 vertical-event/interval-union algorithm and
the distinct compressed Cartesian midpoint-coverage algorithm; all endpoint
integer areas agree.

| t | Strong/first corner D | Other corner D | Four-corner sum | Whole-plane D |
|---|---:|---:|---:|---:|
| 1/4 | .067120965521119 | same at all four corners | .268483862084476 | 1.024812854025325 |
| 499/1000 | SW,NE: .129916563043447 | SE,NW: .062198302395511 | .384229730877916 | 1.183831162339845 |

These decimal displays summarize rational enclosures, which establish both
aggregate sums strictly below K=500003/1000000. At the upper endpoint the
two strong corners reproduce the true R80 local half-mass failure, while
the two weak corners leave enough reserve. Using the equal geometric quota
K/4 only as bookkeeping, total corner excess is approximately .0098316261,
total reserve .1256048952, and reserve minus excess .1157732691. At t=1/4
the excess is zero and reserve is approximately .2315191379. The diagnostics
store the exact directed intervals; these equalities alone are not a new
payment theorem.

The upper-endpoint per-corner credit illustrates the retained mechanism:
each strong corner has I=.062001, T=.036105, G about .1040205630; each weak
corner has I=.062001, T=0, G about .0001973024. Thus the two strong local
gains are compatible with a much smaller net gain in the other corners.
All values refer to one and the same outside measure at that endpoint.
No independently optimized corner measures are pasted together.

## Runtime, reproduction and scope

The one main batch used Python 3.13.18, PID 26771, exit 0, helper wall .138s
and script .049s. The saved-geometry regional postprocess used PID 27824,
exit 0, helper wall .089s. Both are observed absent in the component manifest.
There were no pool workers, BLAS calculations, random seed, numerical solver
tolerance, new agent, MATLAB, commit or push.

From the topic directory, using fresh output/input destinations:

```
/Users/zzh/.venvs/main/bin/python3 computation/round081-compute/shared_side_family.py --family computation/inputs/r081-compute/family.json --output-dir NEW-RESULT-DIR --endpoint-input-dir NEW-ENDPOINT-INPUT-DIR
```

The postprocess has the frozen production path as its explicit input and
refuses to overwrite its diagnostic output; reproduce it in a copied topic
directory or choose a new output path in a copied script. It is not an
independent audit replay. Source, family, generated exact endpoint inputs,
all results and both actual run receipts are included in the component.

Frozen computational dependency: R75
`computation/round075/certify_three_atom_collapse.py`, SHA256
`a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`.
The separate theory agent has proposed a conditional uniform payment for
this fixed family; that proof and its review are separate artifacts.
These two endpoint computations neither prove its continuum conclusion
nor settle the general four-corner aggregate budget or C2.
