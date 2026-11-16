# R74 independent review: actual corner payment and the limits of localization

2026-11-16 (Asia/Shanghai). Sol independent review, strict threshold 1,
centered axis-parallel squares, common radius. New files only; no commit or
push. The unrestricted finite-atomic corner conjecture and exact C2 remain
open. This review distinguishes analytic proof reading, new exact executions,
and producer executions that were not repeated.

## 1. A positive structural theorem

I read the entire proof `proofs/r074-single-effective-outside-corner-payment.md`
in the geometry component, including its all-radius and residual arguments.
Its mathematical Sections 1–5 pass independent review. The final proof SHA256
is `d5e1dd01b5a390ed7cd832909f244de0ddf9175ecb0c96669e15dd365cb3fa39`;
its component manifest is
`computation/results/round074-two-atom/component-manifest.json`, SHA256
`4e5490e456c71c7a7713528a2cb25c2337953e46adfdf5d89c90f5501f2309f1`.

The decisive overlap lemma has no capacity hypothesis. For an interval of
length u, shrinking a centered interval's radius from u to v in [0,u] loses
at most u-v of overlap. If its outer and inner overlap lengths in two axes
are A,B and a,b, set d=u-v. Then a>=(A-d)_+, b>=(B-d)_+. If A<=d or B<=d,
the annulus area is at most du<=u²-v². Otherwise its excess below u²-v²
is d(2u-A-B)>=0. This proves the square-intersection annulus bound for every
outside position, including boundary contacts and empty intersections.

For m delta_0 plus M delta_q, the all-radius mixed receiver is the
intersection of the two open squares of radius u=sqrt(m+M)/2. In the
southwest orthant it lies in (-u,0)². Subtracting the actual M-singleton halo
uses v=sqrt(M)/2 and costs at most u²-v²=m/4 by the preceding lemma. The
old singleton contributes at most another m/4. The same M occurs in both
radii; this is a true residual payment, not an independently selected
inner/outer radius bound. Coincident positions and M=0 cause no gap.

For arbitrary finite internal mass northeast of the receiver orthant's
vertex, moving all its mass to that vertex dominates capture on every
square centered in that orthant. It enlarges full halo while retaining the
same single residual atom. This collapse need not preserve caps: the target
two-atom theorem needs none. Thus the orthant deletion is at most m/2.

Apply this only after the actual global caps establish strict witness
r<1/2 and localization to P_SW=(-1,1/2-e)², e=1/n². Its old mass is
m=floor(n/2)²/n². If each of the four reflected windows has at most one
positive outside location, with arbitrarily many additional sources outside
those windows, all four local payments apply to the same global input.
Consequently, for every n>=2 in this structural class,

    |D| <= 3-3/n^4+2 floor(n/2)^2/n^2 < 7/2.

No uniqueness, connectedness or nonterminal assumption is needed for this
upper bound. The original global square caps are still required for
localization. Multiple effective outside positions in one window remain
unhandled; applying the mixed-label payment separately would charge the
same old mass repeatedly.

The producer's 255 one-dimensional and 4,335 two-dimensional regressions,
and its 10-atom closure fixture, were inspected as source/results, not rerun.
The new independent SymPy execution in `payment-identities.json` checks six
zero identities for the annulus branches, same-M cancellation, container
area and grid parity. These identities supplement, rather than replace,
the continuum argument above. Python 3.13.18, SymPy 1.14.0; PID 93929 exited
0 and was confirmed absent.

## 2. Localization preserves halo membership, not independent feasibility

I read the entire core proof
`proofs/r074-corner-localization-and-borel-half-mass-obstruction.md`.
Its localization, compatibility and Borel arguments pass independent review.
The final proof SHA256 is
`406a08d287ff3478ac9466dbe7323283abffc77409a2130ffc0e19a5eeeced49`.
For a receiver in the closed corner [-1/2,-e]², every r<1/2 square lies
strictly in P_SW. Sources on the boundary of P never participate. Restricted
full and residual measures inherit the global caps by domination, so neither
can create an omitted witness at r>=1/2. Hence both halo memberships, and
their difference, are exactly preserved on the corner at all radii.

The feasible restricted input retains the full S_n. The smaller sigma|P is
only a local captured measure and need not be a local maximum block. The
finite fixture in Section 4 below already demonstrates this: for n=3,
sigma(P)=1/9, whereas the origin plus q form an eligible label of mass
77/450>1/9. The full input still has unique maximum S_3 of mass 1.

The core's distinct n=5 pasting counterexample is correct. Each outside
measure (8/25)delta_{(-12/25,6/25)} and its reflected partner at
(-12/25,18/25), separately with full S_5, keeps S_5 uniquely maximal,
connected and nonterminal. The mixed-subset horizontal-span argument covers
every old subset. Their union puts mass 31/25 in one closed unit square,
and its captured label has bbox span 24/25. Thus disjoint effective source
windows still share square-cap constraints.

D4 mass averaging is a different operation. Because sigma_n is invariant,
linear square caps survive averaging admissible inputs. It is the deletion
objective's pointwise preservation or an unjustified Jensen step that fails;
the finite example below establishes precisely that limited claim.

## 3. Independent exact replay of the Borel obstruction

For sigma supported in a box of widths at most 1 and outside density between
0 and 1, all unit-square caps imply all-square caps. Inside a square of side
R>=1, one can choose a contained unit square including every captured old
source. There is no old mass in the difference, whose outside mass is at
most R²-1. R<=1 follows by unit containment. The coordinate interval used
to choose this unit square is nonempty even at equality of widths. This
argument does not apply to an atomic outside measure.

For n=5 and nu=(100/111)1_(P_SW\J) dx, I independently reconstructed the
entire exact certificate from raw old coordinates and P/J intersections,
without importing producer code. `review_r074_borel_corner.py` checks:

- 144 closed unit-origin vertices, all 529 point/open-interval strata and
  their 1,156 endpoint limits. On open strata the old count is fixed and
  the continuous outside area is bilinear; closed vertex counts dominate
  neighboring atomic limits. Outside the checked origin domain an axis
  captures no old source and density<=1 suffices.
- All 1,655 lower cells using literal maximum distances of the 25 old atoms
  to each whole receiver cell and exact minimum P\J overlap. The smallest
  strict margin is 2209/284160000 and the largest saved radius is 131/320.
  Positive margin and strict separation from 1/2 and P's boundary permit
  a common slight radius increase to open witnesses.
- All 2,161 exclusion cells, enumerating all feasible old prefixes from
  k,l<=n rather than trusting saved prefix rows. This executes 3,283 upper
  tests. Every radius and receiver coordinate is lowered in the safe
  direction in the increasing cost. The other 280 cells stay undecided.

The exact independent enclosure is

    175099/2048000 <= |D intersect Gamma_SW| <= 204723/2048000.

The lower end exceeds sigma(P)/2=2/25 by 11259/2048000. The upper end is
less than the main proposed single-corner threshold 631/5000. Density<1
makes the actual residual halo empty, so this is a Borel deletion enclosure,
not a full-halo surrogate. The binding unit origin (-13/25,-13/25) has mass
9/25+(100/111)(444/625)=1 and proves sharpness of this fixed density family.

The producer result pin is
`ac49ff47032f44997cd220ed597e066d4c2ce92ab9478bd11a2f0f9bf962e912`.
The new replay used Python 3.15.0 Fraction, PID 85837, exit 0 and confirmed
absent after wait, about 0.394 seconds wall time. Full result and process
receipt are `single-corner-density-replay.json` and its `-receipt.json`.
This disproves only the Borel half-local-mass extension. It gives no finite
atomic counterexample: a finite approximation's residual halo cannot be
discarded or transferred by weak convergence.

## 4. A finite strict receiver-cell obstruction to pointwise symmetrization

Let n=3, q=(-1/3,-1/4), and nu=(3/50)delta_q. Its D4 average has eight
distinct outside positions, each of mass 3/400. The new exact fixture
`review_r074_symmetry_fixture.py` enumerates literal bbox closures for both
inputs: 46 for the original 10 atoms and 247 for the averaged 17 atoms.
All caps pass; their eligible counts are 28 and 50. Both have unique
maximum S_3 (mask 511), one eligible component, and positive nonterminal
gaps 3293/8910 and 11543/8910.

On the closed cell K=[-201/1000,-199/1000]², area 1/250000, the original
origin-plus-q label has a common open witness radius 101/500 and full excess
4441/562500. Its residual atom misses the cell's halo with a positive
minimum-distance margin 22609/2250000. Thus K is contained in original D.

For the average, the only effective local sources are the origin and q and
its swapped point; all other sources are outside strict P_SW. The script
tests all seven nonempty subsets using the exact minimum radius over K,

    delta_K(U)=max(span_x(U)/2,span_y(U)/2,
                   x_lo-min_x(U),max_x(U)-x_hi,
                   y_lo-min_y(U),max_y(U)-y_hi,0).

Every mass is strictly less than 4delta_K(U)². Caps already exclude
r>=1/2, so no full witness was omitted. The entire K is outside the
averaged full halo, hence outside its deletion. This is a finite admissible,
unique-maximum, connected, nonterminal failure of receiver-cell preservation
under D4 averaging. It makes no total-area comparison, no assertion of
nonconcavity of the total corner objective, and no violation of either
corner bound. Python 3.15.0 Fraction, PID 87451, exit 0, confirmed absent;
about 0.036 seconds wall time. Result and process receipt are saved.

## 5. Scope and reproduction

Only the three small executions described above are new independent runs.
No old R72/R73 complete halo computation, adaptive Borel classification,
producer regression, MIP or MATLAB batch was repeated. The negative finite
searches and localization controls remain producer evidence with their own
recorded scope; they do not prove the unrestricted bound.

From the topic directory, with the producer certificate available:

```sh
python3 -B computation/review_r074_borel_corner.py \
  /path/to/single-corner-density-certificate-result.json \
  /Volumes/data/CodexTemp/T2/r074-borel-replay.json
python3 -B computation/review_r074_symmetry_fixture.py \
  /Volumes/data/CodexTemp/T2/r074-symmetry-fixture.json
/Users/zzh/.venvs/main/bin/python3 -B computation/review_r074_payment_identities.py \
  /Volumes/data/CodexTemp/T2/r074-payment-identities.json
```

The final report `rounds/074.md`, SHA256
`4a40457aa87bcfec7764ed1b20689fd212e4f94088f92dd9dd95d0c51762d6d0`,
was read in full and passes scope and mathematical review. The core manifest
is SHA256 `d7d446418a96dcf9f35346456c5a7c81dbf21cbc262daf20d2d0ff6f88f1fe32`.
Final source pins are also recorded in this component's manifest. Hash
binding identifies versions; the mathematical review and actual execution
scopes are the separate statements above. Root's unified integrator owns
full payload/dependency hashes, staged bytes and link checking; no duplicate
generic source-binding program was introduced.
