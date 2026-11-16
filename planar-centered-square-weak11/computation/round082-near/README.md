# R82 whole-segment geometry and robust constants only

This component is the independent small exact geometry/constants check for
the separately owned R82 shared-side Borel-measure theorem. It does not
compute a halo area, sample outside masses, repeat the base near-only
three-branch proof, or run parameter B&B.

The actual theorem has old atoms of mass 1/4 at the four corners of [0,w]^2,
3/4<w<1. In each rotated corner there is an arbitrary nonnegative finite
measure supported on the near axis segment with distance [.352,.358] and
another supported on the far segment with distance [.85,.86]. Additional
nonnegative rho is carried outside the four open effective source windows.
All original global square caps are assumed. This kernel checks named
support geometry and constants; the Borel area conclusion is analytic.

## One common window for each entire support segment

The code represents w and eight source-distance variables by exact affine
forms over their full closed interval box. The closed verification range
w in [3/4,1] safely includes the actual open width range; the proof still
needs w<1 for strict eligibility and separation.

A fixed side window in the reference frame is
[-.86,.14] x [0,1]. It contains both old side endpoints and the entire near
bundle at one end and far bundle at the other. A fixed corner window is
[-.358,.642] x [-.86,.14]. It contains the old corner and both entire local
bundles. Four rotations give eight unit capacity windows. The four rotated
open P=(-1,w-1/2)^2 windows supply the exact effective source localization.

Every capacity window depends only on w and fixed support endpoints, never
on a selected distance variable. This is checked in code. Thus one window
contains all points of the relevant segments simultaneously; the quantifier
is not a moving window chosen separately for each representative point.

For each of 12 windows, the code checks all four old-point and eight segment
position families, giving 144 incidences and 576 signed affine facet bounds.
Every affine minimum and maximum is obtained exactly from coefficient signs
and interval endpoints. Membership of the expected named families and a
strict excluding facet for every other named family are stored in full.
This proves continuous-range named-support membership, not a finite point
sample. Endpoints in closed unit windows are retained; expected membership
in each open P is certified by strictly positive margins.

Extra rho is not among these twelve named support families. It may enter a
capacity window; its nonnegative mass only strengthens the necessary
inequalities. Its absence from corner witnesses follows from the separate
analytic containment of every r<1/2 corner square in P, together with rho's
support assumption. The kernel does not enumerate rho or prove a separate
Borel area-transfer statement.

The shared windows give necessary constraints
a_i+b_(i-1)<=1/2, B_i=a_i+b_i<=3/4 and sum_i B_i<=2. Each of the eight
bundle masses occurs exactly once in the summed side constraints. These
constraints are not asserted sufficient for the global square caps.

## Independent rational constant checks

The base analytic estimate G<699/10000 is imported from the separate
near-only proof and is not recomputed here. The kernel verifies:

- u<7/16, z<3/8, hence u+2z<19/16<6/5;
- crop enlargement from e=.001 to zero adds 999/1000000;
- distance displacement is at most 3/1000, yielding near comparison cost
  <74499/1000000<.0745<.075;
- sqrt(3)>433/250 gives far full-rectangle upper bound 139/20000<.007;
- the joint hinge threshold is 189/400 and slope upper bound 23/80;
- active-count bounds for k=0,...,4 are
  0,111/400,111/200,233/400,11/100, with maximum 233/400;
- the joint sum bound is 5359/32000<.168;
- combining the coarse .075 near and .007 far bounds with this joint sum
  gives 3171/6400<62/125=.496<1/2.

For a bundle, u+2z is the factor for comparing a full outer rectangle at the
closest support endpoint and an actual residual inner rectangle containing
the entire bundle. It is not a claim that an arbitrary bundle is replaced
by one atom while preserving global caps. The external proof owns the
rectangle inclusions and area comparison, including zero mass and an empty
residual inner rectangle; this kernel checks their geometric ranges and
rational constants only.

## Actual execution and binding

One actual run used Python 3.13.18 standard-library Fraction arithmetic,
PID 38291, exit 0, helper wall .082 seconds and script .0029 seconds. The
manifest records observed process absence. There is no random seed or
numerical tolerance. No second input or run, parameter area evaluation,
agent, MATLAB, commit or push was used.

The previously prepared source hash d77f409b... was never mathematically
executed. Before the actual run, only scope metadata was added to make the
named-family/rho distinction explicit. The actual executed source is
`fa6f1dc5ca9da4bef0a0508532639957c30271b8ff0e4325f6c5325132f83bb3`;
the input remains
`d7d7dea287c0074c782fe9edb740816ec2ce3767fb5a9540597eefdca4977f01`.
The receipt and output bind these exact executed bytes.

The run also binds the reviewed formal proof
`proofs/r082-robust-shared-side-arbitrary-mass-payment.md`, SHA256
`6817649e38760727974aec8fb97d6022f0a3889924186c0d310406df36ad5fbd`.
That proof, rather than this finite kernel, supplies the Borel theorem.

From the integrated topic directory, using a fresh output path:

```
/Users/zzh/.venvs/main/bin/python3 computation/round082-near/check_robust_geometry.py --input computation/inputs/r082-near/geometry-ranges.json --formal-proof proofs/r082-robust-shared-side-arbitrary-mass-payment.md --out NEW-OUTPUT.json
```

The kernel refuses to replace an existing output. The base near-only scalar
author check is a different component and was not duplicated here.
