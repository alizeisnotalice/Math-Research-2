# R82: a uniform true near-halo payment for all masses up to one half

2026-11-16. Centered axis-parallel squares with one common radius, strict
threshold 1. This is an exact two-atom lemma used by the separate shared-side
theorem. It does not assume a maximum block or replace a residual halo.

## 1. Statement and exact rectangles

Put m=1/4, d=71/200, e=1/1000. For 0<=a<=1/2 let

    alpha=m delta_(0,0), beta_a=a delta_(-d,0),
    D_a=H(alpha+beta_a)\H(beta_a),
    C_e=[-1/2,-e]².

Then the actual all-radius deletion satisfies the uniform strict bound

    |D_a intersect C_e|<699/10000=.0699<.07.          (1)

No claim about all eligible source blocks follows from (1) alone. The
separate core proof owns its extension to variable old width, a position
rectangle and arbitrary eight outside masses under shared global caps.
This appendix proves only the base near-cost lemma, without a mass grid or
parameter B&B. R79--R81 frozen files remain unchanged.

Write positive receiver coordinates X=-x,Y=-y and use the crop
C=[e,1/2]². Let

    v=sqrt(m)/2=1/4, z=sqrt(a)/2, u=sqrt(m+a)/2,
    s=u+z, so u²-z²=v²=1/16.

For a>0 the complete full halo consists of the two singleton receivers and
the mixed receiver. An open receiver label is equivalent to the existence
of one common radius strictly between its farthest source distance and
sqrt(label mass)/2. The exact captured label of any strict witness proves
the reverse inclusion. Hence after subtracting the complete singleton
residual, the relevant rectangles in the positive orthant are

    A0=(0,v)²,
    A1=(d-u,u)x(0,u),
    R =(d-z,d+z)x(0,z),
    D_a intersect C=[(A0 union A1)\R] intersect C.      (2)

All radii, including every residual radius, are represented by these
receivers. Their strict boundaries can be closed for area calculations
only because the boundaries are finitely many line segments. A receiver
of zero width is empty, including at a source-span equality. Extra captured
mass helps a label witness; it need not be the exact captured label.

At a=0 the residual is empty and the answer is directly
(v-e)²=.062001, satisfying (1). A formal zero-mass near source is not needed.

For all remaining masses, u>=v>d/2, u<.434<1/2, and z<.354<d-e. Thus the
residual left endpoint d-z is above e. Also d+z>u because u-z=v²/s<=v<d.
Its intersection with A0 union A1 therefore has exactly the area

    (s-d)_+(z-e)_+.                                   (3)

Indeed u>=z places its left endpoint to the right of A1's left endpoint,
and its top below A1's top. If z<e the residual misses the crop altogether.
In that case u<=v+z and hence s<=v+2z<v+2e<d, so the first positive part
vanishes as well. These endpoint cases cannot be inferred by extending the
R81 calculation that assumed z>e.

## 2. The two full-union formulas and their matching boundary

If u<d-e, the near mixed rectangle has left endpoint d-u>e. It covers
exactly the portion of A0 to the right of that endpoint; the remaining old
strip has area (d-u-e)(v-e). Thus, setting

    C0=m/2+dv-ve+e²=.213501, L=d+v+e=.606,

the old-plus-near full area is

    (2u-d)(u-e)+(d-u-e)(v-e)=C0+2z²-Lu.              (4)

If u>=d-e, its left endpoint is at most e and it contains all the cropped
old square. The full-union area is instead (u-e)². Accordingly the exact
deletion cost G(a)=|D_a intersect C_e| is

    C0+2z²-Lu-(s-d)_+(z-e)_+,   if u<d-e;
    (u-e)²-(s-d)_+(z-e)_+,     if u>=d-e.             (5)

At u=d-e the full-area formulas agree. At s=d the residual contact has zero
width, so its contribution is zero. At z=e its cropped height is zero.
These conventions account for closed crop boundaries, open receivers and
all endpoint masses without adding a line of strict halo points.

## 3. Branch I: the residual has not reached the near receiver

Suppose s<=d. Since u(s)=(s+v²/s)/2 is increasing on s>=v, this is
strictly before the coverage change u=d-e: at s=d one has

    u_d=(d+v²/d)/2, z_d=(d-v²/d)/2,
    u_d<d-e.

The full cost is (4). As a function of z with u=sqrt(v²+z²), its derivative
is

    z(4-L/u)>=z(4-L/v)>=0.

Here 4-L/v=197/125>0. Thus the maximum in this whole branch, including
a=0 and z<e, occurs at the rational contact s=d. Direct substitution gives

    G_contact=691638777/10082000000<699/10000.         (6)

## 4. Branch II: contact, before the old strip disappears

Now suppose s>d and u<d-e. Contact implies z>e by the bound following (3).
Use the rational change of variables

    u=(s+v²/s)/2, z=(s-v²/s)/2.

Substituting in (5) gives exactly

    G2(s)=22737/125000-249s/2000-961/(32000s)
             +1/(512s²).                             (7)

Let P(s)=(7/100-G2(s))s². It is the polynomial

    P(s)=249s³/2000-13987s²/125000+961s/32000-1/512.

At d=71/200 its first Taylor data are

    P0=14101223/80000000000,
    P1=-937829/400000000,
    P2=41393/1000000>0.

For x=s-d>=0 the exact identity is

    P(s)=P0+P1 x+(P2/2)x²+(249/2000)x³.

Dropping the nonnegative cubic and completing the square proves

    P(s)>=P0-P1²/(2P2)
      =291048892263/2649152000000000>1/10000.         (8)

This is a polynomial positivity certificate on the whole half-line s>=d,
not a finite collection of mass slices. Throughout the actual mass domain,
u<=sqrt(3)/4<.433013 and z<=sqrt(2)/4<.353554, so s<1. Combining (7)-(8),

    7/100-G2(s)=P(s)/s²>1/10000,

which proves G2(s)<699/10000 in this branch.

## 5. Branch III: the near receiver contains the old cropped square

If u>=d-e=177/500, then s>3/5. Indeed u(s) is increasing on s>=v, and

    u(3/5)=169/480<177/500.

In particular s>d and z>e. The second formula in (5) simplifies to

    G3(s)=-s²/4+71s/400+31073/500000
                -357/(32000s)+1/(1024s²).           (9)

For s>=3/5 its derivative obeys

    G3'(s)=-s/2+d/2+(d+2e)v²/(2s²)-v⁴/(2s³)
       <=(d-3/5)/2+(d+2e)v²/[2(3/5)²]
       =-1757/19200<0.

Consequently the cost is below the formal value at s=3/5 (evaluation of the
rational function there is valid even though that point precedes this
geometric branch):

    G3(s)<=G3(3/5)=9038149/144000000<699/10000.        (10)

Equations (6), (8) and (10), together with the explicit a=0 calculation,
prove (1) for the entire closed mass interval. There is no uncontrolled
transition at residual contact or at complete coverage of the old square.

## 6. Verification and use in the shared-side theorem

`computation/round082-dual/check_near_uniform_constants.py` reads
`computation/inputs/r082-dual/near-uniform-constants.json`. It uses only
integer/Fraction arithmetic to derive the Laurent-polynomial identities
(4), (7), (9), their matching boundary, the contact value, the Taylor data,
positive lower bound (8), derivative bound and all rational comparisons.
It retains a=0 and z<e explicitly. The result and executed-command,
input/script/software/PID/timing receipt are saved in
`computation/results/round082-dual/checks.json`.

There is no mass sampling, root solving, interval subdivision, halo union
integration, or parameter B&B. The analytic arguments above cover each
entire branch. The earlier source-specific R81 cost formula is extended
here only after verifying the previously excluded zero and high-mass cases.

The separate core R82 proof owns the Lipschitz/crop transfer and shared-side
budget for variable width, near/far distances and arbitrary eight masses.
Those conclusions require their source geometry and common global cap
constraints. Neither this local two-atom lemma nor that structured class
settles arbitrary internal mass profiles or the general maximum-block
selection problem. No MATLAB, new agent, commit or push is used.
