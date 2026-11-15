# R68: a residual-halo-free central rectangle

2026-11-15. Finite positive atomic inputs, coincident atoms merged, threshold1 and strict centered axis-parallel common-radius square witnesses. Let L>0 be the maximum eligible mass in A, and let S be a block of massL. Write H=H(A), R=H(A\S), D=H\R. After translation, bbox(S)=[0,w]×[0,h], with w,h<sqrt(L). The theorem below is a general geometric restriction. It does not establish a3.5 deletion or erosion-union selection bound.

## 1. Maximum-block exclusion applies to every residual witness

For every nonempty source set B outside S,

    d(S∪B)²≥L+m(B).                                    (1)

Otherwise S∪B would be eligible with mass greater than the defining maximumL. This is the full-set version of the previously established maximal-block mass cap, not a new independent principle.

Now take any actual strict witness for the residual input at center x, with half-side r and captured residual set T of mass b. Its strict inequality gives

    0<4r²<b≤L.                                        (2)

The last bound holds because T itself is an eligible source subset of A: d(T)≤2r and b>4r². It does not require T to be bbox-closed in the full input or to avoid S geometrically; only its source indices are residual indices.

Put

    f_L(b)=sqrt(L+b)−sqrt(b)/2.

Then x must lie outside the closed rectangle

    K_b=[w−f_L(b), f_L(b)]×[h−f_L(b), f_L(b)].          (3)

To prove this, suppose x lies in K_b and put rho=sqrt(L+b). In the horizontal coordinate, let u and v be the least and greatest T coordinates. The horizontal span of S∪T is exactly

    max{w, v−u, v, w−u}.

Each term is strictly below rho. Indeed w<sqrt(L)<rho; v−u≤2r<sqrt(b)<rho; and using x_1∈[w−f_L(b),f_L(b)] and r<sqrt(b)/2 gives

    v≤x_1+r<rho,
    u≥x_1−r>w−rho.

The vertical argument is identical. Hence d(S∪T)<rho, contradicting(1). The strict witness slack handles the boundary of K_b as well; this is an actual set exclusion, not merely an almost-everywhere claim.

## 2. A universal central rectangle and exact local states

For every b≥0,

    f_L(b)≥sqrt(3L)/2,                                 (4)

with equality exactly at b=L/3. For a direct sign proof, both sides of

    2sqrt(L+b)≥sqrt(b)+sqrt(3L)

are nonnegative. Squaring reduces it to

    L+3b≥2sqrt(3Lb),

which is (sqrt(L)−sqrt(3b))²≥0. Thus no inadmissible squaring step or numerical root sign is involved.

Since f_L(b)≥sqrt(3L)/2, the rectangle

    K=[w−sqrt(3L)/2, sqrt(3L)/2]
         ×[h−sqrt(3L)/2, sqrt(3L)/2]                    (5)

is contained in every K_b. Every point of R has a residual strict witness, so(3) proves the general theorem

    R∩K=empty.                                        (6)

In particular, for the full erosion unions W and W′ defined in R67,

    W′∩K=W∩K=H∩K=D∩K.                                (7)

This follows from D⊂W′⊂W⊂H and(6). The central region has no old-halo surplus in either construction. All overpayment in W or W′ is confined to the complement of K. This statement uses the actual remaining halo, not a proxy subtractor.

The coefficient sqrt3/2 in this universal rectangle cannot be increased for arbitrary maximum blocks. Take S to be a singleton at0 of mass1, and a residual singleton at(2/sqrt3,0) of mass1/3. Their union has mass4/3 and squared span4/3, so it is not strictly eligible and L=1. The left edge of the open residual halo is exactly sqrt3/2. Any larger coefficient makes the corresponding central rectangle intersect that halo. This is sharpness of the exclusion constant in the unrestricted class, not sharpness of a deletion bound or a statement in the fat/rank-four branch.

## 3. Exact peripheral area bound in the fat branch

Let s=sqrt(L) and let C_1=(-s/2,w+s/2)×(-s/2,h+s/2). The R67 support containment places W and W′ inside C_1. K and C_1 have the same center. Their side lengths are respectively

    (sqrt3 s−w, sqrt3 s−h),
    (s+w, s+h).

If w,h>3s/4, then K is contained in C_1; in fact the weaker coordinate conditions w,h>(sqrt3−1)s/2 suffice. In the fat branch, (6) gives

    |R∩C_1|≤|C_1\K|
      =(s+w)(s+h)−(sqrt3 s−w)(sqrt3 s−h)
      =(1+sqrt3)s(w+h)−2L
      <2sqrt3 L.                                      (8)

The same upper bound applies to the two overpayment regions W∩R and W′∩R. It is not a3.5 bound for W or Delta: neither the true occupied area in K nor the peripheral deletion is paid by (8). Adding independent worst-case estimates for them would discard the geometry again.

For arbitrary smaller w or h, K need not be contained in C_1. The correct general peripheral expression is instead

    (s+w)(s+h)
      −min(s+w,sqrt3 s−w) min(s+h,sqrt3 s−h),             (9)

obtained by intersecting the two same-center rectangles. Using the simplified fat expression outside its domain could even give a negative purported bound and is invalid.

The useful positive interface is (7) together with the explicit four peripheral portions C_1\K. A subsequent payment argument must couple the captured S mass and the residual mass in these portions, or use a geometric support reduction. The six-source attachment in this round shows why simply assigning all peripheral surplus to H(S), even after restricting to reachable complete N-states, is false.

## 4. Scalar verification scope

Root's actual Wolfram run independently verified nine scalar statements: the minimum and equality case of f, scaling, the area expansion, fat-domain containment and positivity, and the strict2sqrt3L consequence. Mathematica15.0.0 ARM, kernel0.162705seconds, wall1.238075seconds, exit0, owned PID12617 checked absent. The three-file symbolic supplement is bound by `root-symbolic-manifest.json`, SHA `8573f08ad72255f68fd4d8f70801c7b1ad0ed77e56b2f474f07f9bf39fc15542`. These are supplementary scalar checks. The source-to-receiver argument in section1 and exact local state equality(7) are proved above and are not supplied by the computer algebra output.
