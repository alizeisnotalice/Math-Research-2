# R69: pure mixed area alone can exceed one fifth of the deleted mass

2026-11-15. Finite positive atomic inputs, strict threshold1, centered axis-parallel squares with one common radius in both coordinates. Let S be a maximum eligible source block of mass L. Write I=H(S), R=H(A\S), and P=H(A)\(I union R). The example below proves |P|>L/5 in the remaining fat/rank-four/support-tree-failure domain. It does not refute a bound for |P|-|I intersect R|, a bound for true deletion, or general3.5 selection.

## Exact positive input and maximum block

Use S with the following four sources:

|Index|Position|Mass|
|---|---|---|
|0|(0,0)|2499/5000|
|1|(99/100,99/100)|1/2|
|2|(0,99/100)|1/10000|
|3|(99/100,0)|1/10000|

Add q=(-3/10,0) of mass1/2. Then m(S)=1, m(A)=3/2, and S has both widths99/100. Thus S is strictly eligible and fat. A source set containing q and any right-column source has horizontal span129/100, whose square exceeds the entire mass3/2. Such a set is ineligible. A set containing q and only left-column sources has mass at most9999/10000<1. All proper subsets of S have mass below1. Therefore S is the unique maximum eligible block, L=1.

S is an eligible edge containing its four sources, and the pair{0,q} is eligible, so the whole input is eligible-connected. It is nonterminal because

    (129/100)^2 > (35/33)*(3/2).

Every support-tree side edge has cost99/100 and every diagonal cost9999/10000 at L=1. Hence the minimum tree cost is297/100>5/2. The four-corner support union has area(199/100)^2>7/2. This example is not eliminated by the preceding tree or filled-support sufficient conditions; neither of these quantities is a claimed lower bound for the actual halo.

## A whole positive-area patch in P

Let r=4999/10000 and t=177/500. Define the two open rectangles

    U_plus = (-r,r-3/10) x (t,r),
    U_minus = (-r,r-3/10) x (-r,-t),
    U=U_plus union U_minus.

At every x in U, the single common half-side r captures source0 and q. Their total mass is4999/5000, while

    4r^2=24990001/25000000 < 4999/5000.

Thus U is contained in H(A). The same radius works throughout both open rectangles, so this is not an independent-axis or center-dependent-radius calculation.

Every x in U has vertical distance greater than t from q and source0. Since

    4t^2=62658/125000 > 1/2 > 2499/5000,

neither singleton can supply a strict witness there. The residual input consists only of q, hence U avoids R.

To exclude every internal S witness, first consider either right-column source. Its horizontal distance from U exceeds7901/10000, so the area of any square capturing it exceeds

    4(7901/10000)^2 > 1=m(S).

No internal witness can capture a right-column source. After these are excluded, only the left-column mass2499/5000+1/10000=4999/10000 is available. Capturing the top-left source from either rectangle requires radius greater than4901/10000, and

    4(4901/10000)^2 > 4999/10000.

That is impossible for a strict S witness. The only remaining possible internal nonempty source set is the origin singleton, already excluded by the vertical distance t. Consequently U avoids I as well, and U is a subset of P.

Its exact area is

    |U|=2*(2r-3/10)*(r-t)
        =5105041/25000000
        =1/5+105041/25000000 >1/5.

This is a strict rational positive-area certificate for |P|>L/5. It only disproves paying pure mixed area by a universal L/5 allowance before accounting for internal/residual overlap or the actual internal-halo slack.

## Actual arithmetic and portable reproduction

Root ran `check_pure_mixed_obstruction.py` with Python3.15.0 and Fraction arithmetic. It checked all31 nonempty subsets, the unique maximum, connectivity witnesses, nonterminal and tree conditions, every entire-patch inequality above, and the exact patch area. The process exited0, wall0.059690seconds, and owned PID23690 was checked absent. The source SHA is `4a88fa61e7d9880010aca9c0714e8901f2448d4bc094522a4a287e210c72360f`. The script and generated exact input, output and process receipt are retained together. No complete halo area or optimization was run for this construction.

For a fresh reproduction, copy the standalone script to a new directory below `/Volumes/data/CodexTemp/T2`, then run `python3 -B` on that copy. It writes the two JSON outputs beside itself, requires mounted data, and rejects disabled assertions. No historical workspace dependency is imported. The analytic coverage and exclusion argument above supplies the proof; the exact calculation checks its rational signs and source-domain conditions.
