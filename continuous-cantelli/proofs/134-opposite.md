# The opposite wing in the saturated tilt budget

Let f be a bounded nonnegative continuous exact graph. Suppose on the negative side, as x decreases to zero through positive values,

    f(−x)=x²/[2(ℓ+b log ℓ+d+o(1))],   ℓ=log(1/x).

No rate or smoothness is required of the o(1). Write L=log z and g₀=1/√(2π). As z tends to positive infinity the entire opposite wing of the exact tilt has mass

    ν_z((−∞,0))=(g₀/z)[1−1/(2L)+(b log L+d+γ_E)/(2L²)+o(L⁻²)],

where γ_E is Euler's constant. Thus for any fixed r>0 the expression on the right is an unavoidable contribution to the deficit 1−ν_z((0,r)). More precisely the deficit equals this opposite-wing mass plus ν_z([r,∞)). If f>0 on [r,∞) and liminf at +∞ exceeds one, the second mass is exponentially small in z². In this case the deficit on the favored positive side determines the first two contact parameters on the opposite side, even if no contact expansion of the favored side has been assumed. The reflected statement applies independently.

**Proof.** Set u=z|x|. Relative to g₀/z, the integrand on u>0 is bounded by exp(−u). On L⁻³≤u≤6 log L the denominator D=L−log u+b log(L−log u)+d+o(1) equals L−log u+b log L+d+o(1) uniformly. There g(u/z)/g₀ differs from one by O((log L/z)²). The two omitted tails have total integral o(L⁻²); the same holds for the integrated expansion terms by integrability near zero of u²|log u| and exponential decay at infinity. Expanding exp(−u²/(4D)) through second order gives the integrated constant, first and second coefficients from

    1 − u²/(4L) + [u²(−log u+b log L+d)/4+u⁴/32]/L².

The uniform denominator error o(1) here contributes o(L⁻²) after integration against u² exp(−u). All other integrated Taylor remainders are O((1+log L)²/L³), using the finite exponential moments of polynomials in u and |log u| with the appropriate u² weight. Since ∫u²exp(−u)du=2, ∫u⁴exp(−u)du=24 and ∫u²log(u)exp(−u)du=3−2γ_E, the constants reduce to the formula above. For the positive far wing under the additional gap hypothesis choose M>r with f≥1+ε beyond M. Gaussian completion gives an upper bound exp(−εz²/2) there, and the positive minimum on [r,M] gives exp(Mz−cz²) there. Exactness provides the claimed equality for the deficit.
