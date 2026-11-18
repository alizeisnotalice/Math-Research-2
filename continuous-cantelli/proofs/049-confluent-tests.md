# Finite confluent tests and uniform transform approximation with a fixed density gap

Fix b>0. For any positive integer m, finite set of real nodes and prescribed finite derivative orders at each, there are bounded even smooth variance graphs v_R from round 44 with the following properties as R tends to zero. The output Y_R has the same first 2m moments as N(0,1+b); its real MGF and the prescribed finite derivatives match at all the nodes. Every positive-order output density derivative matches its Gaussian value at zero. Yet the density value error at zero tends to p_1(0), while the full output law tends to N(0,1+b) in total variation. All real derivatives of each fixed order of the normalized Cantelli transform converge uniformly on the entire real axis to the corresponding derivatives of the constant function 1. The constants may depend on the prescribed finite tests and b.

## Replacing finite value tests by derivative tests

Use the exterior inverse construction of round 44. By symmetry group ±t together and include all orders up to each requested order. At zero keep orders 2,4,...,2m, increasing m if necessary. Order zero and the odd derivatives there are automatic. Along an even bump ψ the derivative of the output MGF with respect to the bump amplitude at the constant variance b is the entire function
    (t²/2)e^(bt²/2) ∫ψ(x)g(x)cosh(tx)dx.
At each t>0 differentiation up to its assigned order gives an invertible triangular transformation of the rows with integrands x^j cosh(tx) (even j) and x^j sinh(tx) (odd j). At t=0, the rows of orders 2,...,2m in even steps are invertible triangular combinations of 1,x²,...,x^(2m−2). All these functions are independent on (1,2): if a linear combination were zero, analytically continue in x, divide by exp(t_max x) and let x→+∞. The polynomial multiplying the largest positive exponential must vanish identically, removing all its coefficients; repeat and finally compare the polynomial at frequency zero. Consequently one can select fixed independent bump locations and supports as before. The finite derivative residuals on the core are O(R) uniformly, by compactness in x and bounded variance. Thus the small exterior adjustments still satisfy ||v_R−b||∞=O(R) on their fixed support. The unchanged fixed-point construction supplies all the jet assertions. One can also impose the conditions on Y_R's normalized transform because multiplication of its MGF by a fixed nonvanishing entire factor is triangular on each finite jet.

## Actual discrepancy of the laws

On the core |x|≤R the total variation cost under mixing is at most its source probability, O(R), with TV normalized as a supremum over events. The transition adds O(R⁵). On the exterior controls use ∫|∂_s p_s| ≤ C/b for s near b and the O(R) amplitude estimate. All other sources have identical conditional laws in the constant model. This proves TV=O(R).

The density at zero needs its smooth representative. The core contributes g(0)+∫_(|x|≤R)g(x)p_(x⁴w(x))(−x)dx, whose off-diagonal part is exponentially small (w lies in a fixed positive compact interval). In the baseline the core contributes O(R). In the transition the estimate sup_s p_s(x)≤C/|x| bounds its actual contribution by O(R⁴); the removed baseline there costs O(R⁵). The fixed outer bumps affect the density at zero by O(R). Hence
    D_R(0)−p_(1+b)(0)=g(0)+O(R).
At every fixed nonzero y the core contribution converges to zero by separation in x; exterior controls and the shrinking transition also contribute vanishing discrepancies. The discrepancy at the center persists despite full-axis transform convergence.

## Uniformity on the full real parameter axis

Normalize F_R(t)=e^(−(1+b)t²/2) E e^(tY_R). Its integrand relative to the source law is exp(tx−(1+b−v_R(x))t²/2). Outside the fixed set |x|≤2 the integrands coincide. On the full perturbation support ensure f_R=1+b−v_R≥1/2 for small R. For each fixed derivative order k, all t-derivatives through order k and their derivatives in v are uniformly bounded there by constants depending on b,k: their absolute values are bounded by a polynomial in |t| times exp(2|t|−t²/4). The core and transition have source mass O(R); on the exterior use the O(R) bound on the variance change. Consequently
    sup_(t real) |∂_t^k(F_R(t)−1)| ≤ C_k R.
The same argument on any fixed horizontal strip holds with another constant, since the extra real-part contribution from the imaginary coordinate is bounded on the strip. As v_R(0)=0 and v_R=b at x=R+R⁵, these smooth graphs have oscillation at least b and are not an equicontinuous family at zero.

For a growing sequence of finite node and order lists, apply this construction separately to list j and make R_j small enough that its first j normalized derivatives on |Im t|≤j have error at most 2^(-j). This diagonal selection has no uniform bound on the derivatives of v_R; it gives no continuous nonconstant exact limiting graph.
