# Return of distant near-maxima forced by slow logarithmic contact

Suppose E={a} is the only highest point in some fixed interval [a−r₀,a+r₀]. On its punctured neighborhood assume
$$
 A_-[\log(1/|u|)]^{-\beta}\le v(a+u)\le A_+[\log(1/|u|)]^{-\beta},
 \qquad 0<\beta<1,\quad 0<A_-\le A_+<\infty .
$$
No limiting coefficient is assumed. Define m_D=min{v(x):r₀≤|x−a|≤D+1}. For a continuous exact graph the following is necessary:
$$
 \limsup_{D\to\infty}D^{2\beta/(1-\beta)}m_D\le \frac{1}{2K},
 \qquad K=\frac{1-\beta}{\beta}\left(\frac{\beta}{2A_+}\right)^{1/(1-\beta)} .
$$
If the highest point is globally unique, the witnessing positions must escape to infinity. In particular a uniform positive exterior variance is unnecessary for exclusion: even exterior variances tending to zero more slowly than this power are insufficient. Other attained maxima outside the interval make the bound automatic.

For fixed D>0 evaluate the complex Gaussian output at a+iD with added variance ε>0 as in round 179. For purposes of proving the limsup it suffices to consider a sequence on which m_D≥cD^(-2β/(1−β))>0; on any such sequence the ε→0 passage is valid for each fixed D. Near a the absolute majorant is a constant times a power of log(1/|u|) times exp(C_D[log(1/|u|)]^β), integrable because β<1. The compact annulus away from a has positive minimum variance. Beyond distance D+1 the negative real exponent controls the prefactor uniformly even at v=0. A zero-age source there has off-diagonal limit zero. Thus the unsmoothed complex integral equals p_H(a+iD).

On |u|≤D^(-2), the cosine of Du/v(a+u) is ≥1−o(1) uniformly, since u(log(1/|u|))^β increases there in |u|. In this region retain e^(−L)≤u≤2e^(−L) with L=[βD²/(2A_+)]^(1/(1−β)). Its real contribution is at least
$$
 \exp\left[-L+(1-o(1))D^2 L^\beta/(2A_+)\right]
 =\exp[(K-o(1))D^{2/(1-\beta)}].
$$
Here the spatial input density stays bounded below and the prefactor cannot reduce the exponential rate. The remaining part of this inner region has nonnegative real contribution. For D^(-2)<|u|<r₀, the absolute contribution is at most exp(O(D²(log D)^β)) using the lower variance bound and compactness away from zero. For |u|>D+1 it is uniformly bounded. The target grows at most exp(O(D²)). The only remaining sources have absolute integral at most (2πm_D)^(-1/2)exp(D²/(2m_D)). If m_D D^{2β/(1−β)} exceeds 1/(2K) by any fixed positive amount along a sequence, these sources also have strictly insufficient exponential rate. This proves the bound.

## An exterior alternative for arithmetic Cantor clusters

For either geometry of round 186, its positive-minus-negative real Fourier lower bound for the compact source U actually holds for all β>0, before comparing with a strip. In the strict semicircle case one can use rate J(1,A_+), optimizing the length exp(−u k^α) of the single positive core; for the middle thirds case use J(1−d,A_+) under the stated strict coefficient ratio. Here α=2/(1+β), and J is defined in round 186. Thus the absolute Fourier contribution is at least exp(−(J+o(1))k^α) at the selected frequencies. For a strictly acute-phase set only a sequence of frequencies is required. For any K'>J, the exact full transform and the Gaussian input tail bound imply that some point outside U with |x|≤√(2K')k^(α/2) must satisfy v(x)≤2K' k^(α−2) for all sufficiently large selected k. Otherwise the entire exterior inside that radius is dominated by exp(−K' k^α), as is its Gaussian tail up to an absolute factor, contradicting cancellation of the U contribution. If E contains all maxima these exterior points escape. This also gives exclusion for every β>0 when v has a positive exterior lower bound. It quantifies the necessary tail behavior for β≤1 where a fixed strip estimate alone can be insufficient.

For β=1 there is an additional purely local conclusion whenever the distance D₀ from E to all other zeros exceeds 2J (use ∞ if no other zeros). The compact source density from round 186 extends holomorphically to every strip of half-width d₀<D₀/2. For any real center y either its distance to E or its distance to the exterior zeros is ≥D₀/2. In the former case use the included sources directly and in the latter use target minus omitted sources: the zero-variance points and the means of sufficiently small variances are spatially separated on compact sets; remote tails are separated as well. Each representation is holomorphic locally in the corresponding vertical strip segment. The two extensions agree on overlaps extending up from the real line. Gaussian decay at large real parts follows from boundedness of U. Thus the Fourier upper bound exp(−d₀k) contradicts the lower bound when J<d₀. For a unique compact zero cluster this uses no condition at infinity.
