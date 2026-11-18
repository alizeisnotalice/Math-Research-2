import json
import platform

import sympy as sp


def check_nonzero_zero():
    length, center = sp.symbols("L A", positive=True)
    identity = length + center**2 / length - 2 * center - (length - center)**2 / length
    assert sp.simplify(identity) == 0
    return True


def check_origin_zero():
    log_t, log_log_t, coefficient = sp.symbols("H J c", positive=True)
    scaled_coordinate = sp.Symbol("b", real=True)
    quadratic = log_t * scaled_coordinate - coefficient * (log_t + log_log_t) * scaled_coordinate**2 / 2
    center = log_t / (coefficient * (log_t + log_log_t))
    top = log_t**2 / (2 * coefficient * (log_t + log_log_t))
    assert sp.simplify(sp.diff(quadratic, scaled_coordinate).subs(scaled_coordinate, center)) == 0
    assert sp.simplify(quadratic.subs(scaled_coordinate, center) - top) == 0
    assert sp.simplify(
        top - (log_t - log_log_t) / (2 * coefficient)
        - log_log_t**2 / (2 * coefficient * (log_t + log_log_t))
    ) == 0
    return {"center": str(center), "maximum": str(top)}


def check_flat_density():
    relative_scale, slope, z = sp.symbols("s r z", real=True)
    gaussian_scale = 1 + slope * relative_scale * z
    integrand = sp.exp(z**2 * (1 - gaussian_scale**-2) / 2) / (
        gaussian_scale * (1 + relative_scale * z)
    )
    expansion = sp.series(integrand, relative_scale, 0, 3).removeO().expand()
    gaussian_moments = {z**2: 1, z**4: 3, z**6: 15}
    even_expansion = ((expansion + expansion.subs(z, -z)) / 2).expand()
    averaged_expansion = even_expansion.xreplace(gaussian_moments)
    expected = 1 + (1 - slope)**2 * relative_scale**2
    assert sp.simplify(averaged_expansion - expected) == 0
    return str(expected)


print(json.dumps({
    "python": platform.python_version(),
    "sympy": sp.__version__,
    "nonzero_zero_square_identity": check_nonzero_zero(),
    "origin_zero_quadratic": check_origin_zero(),
    "flat_density_gaussian_polynomial": check_flat_density(),
}, indent=2))
