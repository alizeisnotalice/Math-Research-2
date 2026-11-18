"""Exact Walsh and joint-gate checks for the equality family, using integer arithmetic."""
from fractions import Fraction
import json
import platform


def walsh(values):
    out = list(values)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            for pos in range(start, start + width):
                a, b = out[pos], out[pos + width]
                out[pos], out[pos + width] = a + b, a - b
        width *= 2
    return out


def check(width):
    addresses = 1 << width
    values = [(-1) ** (x + w.bit_count()) for w in range(addresses) for x in (0, 1)]
    spectrum = walsh(values)
    assert {i for i, value in enumerate(spectrum) if value} == {2 * addresses - 1}
    gate_masses = []
    weighted_terms = []
    for a in range(addresses):
        mask = [int(w == a) for w in range(addresses)]
        assert sum(mask) == 1 and all(abs(z) == 1 for z in walsh(mask))
        h = [(-1) ** a.bit_count() * z for z in mask]
        assert sum(z != 0 for z in walsh(h)) == addresses
        gate_masses.append(Fraction(sum(mask), addresses))
        weighted_terms.append(Fraction(addresses, 4 ** width))
    assert sum(gate_masses) == sum(weighted_terms) == 1
    private_count = addresses
    branching_shared = addresses - 1
    assert private_count + branching_shared + 1 == 2 ** (width + 1)
    return dict(m=width, private_occurrences=private_count, branching_shared=branching_shared,
                gate_rank=width, nonzero_walsh_frequency=2 * addresses - 1,
                sparsity=1, weighted_sum=1)


print(json.dumps(dict(python=platform.python_version(), cases=[check(m) for m in range(9)]),
                 indent=2))
