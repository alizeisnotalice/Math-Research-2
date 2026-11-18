"""True Walsh checks of the sharp 3/2 balanced-support family."""
import itertools
import json
import sys

from verify_nonfull import eigenphase, gamma_uniform_fixed_base, spectrum, stats


def function(point, address_bits, extended):
    leaves = 1 << address_bits
    index = point & (leaves - 1)
    terminal_bit = (point >> (address_bits + index)) & 1
    quadratic = 0 if not extended else eigenphase(point >> (address_bits + leaves + 2 * index), 2)
    return (-1) ** (terminal_bit ^ quadratic)


results = []
for k in (1, 2):
    leaves = 1 << k
    seed_dim = k + leaves
    root_truth = tuple(function(point, k, False) for point in range(1 << seed_dim))
    assert stats(root_truth) == (leaves**2, leaves + k)
    cap, bases = gamma_uniform_fixed_base(spectrum(root_truth), 1 << seed_dim)
    assert cap == 1
    parent = tuple(function(point, k, True) for point in range(1 << (seed_dim + 2 * leaves)))
    assert stats(parent) == (4 * leaves**2, 3 * leaves + k)
    for offset in range(1 << (2 * leaves)):
        slice_values = parent[offset << seed_dim:(offset + 1) << seed_dim]
        assert stats(slice_values) == (leaves**2, leaves + k)
    results.append(dict(k=k, uniform_basis_margin_mult=str(cap), bases=bases,
                        true_addresses=1 << (2 * leaves), seed_s_R=stats(root_truth), parent_s_R=stats(parent)))
print(json.dumps({"python":sys.version.split()[0], "cases":results}, indent=2))
