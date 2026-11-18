from collections import Counter
from functools import lru_cache
from random import Random
import platform


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                a, b = values[index], values[index + width]
                values[index], values[index + width] = a + b, a - b
        width *= 2
    return values


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            i = x.bit_length()
            if i not in pivots:
                pivots[i] = x
                break
            x ^= pivots[i]
    return len(pivots)


def anf(n, table):
    a = [(table >> point) & 1 for point in range(1 << n)]
    for i in range(n):
        for point in range(1 << n):
            if point >> i & 1:
                a[point] ^= a[point ^ (1 << i)]
    return tuple(point for point, coefficient in enumerate(a) if coefficient)


def all_controls(point, controls):
    return all((point >> i) & 1 for i in controls)


def run_layers(layers, point):
    for layer in reversed(layers):
        for target, terms in layer:
            value = sum(all_controls(point, controls) for controls in terms) & 1
            point ^= value << target
    return point & 1


def monomial_circuit(target, controls, scratch):
    if len(controls) <= 2:
        return [(target, controls)]
    anc = scratch[0]
    A = monomial_circuit(anc, controls[:-1], scratch[1:])
    B = [(target, (anc, controls[-1]))]
    return A + B + A + B


def verify_fixed_register(n, table):
    k = max(0, n - 2)
    x_offset = k + 1
    scratch = tuple(range(1, x_offset))
    gates = []
    for monomial in anf(n, table):
        controls = tuple(x_offset + j for j in range(n) if (monomial >> j) & 1)
        gates.extend(monomial_circuit(0, controls, scratch))
    assert all(len(controls) <= 2 and all(target < control for control in controls)
               for target, controls in gates)
    for initial in range(1 << (n + k + 1)):
        point = initial
        for target, controls in gates:
            point ^= all_controls(point, controls) << target
        expected = initial ^ ((table >> (initial >> x_offset)) & 1)
        assert point == expected
    return 1 << (n + k + 1)


def verify_fresh_addresses(n, table):
    masks = anf(n, table)
    height = (n - 1).bit_length()
    size = 1 << height
    layers = [[(0, tuple((j + 1,) for j in range(len(masks))))]]
    frontier = tuple(range(1, len(masks) + 1))
    cursor = 1 + len(masks)
    for _ in range(height):
        layer = []
        child = []
        for target in frontier:
            child.extend((cursor, cursor + 1))
            layer.append((target, ((cursor, cursor + 1),)))
            cursor += 2
        layers.append(layer)
        frontier = tuple(child)
    assert len(frontier) == len(masks) * size
    for x in range(1 << n):
        point = 0
        for j, mask in enumerate(masks):
            factors = [(x >> i) & 1 for i in range(n) if (mask >> i) & 1]
            factors.extend([1] * (size - len(factors)))
            for i, value in enumerate(factors):
                point |= value << frontier[j * size + i]
        assert run_layers(layers, point) == (table >> x) & 1
    return cursor


def derivative_barrier():
    pi = (6,7,12,11,2,0,5,4,13,15,8,14,1,9,10,3)
    size = 16
    assert set(pi) == set(range(size))
    for v in range(1, size):
        for u in range(1, size):
            assert any((u & (pi[y] ^ pi[y ^ v])).bit_count() % 2 == 0
                       for y in range(size))
        for w in range(v+1, size):
            assert any(pi[y] ^ pi[y ^ v] ^ pi[y ^ w] ^ pi[y ^ v ^ w]
                       for y in range(size))
    truth = [(x & pi[y]).bit_count() % 2 for y in range(16) for x in range(16)]
    raw = walsh(1 - 2 * bit for bit in truth)
    assert len(raw) == 256 and set(map(abs, raw)) == {16}
    checks = 0
    for v in range(1, 256):
        derivative = [value ^ truth[point ^ v] for point, value in enumerate(truth)]
        for w in range(v+1, 256):
            assert any(value == derivative[point ^ w] for point, value in enumerate(derivative))
            checks += 1
    return checks


def shift_skeleton(depth):
    cursor = 1
    leaves = [0]
    layers = []
    for _ in range(depth):
        children = []
        layer = []
        for target in leaves:
            A, l, r = cursor, cursor + 1, cursor + 2
            layer.append((target, ((l,), (A,l), (A,r))))
            children.extend((l,r))
            cursor += 3
        layers.append(layer)
        leaves = children
    return layers, leaves, cursor


def spectrum(layers, n):
    raw = walsh(1 - 2*run_layers(layers, point) for point in range(1 << n))
    support = [freq for freq, value in enumerate(raw) if value]
    assert all((abs(raw[freq]) & (abs(raw[freq])-1)) == 0 for freq in support)
    return len(support), rank((1 << n) | freq for freq in support), rank(support)


def verify_selector_padding():
    checks = []
    hard = None
    for h in range(3):
        layers, leaves, cursor = shift_skeleton(h)
        M = len(leaves)
        assert spectrum(layers, cursor) == (M*M, 2*M-1, 2*M-1)
        temp, X3 = cursor, cursor + 1
        X1, X2 = cursor + 2, cursor + 3
        layers += [[(leaves[0], ((temp,X3),))],[(temp, ((X1,X2),))]]
        actual = spectrum(layers, cursor+4)
        if h == 0:
            hard = actual
        s_h, R_h, _ = hard
        assert actual == (M*(M-1+s_h), 2*(M-1)+R_h, 2*(M-1)+R_h)
        for x in range(8):
            point = ((x & 1) << X1) | (((x >> 1) & 1) << X2) | (((x >> 2) & 1) << X3)
            assert run_layers(layers, point) == (x == 7)
        checks.append(actual)
    return checks


@lru_cache(None)
def signed_weight(x):
    x = abs(x)
    if x <= 1:
        return x
    if x % 2 == 0:
        return signed_weight(x // 2)
    return 1 + min(signed_weight((x-1)//2), signed_weight((x+1)//2))


def digit_example():
    n = 14
    rng = Random(85025)
    truth = [rng.randrange(2) for _ in range(1 << n)]
    if sum(truth) % 2 == 0:
        truth[0] ^= 1
    raw = walsh(1-2*bit for bit in truth)
    assert all(value % 4 == 2 for value in raw)
    return dict(sorted(Counter(signed_weight(value) for value in raw).items()))


if __name__ == "__main__":
    print("Python", platform.python_version())
    total = tables = 0
    rng = Random(85023)
    for n in range(1,5):
        tests = range(1 << (1 << n)) if n <= 3 else [rng.randrange(1 << 16) for _ in range(64)]
        for table in tests:
            total += verify_fixed_register(n, table)
            verify_fresh_addresses(n, table)
            tables += 1
    print("synthesized truth tables and exact dirty-register states:", tables, total)
    print("explicit flat-spectrum obstruction, directional pairs:", derivative_barrier())
    print("fresh-selector padding (s,R,d):", verify_selector_padding())
    print("14-bit random odd table signed weights:", digit_example())
