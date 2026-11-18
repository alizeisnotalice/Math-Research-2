"""Exact checks of the shared/private selector-layer ranks and weighted costs."""

from fractions import Fraction
import json
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "round042"))
from audit_independent_tree_layers import rank, walsh


def private(left, right):
    return "private", None, left, right


def shared(mask, left, right):
    return "shared", mask, left, right


def assign(tree, counter):
    if isinstance(tree, int):
        return tree
    kind, index, left, right = tree
    if kind == "private":
        index = counter[0]
        counter[0] += 1
    return kind, index, assign(left, counter), assign(right, counter)


def parity(value):
    return value.bit_count() % 2


def evaluate(tree, x, w):
    while not isinstance(tree, int):
        kind, index, left, right = tree
        go_right = (x >> index) & 1 if kind == "private" else parity(w & index)
        tree = right if go_right else left
    return 1 - 2 * ((tree >> w) & 1)


def outcomes(tree, w):
    if isinstance(tree, int):
        return {evaluate(tree, 0, w)}
    kind, index, left, right = tree
    if kind == "shared":
        return outcomes(right if parity(w & index) else left, w)
    return outcomes(left, w) | outcomes(right, w)


def average(tree, w, scale):
    if isinstance(tree, int):
        return scale * evaluate(tree, 0, w)
    kind, index, left, right = tree
    if kind == "shared":
        return average(right if parity(w & index) else left, w, scale)
    numerator = average(left, w, scale) + average(right, w, scale)
    assert numerator % 2 == 0
    return numerator // 2


def profile(tree, reachable, data_size, scale):
    if not reachable or isinstance(tree, int):
        return 0, 0, Fraction(0), 0
    kind, index, left, right = tree
    probability = Fraction(len(reachable), data_size)
    if kind == "shared":
        reachable_left = {w for w in reachable if not parity(w & index)}
        reachable_right = reachable - reachable_left
        profiles = [profile(left, reachable_left, data_size, scale),
                    profile(right, reachable_right, data_size, scale)]
        nodes, first_private, cost, weighted_sparsity = tuple(map(sum, zip(*profiles)))
    else:
        n_left, _, c_left, s_left = profile(left, reachable, data_size, scale)
        n_right, _, c_right, s_right = profile(right, reachable, data_size, scale)
        values = []
        for w in range(data_size):
            if w not in reachable:
                values.append(0)
                continue
            difference = average(left, w, scale) - average(right, w, scale)
            assert difference % 2 == 0
            values.append(difference // 2)
        layer_sparsity = sum(value != 0 for value in walsh(values))
        if layer_sparsity:
            assert layer_sparsity * probability >= 1
        effective = any(len(outcomes(tree, w)) > 1 for w in reachable)
        first_private = int(effective)
        nodes = n_left + n_right + first_private
        cost = 2 * (c_left + c_right)
        if layer_sparsity:
            cost += 1 / probability
        weighted_sparsity = layer_sparsity + 2 * (s_left + s_right)
    assert 4 * probability * cost >= (nodes + first_private) ** 2
    return nodes, first_private, cost, weighted_sparsity


def check(tree, data_dimension):
    counter = [0]
    tree = assign(tree, counter)
    selector_size = 1 << counter[0]
    data_size = 1 << data_dimension
    spectrum = walsh(
        evaluate(tree, x, w) for w in range(data_size) for x in range(selector_size)
    )
    support = [i for i, value in enumerate(spectrum) if value]
    average_spectrum = walsh(average(tree, w, selector_size) for w in range(data_size))
    mean_sparsity = sum(value != 0 for value in average_spectrum)
    nodes, first_private, cost, layers = profile(
        tree, set(range(data_size)), data_size, selector_size
    )
    assert len(support) == layers + mean_sparsity
    assert nodes == rank(index & (selector_size - 1) for index in support)
    assert cost <= len(support)
    return {"all_private_nodes": counter[0], "effective_private_nodes": nodes,
            "first_private_nodes": first_private, "weighted_cost": str(cost),
            "sparsity": len(support), "bound_attained": (nodes + first_private) ** 2 == 4 * len(support)}


def main():
    plus, minus = 0, 15
    balanced = private(private(plus, minus), private(plus, minus))
    parallel = shared(1, shared(2, balanced, balanced), shared(2, balanced, balanced))
    interleaved = shared(
        1,
        private(shared(1, plus, private(plus, minus)),
                shared(2, minus, private(minus, plus))),
        private(shared(3, plus, minus), shared(2, private(minus, plus), plus)),
    )
    conditional = private(shared(2, plus, plus), shared(2, plus, minus))
    result = {
        "python": platform.python_version(),
        "shared_prefix_sharp_example": check(parallel, 2),
        "correlated_repeated_queries_and_unreachable_node": check(interleaved, 2),
        "marked_signal_support_strictly_smaller_than_reachable_path": check(conditional, 2),
    }
    assert result["shared_prefix_sharp_example"]["bound_attained"]
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
