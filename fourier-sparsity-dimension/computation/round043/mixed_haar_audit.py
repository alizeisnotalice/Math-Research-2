#!/usr/bin/env python3
"""Exact conditional Walsh and weighted Haar checks for mixed selector trees."""

import argparse
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
import json
from math import lcm
from pathlib import Path
import platform
import time


Y = 4
LEAF_MASKS = (0, 15, 3, 1, 2)


def private(left, right):
    return ("private", None, left, right)


def shared(predicate, left, right):
    return ("shared", predicate, left, right)


def parity_mask(normal):
    return sum((((normal & y).bit_count() & 1) << y) for y in range(Y))


def fwht(values):
    out = list(values)
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2 * step):
            for i in range(start, start + step):
                j = i + step
                a, b = out[i], out[j]
                out[i], out[j] = a + b, a - b
        step *= 2
    return out


def gf2_rank(values):
    basis = {}
    for value in values:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


@dataclass
class Node:
    kind: str
    bit: int
    left: int
    right: int
    y_path: int
    private_path: tuple


class Compiled:
    def __init__(self, shape):
        self.nodes = []
        self.private_count = 0
        self.leaf_count = 0
        self.compile(shape, (1 << Y) - 1, ())
        self.Z = 1 << self.private_count
        self.routes = [
            [self.terminal(v, point % Y, point // Y) for point in range(self.Z * Y)]
            for v in range(len(self.nodes))
        ]
        self.descendants = {}
        self.fibres = {}
        for v, node in enumerate(self.nodes):
            down = [
                u for u in self.below(v)
                if self.nodes[u].kind == "private" and self.nodes[u].y_path
            ]
            self.descendants[v] = down
            prefix = len(node.private_path)
            for u in down:
                descendant = self.nodes[u]
                path = descendant.private_path[prefix:]
                shifts = [(1 << descendant.bit, 1)]
                for bit, branch in path:
                    shifts += [(k ^ (1 << bit), s * (1 - 2 * branch)) for k, s in shifts]
                self.fibres[v, u] = shifts

    def compile(self, shape, y_path, private_path):
        v = len(self.nodes)
        if shape is None:
            self.nodes.append(Node("leaf", self.leaf_count, -1, -1, y_path, private_path))
            self.leaf_count += 1
            return v
        kind, pred, left, right = shape
        if kind == "private":
            bit = self.private_count
            self.private_count += 1
        else:
            bit = pred
        node = Node(kind, bit, -1, -1, y_path, private_path)
        self.nodes.append(node)
        if kind == "private":
            node.left = self.compile(left, y_path, private_path + ((bit, 0),))
            node.right = self.compile(right, y_path, private_path + ((bit, 1),))
        else:
            node.left = self.compile(left, y_path & ~pred, private_path)
            node.right = self.compile(right, y_path & pred, private_path)
        return v

    def terminal(self, v, y, z):
        while self.nodes[v].kind != "leaf":
            node = self.nodes[v]
            address = z if node.kind == "private" else node.bit
            bit = node.bit if node.kind == "private" else y
            v = node.right if ((address >> bit) & 1) else node.left
        return self.nodes[v].bit

    def below(self, v):
        yield v
        node = self.nodes[v]
        if node.kind != "leaf":
            yield from self.below(node.left)
            yield from self.below(node.right)


def expected_means(tree, values):
    """Each recursive conditional mean is scaled by the total private domain size."""
    nodes = tree.nodes
    means = [[] for _ in nodes]
    profiles = {}
    for v in reversed(range(len(nodes))):
        node = nodes[v]
        path = node.y_path
        if node.kind == "leaf":
            means[v] = [tree.Z * val for val in values[node.bit]]
        elif node.kind == "shared":
            left, right = means[node.left], means[node.right]
            means[v] = [right[y] if (node.bit >> y) & 1 else left[y] for y in range(Y)]
        else:
            left, right = means[node.left], means[node.right]
            assert all((a - b) % 2 == 0 for a, b in zip(left, right))
            means[v] = [(a + b) // 2 for a, b in zip(left, right)]
            profile = fwht((left[y] - right[y]) // 2 if (path >> y) & 1 else 0 for y in range(Y))
            if any(profile):
                profiles[v] = profile
                assert sum(w != 0 for w in profile) * path.bit_count() >= Y
    return means, profiles


def audit_case(tree, leaf_masks, stats):
    values = [[1 - 2 * ((mask >> y) & 1) for y in range(Y)] for mask in leaf_masks]
    means, profiles = expected_means(tree, values)
    denominator = lcm(*range(1, Y + 1))
    active_counts = [0] * len(tree.nodes)
    shared_fork_counts = [0] * len(tree.nodes)
    scaled_cost = [0] * len(tree.nodes)
    for v in reversed(range(len(tree.nodes))):
        node = tree.nodes[v]
        path = node.y_path
        if not path:
            stats["dead_path_checks"] += 1
        raw = [
            values[leaf][x % Y] if (path >> (x % Y)) & 1 else 0
            for x, leaf in enumerate(tree.routes[v])
        ]
        actual = fwht(raw)
        stats["compared_coefficients"] += len(actual)
        assembled = fwht(means[v][y] if (path >> y) & 1 else 0 for y in range(Y))
        assembled += [0] * (len(raw) - Y)
        active_bits, cost = 0, 0
        for u in tree.descendants[v]:
            if u not in profiles:
                continue
            fibres = tree.fibres[v, u]
            cost += len(fibres) * (denominator // tree.nodes[u].y_path.bit_count())
            for k, sign in fibres:
                active_bits |= k
                for beta, coefficient in enumerate(profiles[u]):
                    index = k * Y + beta
                    assert assembled[index] == 0 and coefficient % len(fibres) == 0
                    assembled[index] = sign * (coefficient // len(fibres))
        assert actual == assembled
        private_frequencies = {i // Y for i, entry in enumerate(actual) if entry}
        support_size = sum(entry != 0 for entry in actual)
        n = gf2_rank(private_frequencies)
        active_counts[v] = n
        scaled_cost[v] = cost
        g = 0
        if node.kind != "leaf":
            left_n, right_n = active_counts[node.left], active_counts[node.right]
            g = shared_fork_counts[node.left] + shared_fork_counts[node.right]
            if node.kind == "shared" and left_n and right_n:
                g += 1
                stats["shared_double_effective_forks"] += 1
        shared_fork_counts[v] = g
        assert n == active_bits.bit_count()
        assert all(
            bool((active_bits >> bit) & 1)
            == any(raw[x] != raw[x ^ (Y << bit)] for x in range(len(raw)))
            for bit in range(tree.private_count)
        )
        assert support_size * denominator >= cost * Y
        if n:
            lower = denominator * (n + 1) ** 2
            upper = 4 * path.bit_count() * cost
            assert upper >= lower
            stats["positive_subtree_budget_checks"] += 1
            stats["tight_subtree_budgets"] += upper == lower
            stronger = denominator * (n + g + 1) ** 2
            assert upper >= stronger
            stats["tight_strengthened_subtree_budgets"] += upper == stronger
            stats["positive_subtree_budget_checks_with_shared_forks"] += g > 0
        if v == 0:
            assert 4 * support_size >= (n + g + 1) ** 2
            stats["root_rank_budget_checks"] += 1
            stats["tight_root_rank_budgets"] += 4 * support_size == (n + 1) ** 2
            stats["tight_strengthened_root_rank_budgets"] += 4 * support_size == (n + g + 1) ** 2
        if node.kind == "leaf":
            continue
        left_c, right_c = scaled_cost[node.left], scaled_cost[node.right]
        if node.kind == "shared":
            assert n == left_n + right_n and cost == left_c + right_c
        elif path:
            mark = v in profiles
            assert cost == 2 * (left_c + right_c) + mark * (
                denominator // path.bit_count()
            )
            assert n == left_n + right_n + bool(mark or left_n or right_n)
            if bool(left_n) != bool(right_n):
                stats["single_effective_child_checks"] += 1
                assert mark
            if left_n and right_n and not mark:
                stats["equal_means_two_effective_children_checks"] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()

    L = None
    qx, qy, qxy = map(parity_mask, (1, 2, 3))
    and_mask = 1 << 3
    templates = {
        "all_private": private(private(L, L), private(L, L)),
        "root_with_one_private_child": private(private(L, L), L),
        "repeated_parity_dead_paths": shared(
            qxy, shared(qxy, private(L, L), private(L, L)), private(L, L)
        ),
        "same_parity_through_private_branches": private(
            shared(qxy, L, L), private(shared(qxy, L, L), L)
        ),
        "different_shared_parities": private(
            shared(qx, private(L, L), L), shared(qy, L, private(L, L))
        ),
        "unbalanced_shared_predicate": shared(
            and_mask, private(private(L, L), L), private(L, L)
        ),
        "conditional_unbalanced_parity": shared(
            and_mask, shared(qx, L, private(L, L)), private(L, L)
        ),
    }
    output = {}
    total = Counter()
    for name, shape in templates.items():
        tree = Compiled(shape)
        stats = Counter()
        for family in product(LEAF_MASKS, repeat=tree.leaf_count):
            audit_case(tree, family, stats)
            stats["leaf_families"] += 1
        cuts = [
            (tree.nodes[node.left].y_path.bit_count(), tree.nodes[node.right].y_path.bit_count())
            for node in tree.nodes if node.kind == "shared"
        ]
        output[name] = {
            "tree": shape,
            "private_coordinates": tree.private_count,
            "shared_child_path_counts_out_of_4": cuts,
            "counts": dict(stats),
        }
        total.update(stats)
        print(name, stats["leaf_families"], flush=True)
    assert total["dead_path_checks"] and total["equal_means_two_effective_children_checks"]
    assert total["single_effective_child_checks"] and total["shared_double_effective_forks"]
    result = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_variables": 2,
        "leaf_truth_masks": [hex(mask) for mask in LEAF_MASKS],
        "templates": output,
        "totals": dict(total),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
