#!/usr/bin/env python3
"""Portable archive rerun of the all-labelled-tree cut-formula scan.

This is the archived adaptation of ``exhaustive_cut_formula.py``.  The sole
semantic change is replacing its private absolute import path by the retained
sibling ``r024-global-normalizer/exact_kernel.py``.
"""
from itertools import product
from math import factorial, comb
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


t = load("frozen", HERE / "tree_peak_dp_r024_frozen.py")
e = load("subset", HERE.parent / "r024-global-normalizer" / "exact_kernel.py")


def prufer(code):
    n = len(code) + 2
    deg = [1] * n
    for x in code:
        deg[x] += 1
    adjacency = [set() for _ in range(n)]
    for x in code:
        u = next(i for i, d in enumerate(deg) if d == 1)
        adjacency[u].add(x)
        adjacency[x].add(u)
        deg[u] -= 1
        deg[x] -= 1
    u = [i for i, d in enumerate(deg) if d == 1]
    adjacency[u[0]].add(u[1])
    adjacency[u[1]].add(u[0])
    return adjacency


def rooted(adjacency, vertices, root):
    def size(v, parent):
        return 1 + sum(size(w, v) for w in adjacency[v]
                       if w != parent and w in vertices)

    product_hooks = 1

    def visit(v, parent):
        nonlocal product_hooks
        for w in adjacency[v]:
            if w != parent and w in vertices:
                visit(w, v)
                product_hooks *= size(w, v)

    visit(root, -1)
    return factorial(len(vertices) - 1) // product_hooks


def cut_sum(adjacency):
    n = len(adjacency)
    answer = 0
    for u in range(n):
        for v in adjacency[u]:
            if u >= v:
                continue
            left, stack = {u}, [u]
            while stack:
                x = stack.pop()
                for w in adjacency[x]:
                    if {x, w} == {u, v}:
                        continue
                    if w not in left:
                        left.add(w)
                        stack.append(w)
            right = set(range(n)) - left
            answer += comb(n - 2, len(left) - 1) * (
                sum(rooted(adjacency, left, x) for x in left)
                * sum(rooted(adjacency, right, x) for x in right)
                - rooted(adjacency, left, u) * rooted(adjacency, right, v)
            )
    return answer


for n in range(3, 7):
    completed = 0
    for code in product(range(n), repeat=n - 2):
        adjacency = prufer(code)
        _, by_cut = e.exact_counts([sum(1 << w for w in row)
                                    for row in adjacency])
        assert by_cut == cut_sum(adjacency)
        completed += 1
    print(n, completed, "passed")
