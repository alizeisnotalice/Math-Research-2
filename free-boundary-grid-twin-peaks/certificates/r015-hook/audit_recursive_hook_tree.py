#!/usr/bin/env python3
"""Build the r015 recursive H-skeleton tree and check it exactly for n<=128."""

from math import factorial, log
import json
from pathlib import Path
import sys


def vid(n, x, y):
    return y * n + x


def add_edge(edges, n, a, b):
    assert abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1
    edges.add(tuple(sorted((vid(n, *a), vid(n, *b)))))


def build_block(n, edges, x0, y0, w, h):
    """Add a tree on a rectangle rooted at its top-left corner."""
    if w <= 2 or h <= 2:
        snake = []
        for yy in range(y0, y0 + h):
            xs = range(x0, x0 + w) if (yy - y0) % 2 == 0 else range(x0 + w - 1, x0 - 1, -1)
            snake.extend((xx, yy) for xx in xs)
        for a, b in zip(snake, snake[1:]):
            add_edge(edges, n, a, b)
        return

    c, d = w // 2, h // 2
    # H: left column, central column, and central row.  It is a tree.
    for y in range(y0, y0 + h - 1):
        add_edge(edges, n, (x0, y), (x0, y + 1))
        add_edge(edges, n, (x0 + c, y), (x0 + c, y + 1))
    for x in range(x0, x0 + w - 1):
        add_edge(edges, n, (x, y0 + d), (x + 1, y0 + d))

    children = [
        (x0 + 1, y0, c - 1, d, (x0, y0)),
        (x0 + 1, y0 + d + 1, c - 1, h - d - 1, (x0, y0 + d + 1)),
        (x0 + c + 1, y0, w - c - 1, d, (x0 + c, y0)),
        (x0 + c + 1, y0 + d + 1, w - c - 1, h - d - 1, (x0 + c, y0 + d + 1)),
    ]
    for cx, cy, cw, ch, attach in children:
        if cw <= 0 or ch <= 0:
            continue
        assert cw <= w // 2 and ch <= h // 2
        root = (cx, cy)
        add_edge(edges, n, attach, root)
        build_block(n, edges, cx, cy, cw, ch)


def check(n):
    edges = set()
    build_block(n, edges, 0, 0, n, n)
    N = n * n
    assert len(edges) == N - 1
    adj = [[] for _ in range(N)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    parent = [-2] * N
    parent[0] = -1
    stack = [0]
    order = []
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj[v]:
            if parent[u] == -2:
                parent[u] = v
                stack.append(u)
    assert len(order) == N
    sub = [1] * N
    for v in reversed(order[1:]):
        sub[parent[v]] += sub[v]
    hook_product = 1
    for x in sub:
        hook_product *= x
    assert hook_product <= 1 << (24 * N)
    labels = factorial(N) // hook_product
    assert labels * hook_product == factorial(N)
    log_hook = sum(log(x) for x in sub)
    assert log_hook <= 17 * N
    return {"n": n, "vertices": N, "edges": len(edges),
            "log_hook_over_N": log_hook / N,
            "hook_integer": labels > 0}


def main():
    rows = [check(n) for n in range(1, 129)]
    worst = max(rows, key=lambda x: x["log_hook_over_N"])
    ans = {"range": "n=1,...,128", "all_tree_connectivity_and_integer_hook_checks_passed": True,
           "worst_log_hook_over_N": worst, "status": "passed"}
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("recursive-hook-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
