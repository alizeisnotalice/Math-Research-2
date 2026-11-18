"""Whole-box spindle lower forms, rational LP recombination and coordinate pruning."""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from itertools import product
import ctypes
import os
import sys
import warnings
import numpy as np
from flint import arb, ctx
from scipy.optimize import linprog, OptimizeWarning

ctx.prec = 144
P = 2**40
S = 2**44
goal = int(2**57 * Fraction(".836490455462")) + 1
warnings.filterwarnings("ignore", category=OptimizeWarning)
OPTIONS = dict(dual_feasibility_tolerance=1e-10, primal_feasibility_tolerance=1e-10)


@lru_cache(None)
def proposals():
    sys.path.insert(0, os.environ["FOUR_SPINDLE_GEOMETRY"])
    from gen8 import hull_intervals

    return hull_intervals


def hull_intervals(mid):
    return proposals()(mid)


@lru_cache(None)
def basis_moves():
    corners = np.array(list(product((-1, 1), repeat=8)), dtype=np.int32)
    rng = np.random.default_rng(20261118)
    eye = np.eye(8, dtype=np.int32)
    return np.r_[
        np.zeros((1, 8), dtype=np.int32),
        eye,
        -eye,
        corners[rng.permutation(len(corners))],
    ]


@lru_cache(None)
def base():
    dll = ctypes.CDLL(os.environ["FOUR_SPINDLE_LIBRARY"])
    v = lambda t: np.ctypeslib.ndpointer(dtype=t, ndim=1, flags="C_CONTIGUOUS")
    dll.forms_for_box.argtypes = (
        [v(np.int64)] * 4 + [v(np.int32)] * 2 + [ctypes.c_int, v(np.int64)]
    )
    dll.box_try.argtypes = [v(np.int64)] * 4 + [v(np.int32), ctypes.c_int64]
    dll.rotate_material.argtypes = [
        v(np.int64),
        ctypes.c_int64,
        ctypes.c_int64,
        ctypes.c_int,
        v(np.int64),
    ]
    t = np.load(Path(__file__).with_name("four-body-spindle-material.npz"))
    q = t["points"]
    cnt = t["counts"]
    assert q.dtype == np.int64
    off = np.r_[0, np.cumsum(cnt)]
    blocks = [q[off[k] : off[k + 1]].ravel() for k in range(4)]
    return dll, blocks, cnt, t["original_hi"], basis_moves()


@lru_cache(50000)
def rotated(b, total):
    dll, bs, cnt, hi, move = base()
    th = arb(total) / (2 * S)
    s, c = th.sin_cos()
    C = round(float(c.mid()) * P)
    sin_n = round(float(s.mid()) * P)
    assert abs(c * P - C) < 1 and abs(s * P - sin_n) < 1
    q = bs[b]
    out = np.empty_like(q)
    dll.rotate_material(q, C, sin_n, len(q) // 2, out)
    return out


@lru_cache(40000)
def vertices(n, total):
    th = arb(total) / (2 * S)
    R = 1 / (2 * (arb.pi() / (2 * n)).cos())
    v = []
    for j in range(n):
        s, c = (th + 2 * arb.pi() * j / n).sin_cos()
        x, y = R * c * P, R * s * P
        i, jj = round(float(x.mid())), round(float(y.mid()))
        assert abs(x - i) < 1 and abs(y - jj) < 1
        v.append((i, jj))
    return np.asarray(v, dtype=np.int64)


def materials(l, h):
    mid = (l[:8] + h[:8]) / (2 * S)
    seq = hull_intervals(mid)
    cts = np.array(
        [
            a
            for prev, (a, b, v) in zip(seq[-1:] + seq[:-1], seq)
            if prev[2][5:] != v[5:] or prev[2][4] != v[4]
        ]
    )
    c, s = np.cos(cts), np.sin(cts)
    dll, static, cnt, hi, move = base()
    parts = []
    counts = cnt.copy()
    for j, n in enumerate((0, 3, 5, 7)):
        i = 3 * j - 4
        total = int(l[i]) + int(h[i]) if j >= 2 else 0
        w = total / (2 * S)
        if n:
            R = 0.5 / np.cos(np.pi / (2 * n))
            ang = (cts - w) % (2 * np.pi)
            iv = np.rint(ang * n / (2 * np.pi)).astype(int) % n
            delta = (ang + np.pi / n) % (2 * np.pi / n) - np.pi / n
            ic = np.rint((ang + np.pi) * n / (2 * np.pi)).astype(int) % n
            corner = np.abs(delta) <= np.pi / (2 * n)
            x = np.where(
                corner,
                R * np.cos(w + 2 * np.pi * iv / n),
                c + R * np.cos(w + 2 * np.pi * ic / n),
            )
            y = np.where(
                corner,
                R * np.sin(w + 2 * np.pi * iv / n),
                s + R * np.sin(w + 2 * np.pi * ic / n),
            )
        else:
            x, y = c / 2, s / 2
        cand = np.unique(
            np.rint(np.column_stack((x, y)) * (P - 20)).astype(np.int64), axis=0
        )
        if n:
            vv = vertices(n, total).astype(object)
            delta = cand.astype(object)[:, None, :] - vv[None, :, :]
            inside = np.all(np.sum(delta * delta, axis=2) < (P - 2) ** 2, axis=1)
        else:
            coords = cand.astype(object)
            inside = np.sum(coords * coords, axis=1) < (P // 2) ** 2
        new = cand[inside].ravel()
        old = rotated(j, total) if j >= 2 else static[j]
        parts.extend([old, new])
        counts[j] += len(new) // 2
    return dll, np.concatenate(parts), counts, hi, move


def dense_forms(l, h, m):
    dll, p, cnt, hi, move = materials(l, h)
    axis = dll.box_try(l, h, hi, p, cnt, goal)
    if axis < 0:
        return None, axis
    arr = np.zeros((m, 9), dtype=np.int64)
    v = dll.forms_for_box(l, h, hi, p, cnt, move[:m].ravel(), m, arr.ravel())
    assert v >= 0
    return np.unique(arr[:v], axis=0), axis


from math import ceil, cos, sin


def owned(l, h):
    mid = (l[:8] + h[:8]) / (2 * S)
    seq = hull_intervals(mid)
    runs = []
    for a, b, v in seq:
        if runs and runs[-1][2] == v:
            runs[-1][1] = b
        else:
            runs.append([a, b, v])
    dll, bs, cnt, hi, mv = base()
    parts = []
    counts = np.zeros(5, dtype=np.int32)
    for j, n in enumerate([0, 3, 5, 7]):
        idx = 3 * j - 3
        tx, ty = mid[idx : idx + 2] if n else (0, 0)
        ca = []
        for a, b, v in runs:
            if v[5] != j:
                continue
            cx, cy, r = v[2:5]
            cx -= tx
            cy -= ty
            for u in (
                np.linspace(a, b, 1 + max(1, ceil((b - a) / 0.035))) if r else (a,)
            ):
                ca.append(
                    (
                        round((cx + r * cos(u)) * (P - 20)),
                        round((cy + r * sin(u)) * (P - 20)),
                    )
                )
        cand = np.unique(np.array(ca, dtype=np.int64).reshape(-1, 2), axis=0)
        if n:
            i = 3 * j - 4
            total = int(l[i]) + int(h[i]) if j >= 2 else 0
            vv = vertices(n, total).astype(object)
            diff = cand.astype(object)[:, None, :] - vv
            assert np.all(np.sum(diff * diff, axis=2) < (P - 2) ** 2)
        else:
            vv = cand.astype(object)
            assert np.all(np.sum(vv * vv, axis=1) < (P // 2) ** 2)
        ca = np.r_[cand, np.array([[0, 0], [1, 0], [0, 1]], np.int64)]
        counts[j] = len(ca)
        parts.append(ca.ravel())
    return dll, np.concatenate(parts), counts, hi, mv


def formset(l, h, m):
    dll, p, cnt, hi, mv = owned(l, h)
    axis = dll.box_try(l, h, hi, p, cnt, goal)
    if axis < 0:
        return None, axis
    # Additional support points provide guide alternatives for multi-contact centers.
    a, axis2 = dense_forms(l, h, m)
    if a is None:
        return None, axis2
    arr = np.zeros((m, 9), dtype=np.int64)
    v = dll.forms_for_box(l, h, hi, p, cnt, mv[:m].ravel(), m, arr.ravel())
    assert v >= 0
    return np.unique(np.r_[arr[:v], a], axis=0), axis


def combined(forms, vals):
    ids = np.flatnonzero(vals > 1e-11)
    if not len(ids):
        return None, 0
    v = vals[ids]
    v = v / v.max() * 2**31
    ws = [max(0, round(z)) for z in v]
    return forms[ids].T.astype(object) @ np.asarray(ws, dtype=object), sum(ws)


def contract(l, h, m):
    a, axis = formset(l, h, m)
    if a is None:
        return 1, axis, l, h, 0
    if not len(a):
        return 0, axis, l, h, 0
    g = a[:, 1:].astype(float) / 2**52
    f = (a[:, 0] - goal).astype(float) / 2**52
    minmax = linprog(
        [0] * 8 + [1],
        A_ub=np.c_[g, -np.ones(len(f))],
        b_ub=-f,
        bounds=[(-1, 1)] * 8 + [(None, None)],
        method="highs",
        options=OPTIONS,
    )
    assert minmax.success
    z, den = combined(a, -minmax.ineqlin.marginals)
    if den and z[0] - sum(map(abs, z[1:])) >= den * goal:
        return 2, axis, l, h, 0
    if minmax.fun > 1e-9:
        return 0, axis, l, h, 0
    newl = l.copy()
    newh = h.copy()
    np_cut = 0
    for j in range(8):
        for sign in (1, -1):
            d = np.zeros(8)
            d[j] = -sign
            lp = linprog(
                d,
                A_ub=g,
                b_ub=-f,
                bounds=[(-1, 1)] * 8,
                method="highs",
                options=OPTIONS,
            )
            if not lp.success or sign * lp.x[j] > 0.95:
                continue
            z, den = combined(a, -lp.ineqlin.marginals)
            if not den:
                continue
            Gj = z[j + 1]
            D = sign * Gj
            if D <= 0:
                continue
            N = (
                goal * den
                - z[0]
                + sum(abs(val) for i, val in enumerate(z[1:]) if i != j)
            )
            a0 = int(l[j])
            b0 = int(h[j])
            num = (a0 + b0) * D + sign * (b0 - a0) * N
            if sign == 1:
                cut = -((-num) // (2 * D))
                if cut < b0:
                    newh[j] = max(a0 - 1, int(cut))
                    np_cut += 1
            else:
                cut = num // (2 * D)
                if cut > a0:
                    newl[j] = min(b0 + 1, int(cut))
                    np_cut += 1
    if (newh[:8] < newl[:8]).any():
        return 3, axis, newl, newh, np_cut
    return 0, axis, newl, newh, np_cut


def work(row):
    i, label, v, limit, depth_limit = row
    l = np.zeros(11, dtype=np.int64)
    r = l.copy()
    l[:8] = v[:8]
    r[:8] = v[8:]
    stack = [(l, r, 0)]
    vis = 0
    succ = Counter()
    left = []
    splits = 0
    md = 0
    contractions = 0
    while stack and vis < limit:
        a, b, k = stack.pop()
        vis += 1
        md = max(md, k)
        ok, axis, L, H, cut = contract(a, b, 65)
        if ok:
            succ[ok] += 1
            continue
        wid = b[:8] - a[:8]
        change = (wid - (H[:8] - L[:8])) * np.array([1, 1, 0.16, 1, 1, 0.11, 1, 1])
        eff = wid * np.array([1, 1, 0.16, 1, 1, 0.11, 1, 1])
        ratio = 1 - change.max() / max(1, eff.max())
        if k >= depth_limit:
            left.append(np.r_[L[:8], H[:8]].tolist())
            continue
        if ratio < 0.90:
            contractions += 1
            stack.append((L, H, k + 1))
            continue
        eff = (H[:8] - L[:8]) * np.array([1, 1, 0.16, 1, 1, 0.11, 1, 1])
        if eff[axis] * 24 < eff.max() or H[axis] - L[axis] < 2:
            axis = int(np.argmax(eff))
        if H[axis] - L[axis] < 2:
            left.append(np.r_[L[:8], H[:8]].tolist())
            continue
        mid = L[axis] + (H[axis] - L[axis]) // 2
        C = L.copy()
        C[axis] = mid
        B = H.copy()
        B[axis] = mid
        splits += 1
        stack.extend([(C, H, k + 1), (L, B, k + 1)])
    left.extend(np.r_[a[:8], b[:8]].tolist() for a, b, k in stack)
    return dict(
        id=i,
        label=label,
        vis=vis,
        splits=splits,
        contractions=contractions,
        closed=dict(succ),
        remaining=left,
        md=md,
    )
