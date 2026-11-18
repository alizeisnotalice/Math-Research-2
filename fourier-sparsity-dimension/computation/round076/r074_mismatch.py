"""Exact verification of the uncertainty-sharp, simultaneously matched branches."""
import json
import platform


def walsh(a):
    w = a[:]
    h = 1
    while h < len(w):
        for off in range(0, len(w), 2*h):
            for i in range(off, off+h):
                x, y = w[i], w[i+h]
                w[i], w[i+h] = x+y, x-y
        h *= 2
    return w


rows = []
for n in range(6, 11):
    M = 1 << n
    for k in range(2, n-1):
        D = {x for x in range(M) if x % (1 << k) == 1}
        f = [1-2*(x == 0) for x in range(M)]
        g = [z * (1-2*(x in D)) for x, z in enumerate(f)]
        wf, wg = walsh(f), walsh(g)
        E = {i for i in range(M) if abs(wf[i]) != abs(wg[i])}
        assert E == {i for i in range(M) if wf[i] != wg[i]}
        assert E == set(range(1 << k)) and len(E) * len(D) == M
        assert 0 not in wf and 0 not in wg
        assert max(map(abs, wg)) < max(map(abs, wf))
        rows.append({"n": n, "k": k, "physical_disagreements": len(D),
                     "spectral_disagreements": len(E), "walsh_maxima": [max(map(abs, wf)), max(map(abs, wg))]})
print(json.dumps({"python": platform.python_version(), "cases": rows}, indent=2))
