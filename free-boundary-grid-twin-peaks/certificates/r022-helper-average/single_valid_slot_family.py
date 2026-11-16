#!/usr/bin/env python3
"""Full-grid family where only one old z-slot keeps exactly three peaks."""
import json


def witness(n):
    assert n >= 5
    # q,p,r are the retained roots; z suppresses the fourth root v.
    q, p, r, z, v = (4, 2), (1, 2), (2, 1), (2, 2), (2, 3)
    V = [(x, y) for y in range(n) for x in range(n)]
    neigh = {x: [(xx, yy) for xx, yy in
                 ((x[0]-1, x[1]), (x[0]+1, x[1]),
                  (x[0], x[1]-1), (x[0], x[1]+1))
             if 0 <= xx < n and 0 <= yy < n]
             for x in V}
    roots_base = [q, p, r, v]
    seen, queue, tail = set(roots_base), roots_base.copy(), []
    for x in queue:
        for y in neigh[x]:
            if y != z and y not in seen:
                seen.add(y); queue.append(y); tail.append(y)
    assert seen == set(V) - {z}
    base = roots_base + tail

    def peakset(word):
        earlier, out = set(), []
        for x in word:
            if not (set(neigh[x]) & earlier):
                out.append(x)
            earlier.add(x)
        return out

    assert peakset(base) == roots_base
    # The slots strictly after r in base.  Only slot 3 (before v) permits z
    # to suppress v before its birth; later slots leave v as a fourth peak.
    rpos = base.index(r)
    valid = []
    for slot in range(rpos + 1, len(base) + 1):
        source = base.copy(); source.insert(slot, z)
        if peakset(source) == [q, p, r]:
            valid.append(slot)
    assert valid == [3]
    source = base.copy(); source.insert(3, z)
    target = source.copy(); target.remove(z); target.insert(2, z)
    assert peakset(source) == [q, p, r]
    assert peakset(target) == [q, p]
    P, R, Z = map(source.index, (p, r, z))
    assert (P, R, Z) == (1, 2, 3)
    return {'n': n, 'base_peaks_after_deleting_z': roots_base,
            'source_peaks': [q, p, r], 'target_peaks': [q, p],
            'valid_old_z_slots_after_r': valid,
            'positions': {'p': P, 'r': R, 'z': Z},
            'weighted_gain_numerator': R-P,
            'weighted_gain_denominator': n*n-R-1}


print(json.dumps({'status': 'passed',
                  'witnesses': [witness(n) for n in range(5, 13)]}, indent=2))
