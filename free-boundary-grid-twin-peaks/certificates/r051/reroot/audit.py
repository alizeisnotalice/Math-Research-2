"""Exhaustively test the reversible maps on induced branches of the 3-by-3 grid."""
from collections import Counter, deque
from math import comb, perm

VERTICES = frozenset(range(9))
ADJ = {u: {v for v in VERTICES if abs(u//3-v//3)+abs(u%3-v%3)==1}
       for u in VERTICES}


def paths(roots, branch, position=None):
    reach = {v: (v,) for v in sorted(roots)}
    waiting = deque(reach)
    while waiting:
        u = waiting.popleft()
        for v in sorted((ADJ[u] & branch)-reach.keys()):
            if position is None or position[u] > position[v]:
                reach[v] = reach[u]+(v,)
                waiting.append(v)
    return reach


def words(branch, prefix=()):
    if len(branch) == len(prefix):
        yield prefix
        return
    fringe = set().union(*(ADJ[v] for v in prefix)) & branch if prefix else branch
    for v in sorted(fringe-set(prefix)):
        yield from words(branch, prefix+(v,))


def change_order(word, path, monotone, branch, contacts, preimages):
    on_path = set(path)
    output = path+tuple(v for v in word if v not in on_path)
    assert set(output) == branch and output[0] in contacts
    exposed = {output[0]}
    for v in output[1:]:
        assert ADJ[v] & exposed
        exposed.add(v)
    position = {v: i for i,v in enumerate(word)}
    slots = ([0]+sorted(position[v] for v in path[:-1]) if monotone
             else [position[v] for v in path])
    assert (slots[0] == 0 if monotone else slots[-1] == 0)
    j = len(path)-1
    recovered = dict(zip(slots, reversed(path) if monotone else path, strict=True))
    remainder = iter(output[j+1:])
    assert tuple(recovered[i] if i in recovered else next(remainder) for i in range(len(word))) == word
    preimages[output,j] += 1


checks = Counter()
for bits in range(1, (1 << 9)-1):
    branch = frozenset(v for v in VERTICES if bits & (1 << v))
    if len(paths({min(branch)}, branch)) != len(branch):
        continue
    checks['connected_branches'] += 1
    candidates = list(words(branch))
    for w in sorted(VERTICES-branch):
        contacts = ADJ[w] & branch
        if not contacts:
            continue
        checks['contact_sets_with_multiplicity'] += 1
        shortest = paths(contacts, branch)
        normal, monotone = Counter(), Counter()
        for word in candidates:
            root = word[0]
            change_order(word, shortest[root], False, branch, contacts, normal)
            checks['geodesic_words'] += 1
            position = {v: i for i,v in enumerate(word)}
            ordered = paths(contacts, branch, position)
            assert root in ordered
            assert len(ordered[root]) >= len(shortest[root])
            change_order(word, ordered[root], True, branch, contacts, monotone)
            checks['monotone_words'] += 1
        for table, capacity in ((normal, perm),(monotone, comb)):
            for (_, length), multiplicity in table.items():
                assert multiplicity <= capacity(len(branch)-1, length)
                checks['inverse_buckets'] += 1

for key, count in checks.items():
    print(key, count)
