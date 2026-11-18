"""Exact equality set and extended-affine orbit at Boolean spectral norm 3."""
from collections import Counter
import json
import platform

n = 4
N = 1 << n
U = (1 << N)-1
masks = [sum(1<<x for x in range(N) if (a & x).bit_count()%2) for a in range(N)]
canonical = sum(1<<x for x in range(N) if (x>>3&1)*((x>>2&1)^((x>>1&1)*(x&1))))
templates = set()
for f in range(1<<N):
    W = [N-2*(f^m).bit_count() for m in masks]
    if sum(map(abs,W)) != 3*N:
        continue
    assert Counter(map(abs,W)) == {0:6,4:8,8:2}
    balanced = sum(2*sum(abs(W[a]) for a in range(N) if (a&v).bit_count()%2)==3*N for v in range(1,N))
    assert balanced == 8
    templates.add(f)

def translate(f,i):
    mask = masks[1<<i]
    shift = 1<<i
    return ((f&mask)>>shift)|((f&(U^mask))<<shift)

def neighbors(f):
    yield f^U
    for i in range(n):
        yield f^masks[1<<i]
        translated=translate(f,i)
        yield translated
        for j in range(n):
            if i!=j:
                m=masks[1<<j]
                yield (translated&m)|(f&(U^m))

seen = {canonical}
todo = [canonical]
while todo:
    f = todo.pop()
    for image in neighbors(f):
        if image not in seen:
            seen.add(image)
            todo.append(image)
assert seen == templates
print(json.dumps({"python":platform.python_version(),"variables":n,"all_truth_tables":1<<N,
 "canonical_ANF":"w*(v+u1*u2)","canonical_truth_mask":canonical,"orbit_generators":
 ["translations","coordinate transvections","affine output summands"],
 "norm_three_functions":len(templates),"orbit_size":len(seen),"set_equality":True,
 "unnormalized_absolute_Walsh_histogram":{"0":6,"4":8,"8":2},"balanced_directions":8},indent=2))
