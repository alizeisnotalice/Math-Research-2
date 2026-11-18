from itertools import permutations
from math import factorial
from audit_early_fibres import comb_tree, local_neighbors


def code_word(word):
    if len(word) < 2:
        raise ValueError("need a word of length at least two")
    increases = [a < b for a, b in zip(word, word[1:])]
    peaks = sum(a and not b for a, b in zip(increases, increases[1:]))
    if peaks > 1:
        return None
    blocks = [[]]
    for i, value in enumerate(word):
        blocks[-1].append(value)
        if 0 < i < len(word) - 1 and increases[i-1] != increases[i]:
            blocks.append([])
    assert len(blocks) <= 4
    direction = increases[0]
    recon = tuple(v for i, block in enumerate(blocks)
                  for v in sorted(block, reverse=not (direction ^ bool(i % 2))))
    assert recon == word
    group = [0] * len(word)
    for i, block in enumerate(blocks):
        for value in block:
            group[value] = i
    return direction, tuple(group)


def audit_orders():
    for s in range(2, 10):
        found = set()
        for word in permutations(range(s)):
            code = code_word(word)
            if code is not None:
                assert code not in found
                found.add(code)
        assert len(found) <= min(factorial(s), 2*4**s)
        print(f"s={s}: admissible={len(found)}")


def audit_teeth():
    found=0
    for d in (3,4,8):
        for L in (8,16,32,64):
            A=comb_tree(L)|{(x,0) for x in range(1,d-1)}
            B={(d-x,y) for x,y in comb_tree(L)}
            teeth=[]
            positions=range(2*L+2,3*L+1,2)
            for row in positions:
                teeth.append(tuple((-x,row) for x in range(1,L+1)))
                teeth.append(tuple((-x,-row) for x in range(1,L+1)))
                teeth.append(tuple((-row,y) for y in range(1,L+1)))
            for side, selected in ((A,teeth),(B,[tuple((d-x,y) for x,y in path) for path in teeth])):
                assert len({v for path in selected for v in path})==len(selected)*L
                for path in selected:
                    assert all(v in side for v in path)
                    for before, v, after in zip(path,path[1:],path[2:]):
                        assert local_neighbors(side,v)=={before,after}
                    found+=1
    print(f"verified disjoint teeth across 12 comb states: {found}")


if __name__=="__main__":
    audit_orders()
    audit_teeth()
