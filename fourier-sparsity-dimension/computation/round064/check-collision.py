"""Integer check of the hash occupancy bound for every nonempty 3-cube subset."""
from itertools import product
from fractions import Fraction
import json,platform

def run(h):
    cases=0
    mins={}
    for k in range(4):
        M=1<<k
        maps=[]
        for col in product(range(M),repeat=h):
            images=[0]
            for a in col:
                images += [x^a for x in images]
            maps.append(images)
        margin=None
        for mask in range(1,1<<(1<<h)):
            src=[x for x in range(1<<h) if mask>>x&1]
            n=len(src)
            loss=Fraction(sum(n-len({f[x] for x in src}) for f in maps),len(maps))
            target=max(n-M,Fraction(n*(n-1),2*M+n-2))
            assert loss>=target
            if n>=2:
                d=loss-target
                margin=d if margin is None else min(d,margin)
            cases+=1
        mins[M]=str(margin)
    return {"python":platform.python_version(),"h":h,"all_nonempty_masks":(1<<(1<<h))-1,
            "M":[1,2,4,8],"cases":cases,"minimum_gap_n_at_least_two":mins}
print(json.dumps(run(3),indent=2))
