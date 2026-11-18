from deletion import rank,spaces
from collections import Counter
import json

def main():
    flats={h for v in spaces(4) for b in range(16) for h in [frozenset(x^b for x in v)]}
    cases=Counter()
    from functools import cache
    @cache
    def groups(s):
        out=[]
        for flat in flats:
            fibre=tuple(x+16 for x in s if x in flat)
            r=rank(fibre)
            isolated=[x-16 for x in fibre if rank(y for y in fibre if y!=x)<r]
            if isolated:
                out.append((16//len(flat),isolated))
        return out
    chars=[sum(1<<x for x in range(16) if (x&a).bit_count()%2) for a in range(16)]
    for table in range(65536):
        w=[16-2*(table^mask).bit_count() for mask in chars]
        for scale,isolated in groups(tuple(i for i,v in enumerate(w) if v)):
            assert all((scale*w[i])%16==0 for i in isolated)
            assert sum(abs(w[i]) for i in isolated)<=16
            cases[scale]+=1
    print(json.dumps(dict(tables=65536,affine_flats=len(flats),nonempty_coloop_tests=cases),indent=2))
if __name__=='__main__':main()
