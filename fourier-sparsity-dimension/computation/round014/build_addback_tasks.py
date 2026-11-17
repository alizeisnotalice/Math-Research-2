#!/usr/bin/env python3
"""Generate all unital finite function codes and all genuinely unstable algebra/core addbacks."""
import argparse,hashlib,itertools,json,sys,time
from pathlib import Path
def span(b):
    s={0}
    for a in b:s|={v^a for v in list(s)}
    return s
def canonical(vs):
    piv={}
    for x in vs:
        for j in sorted(piv,reverse=True):
            if x>>j&1:x^=piv[j]
        if x:
            j=x.bit_length()-1
            for k in piv:
                if piv[k]>>j&1:piv[k]^=x
            piv[j]=x
    return tuple(piv[j] for j in sorted(piv,reverse=True))
def rref(n):
    for d in range(n+1):
        for piv in itertools.combinations(range(n),d):
            free=[(i,j) for i,p in enumerate(piv) for j in range(p+1,n) if j not in piv]
            for v in range(1<<len(free)):
                rows=[1<<p for p in piv]
                for k,(i,j) in enumerate(free):
                    if v>>k&1:rows[i]|=1<<j
                yield rows
def partitions(n):
    def rec(a):
        if len(a)==n:yield a;return
        for v in range(max(a,default=-1)+2):yield from rec(a+[v])
    yield from rec([])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--points',type=int,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--summary',type=Path,required=True);a=ap.parse_args()
    n=a.points;assert 2<=n<=7;beg=time.monotonic();full=(1<<n)-1
    parts=[tuple(sum(1<<j for j,v in enumerate(p) if v==c) for c in range(max(p)+1)) for p in partitions(n)]
    # Four copies of each old point; repetition into 64 points preserves the function algebra.
    embedC={g:sum(((g>>(i%(4*n)%n))&1)<<i for i in range(64)) for g in range(1<<n)}
    wm=sum(1<<i for i in range(64) if ((i%(4*n))//n)&1)
    zm=sum(1<<i for i in range(64) if ((i%(4*n))//(2*n))&1)
    counts={'points':n,'base_codes':0,'all_algebra_pairs':0,'stable_pairs_proved_without_search':0,'tasks':0,'max_lift_dimension':0}
    with a.out.open('w') as f:
        for base in rref(n-1):
            C=span([full]+base);Cb=canonical(C);counts['base_codes']+=1
            legal=[cells for cells in parts if all(cell in C for cell in cells)]
            kapC=max(sum(sum(g&~cell==0 for g in C).bit_length() for cell in cells) for cells in legal)
            for cells in legal:
                counts['all_algebra_pairs']+=1
                core=[g for g in C if all(g&e in C for e in cells)]
                if len(core)==len(C):
                    counts['stable_pairs_proved_without_search']+=1;continue
                Bb=canonical(core)
                basis=canonical([embedC[g] for g in Cb]+[embedC[g]&wm for g in Bb]+[embedC[e]&zm for e in cells])
                assert len(basis)<=16
                rec={'id':counts['tasks'],'points':n,'C_basis':Cb,'E_cells':cells,'B_basis':Bb,'base_capacity':kapC,'lift_basis':basis,'lift_dimension':len(basis),'upper':2*kapC}
                f.write(json.dumps(rec,separators=(',',':'))+'\n');counts['tasks']+=1
                counts['max_lift_dimension']=max(counts['max_lift_dimension'],len(basis))
    counts.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),tasks_sha256=hashlib.sha256(a.out.read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-beg,python=sys.version,scope='All unital codes on the specified finite set; stable pairs proved by selector containment; every other pair dispatched. Not a bent-model counterexample search.')
    a.summary.write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts),flush=True)
if __name__=='__main__':main()
