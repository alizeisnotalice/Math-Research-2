"""Independent full-grid peak verification of the improved boundary cut."""
from collections import Counter
from itertools import permutations
from pathlib import Path
import json
import platform
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'r048/merge-repair'))

from audit import branch_order, graph, merge, peaks, shell


def main():
    n=3
    N=n*n
    adj=graph(n)
    total=Counter()
    for order in permutations(range(N)):
        pos={v:j+1 for j,v in enumerate(order)}
        ps=peaks(adj,pos)
        if len(ps)!=2:
            continue
        p=order[0]
        q,=ps-{p}
        w,B=merge(adj,order,q)
        i=pos[w]
        outer=shell(adj,B)-{w}
        vulnerable={u for u in outer if all(pos[u]<pos[v] for v in adj[u] if v not in B)}
        normal=min((pos[v] for v in outer),default=N+1)
        cutoff=min((pos[v] for v in vulnerable),default=N+1)
        assert cutoff>=normal>i
        total['two_peak']+=1
        if cutoff>normal:
            total['strict_gain']+=1
        deleted=all(any(pos[u]<pos[v] for u in adj[v] if u not in B)
                    for v in range(N) if v not in B and v!=p)
        assert deleted==(not vulnerable)
        if deleted:
            total['empty_vulnerable']+=1
        qscale=len(B)+1
        out={v:j*qscale for v,j in pos.items()}
        for k,v in enumerate(branch_order(adj,B,w),1):
            out[v]=qscale*i+k
        assert max(out[v] for v in B)<qscale*cutoff
        assert peaks(adj,out)=={p}
    assert total['two_peak']==143808 and total['strict_gain']>0
    print(json.dumps(dict(python=platform.python_version(),audit=total),indent=2))

if __name__=='__main__':
    main()
