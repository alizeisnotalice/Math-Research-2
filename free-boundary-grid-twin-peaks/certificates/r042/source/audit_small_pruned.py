import itertools,json,subprocess,sys,os,hashlib
from pathlib import Path

binary=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False)
n=3;N=n*n
adj=[{u for u in range(N) if abs(u%n-v%n)+abs(u//n-v//n)==1} for v in range(N)]
vertical=[''.join(t) for t in itertools.product('12',repeat=n)]
curves=vertical+['h'+h for h in vertical]
regions=[]
for token in curves:
    horizontal=token.startswith('h');h=token[1:] if horizontal else token
    L={v for v in range(N) if (v//n if horizontal else v%n)<int(h[v%n if horizontal else v//n])};R=set(range(N))-L
    regions.append((L,R))
counts=[[0]*len(curves) for _ in curves];union=0;peak_checked=0
for word in itertools.permutations(range(N)):
    pos=[0]*N
    for i,v in enumerate(word):pos[v]=i
    hits=[]
    for i,(L,R) in enumerate(regions):
        peaks=[]
        for H in (L,R):
            peaks.append([v for v in H if all(pos[v]<pos[u] for u in adj[v]&H)])
        if any(len(p)!=1 for p in peaks):continue
        l,r=peaks[0][0],peaks[1][0]
        band=(n+3)//4
        horizontal=curves[i].startswith('h')
        xl=l//n if horizontal else l%n;xr=r//n if horizontal else r%n
        if not(xl<band and xr>=n-band):continue
        if adj[l]&R or adj[r]&L:continue
        hits.append(i)
        fullpeaks=[v for v in range(N) if all(pos[v]<pos[u] for u in adj[v])]
        assert set(fullpeaks)=={l,r}
        assert abs(l%n-r%n)+abs(l//n-r//n)>=2
        peak_checked+=1
    union+=bool(hits)
    for i in hits:
        for j in hits:counts[i][j]+=1
receipts=[]
for i,h in enumerate(curves):
    for j in range(i,len(curves)):
        p=out/f'{i}-{j}.json'
        child=subprocess.Popen([str(binary),str(n),str(p),h,curves[j]])
        rc=child.wait();assert rc==0
        assert int(json.loads(p.read_text())['count'])==counts[i][j],(h,curves[j])
        receipts.append(dict(pid=child.pid,returncode=rc,waited_reaped=True))
S1=sum(counts[i][i] for i in range(len(curves)));S2=sum(map(sum,counts))
result=dict(n=3,permutations=362880,pairs_checked=len(receipts),S1=S1,S2=S2,union_count=union,full_peak_checks=peak_checked,curves=curves,matrix=counts,children=receipts)
(out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('matrix','children')}))
