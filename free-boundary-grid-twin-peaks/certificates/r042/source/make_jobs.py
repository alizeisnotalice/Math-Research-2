import itertools,json,sys
from pathlib import Path
n=int(sys.argv[1]);out=Path(sys.argv[2]);assert n in (3,5)
heights=(n//2,n//2+1)
curves=[''.join(map(str,t)) for t in itertools.product(heights,repeat=n)]
def transforms(h):
    c=''.join(str(n-int(x)) for x in h)
    return (h,h[::-1],c,c[::-1])
classes={}
for i,h in enumerate(curves):
    for j in range(i,len(curves)):
        k=curves[j]
        canon=min(tuple(sorted(p)) for p in zip(transforms(h),transforms(k)))
        z=classes.setdefault(canon,dict(curves=list(canon),diagonal_multiplicity=0,offdiagonal_multiplicity=0,members=[]))
        z['diagonal_multiplicity' if i==j else 'offdiagonal_multiplicity']+=1
        z['members'].append([h,k])
jobs=[dict(id=i,**v) for i,(k,v) in enumerate(sorted(classes.items()))]
assert sum(z['diagonal_multiplicity'] for z in jobs)==len(curves)
assert sum(z['offdiagonal_multiplicity'] for z in jobs)==len(curves)*(len(curves)-1)//2
out.write_text(json.dumps(dict(n=n,curve_count=len(curves),curves=curves,jobs=jobs,symmetry='simultaneous row reversal and left-right reflection; unordered pairs'),indent=2)+'\n')
print('n',n,'curves',len(curves),'orbit_jobs',len(jobs))
