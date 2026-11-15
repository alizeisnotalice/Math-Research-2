import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from exact import tables,vertices
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
start=time.perf_counter();out={'python':sys.version,'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'results':[]}
for n,A2 in ((4,779299174080),(5,7063261481362958280000)):
    cuts=[f'{axis}{j}' for axis in ('v','h') for j in range(2,n-1)]
    counts={}
    for k in range(1,len(cuts)+1):
        for sub in itertools.combinations(cuts,k):
            count=int(subprocess.check_output([str(a.binary),str(n),*sub],text=True))
            counts[','.join(sub)]=count
    union=sum((1 if len(key.split(','))%2 else -1)*value for key,value in counts.items())
    S1=sum(value for key,value in counts.items() if ',' not in key)
    S2=S1+2*sum(value for key,value in counts.items() if len(key.split(','))==2)
    rect={}
    for j in range(2,n-1):
        xy,adj,one,two,h,roots=tables(n,j)
        full=(1<<(n*j))-1
        outer=sum(roots[full][u] for u,p in enumerate(xy) if p[0]<=(j-1)//2)
        rect[j]={'one':one[full],'outer_root':outer}
    for j in range(2,n-1):
        assert counts[f'v{j}']==math.comb(n*n,n*j)*rect[j]['outer_root']*rect[n-j]['outer_root']
    result={'n':n,'A2':A2,'rectangles':rect,'intersections':counts,'union_count':union,
            'union_ratio':str(Fraction(union,A2)),
            'second_moment_bound':str(Fraction(S1*S1,A2*S2))}
    out['results'].append(result);print(json.dumps(result),flush=True)
out['wall_seconds']=time.perf_counter()-start;a.output.write_text(json.dumps(out,indent=2)+'\n')
