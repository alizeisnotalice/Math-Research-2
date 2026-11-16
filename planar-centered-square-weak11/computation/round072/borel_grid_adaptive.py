#!/usr/bin/env python3
"""Exact adaptive actual halo area for the R72 Borel shared-source family.

sigma: n*n atoms, mass 1/n^2, positions (ih,jh), h=(n+1)/n^2.
nu: density n/(n+1) on [-1,2]^2 minus J=[-1/n^2,1]^2.
Relies on the separately proved square cap and empty residual halo.
No finite atomic approximation or subtraction transfer is made.
"""
from bisect import bisect_right
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import heapq
import json
import os
from pathlib import Path
import platform
import sys
import time


def area(c):return(c[1]-c[0])*(c[3]-c[2])


def overlap(c,r,lo,hi):return max(0,min(c+r,hi)-max(c-r,lo))


def classify(cell,n,D,grid,J):
    """Return +1 inner / -1 excluded / 0 unknown, with a replayable witness."""
    bounds=((cell[0],cell[1]),(cell[2],cell[3]))
    mid=(J[0]+J[1])//2
    assert 2*mid==J[0]+J[1]
    near=[max(a,min(b,mid))for a,b in bounds]
    far=[a if abs(a-mid)>=abs(b-mid)else b for a,b in bounds]
    gdist=[sorted(max(abs(p-a),abs(p-b))for p in grid)for a,b in bounds]
    pdist=[sorted(max(a-p,p-b,0)for p in grid)for a,b in bounds]
    # Multiplication by n^2(n+1)D^2 makes every test integral:
    # count(n+1)D^2 > n^3 overlap_units + 4n^2 radius_units^2.
    coeff=(n+1)*D*D
    events=sorted({r for ds in gdist for r in ds if 2*r<D})
    for r in events:
        count=bisect_right(gdist[0],r)*bisect_right(gdist[1],r)
        if not count:continue
        ov=overlap(near[0],r,*J)*overlap(near[1],r,*J)
        margin=count*coeff-n**3*ov-4*n*n*r*r
        if margin>0:
            return 1,{"radius_units":r,"guaranteed_atom_count":count,
                      "maximum_J_overlap_units_squared":ov,"strict_integer_margin":margin}
    events=sorted({r for ds in pdist for r in ds if 2*r<D})
    margins=[];checked=0
    for r in events:
        count=bisect_right(pdist[0],r)*bisect_right(pdist[1],r)
        if not count:continue
        checked+=1
        ov=overlap(far[0],r,*J)*overlap(far[1],r,*J)
        margin=n**3*ov+4*n*n*r*r-count*coeff
        if margin<0:return 0,None
        margins.append(margin)
    return -1,{"checked_possible_distance_events":checked,
               "minimum_nonnegative_integer_margin":min(margins)if margins else None}


def run(n,target,max_depth,max_nodes):
    start=time.perf_counter()
    assert n>=2 and target>0 and max_depth>=1
    # All endpoints, J midpoint and dyadic split points are integer units.
    D=2*n*n*(1<<max_depth)
    grid=[i*(n+1)*D//(n*n)for i in range(n)]
    assert all(p*n*n==i*(n+1)*D for i,p in enumerate(grid))
    J=(-D//(n*n),D)
    domain=(-D//2,(n*n-1)*D//(n*n)+D//2)*2
    queue=[(-area(domain),0,0,domain)];serial=0;processed=0
    unresolved=area(domain);inside=outside=0;leaves=[];counts={-1:0,0:0,1:0}
    while queue and processed<max_nodes and F(unresolved,D*D)>target:
        _,_,depth,cell=heapq.heappop(queue)
        state,witness=classify(cell,n,D,grid,J);processed+=1
        if state:
            aa=area(cell);unresolved-=aa
            if state==1:inside+=aa
            else:outside+=aa
            counts[state]+=1
            leaves.append({"state":state,"depth":depth,"cell_units":cell,"witness":witness})
        elif depth>=max_depth:
            counts[0]+=1
            leaves.append({"state":0,"depth":depth,"cell_units":cell,"witness":None})
        else:
            a,b,c,d=cell;xx=(a+b)//2;yy=(c+d)//2
            assert a<xx<b and c<yy<d
            for child in ((a,xx,c,yy),(xx,b,c,yy),(a,xx,yy,d),(xx,b,yy,d)):
                serial+=1;heapq.heappush(queue,(-area(child),serial,depth+1,child))
    for _,_,depth,cell in queue:
        counts[0]+=1
        leaves.append({"state":0,"depth":depth,"cell_units":cell,"witness":None})
    assert inside+outside+unresolved==area(domain)==sum(area(l["cell_units"])for l in leaves)
    assert unresolved==sum(area(l["cell_units"])for l in leaves if l["state"]==0)
    summary={"schema":"r072-borel-grid-adaptive-v1","n":n,"density_t":str(F(n,n+1)),
             "grid_spacing_h":str(F(n+1,n*n)),"grid_width_w":str(1-F(1,n*n)),
             "J":[str(-F(1,n*n)),"1",str(-F(1,n*n)),"1"],
             "receiver_domain":[str(F(v,D))for v in domain],"integer_coordinate_scale":D,
             "area_lower":str(F(inside,D*D)),"area_upper":str(F(inside+unresolved,D*D)),
             "area_lower_decimal":float(F(inside,D*D)),"area_upper_decimal":float(F(inside+unresolved,D*D)),
             "unresolved_area":str(F(unresolved,D*D)),"target_error":str(target),
             "target_reached":F(unresolved,D*D)<=target,"leaf_counts":{str(k):v for k,v in counts.items()},
             "processed_cells":processed,"max_depth":max_depth,"max_nodes":max_nodes,
             "strict_upper_below_3point5":2*(inside+unresolved)<7*D*D,
             "scope":"actual Borel density outside sources; area Delta=H(sigma+nu) using separately proved empty residual halo; no finite atomic approximation or new C2 bound",
             "assumptions":"all square mass <= max(1,side^2); residual halo empty; strict full witness radius<1/2; halo localized to saved receiver domain",
             "arithmetic":"integer predicates and fractions.Fraction output; no float decisions; displayed decimals only",
             "wall_seconds":time.perf_counter()-start}
    cert={"schema":"r072-borel-grid-cell-certificate-v1","n":n,"integer_coordinate_scale":D,
          "grid_units":grid,"J_axis_units":J,"receiver_domain_units":domain,
          "legend":{"1":"closed cell wholly in full halo, one common strict radius","-1":"closed cell outside full halo at every radius","0":"unresolved, charged only to upper area"},
          "leaves":leaves}
    return summary,cert


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,required=True);p.add_argument('--target-error',default='1/50')
    p.add_argument('--max-depth',type=int,default=16);p.add_argument('--max-nodes',type=int,default=200000)
    p.add_argument('--out-dir',type=Path,required=True);a=p.parse_args()
    assert os.path.ismount('/Volumes/data')
    summary,cert=run(a.n,F(a.target_error),a.max_depth,a.max_nodes)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    cp=a.out_dir/('n'+str(a.n)+'.cells.json.gz')
    with cp.open('wb')as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0)as zipped:
            zipped.write((json.dumps(cert,separators=(',',':'))+'\n').encode())
    summary['certificate_file']=cp.name
    summary['certificate_sha256']=hashlib.sha256(cp.read_bytes()).hexdigest()
    summary['receipt']={'python':sys.version,'platform':platform.platform(),'pid':os.getpid(),
                        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        'parameters':{'n':a.n,'target_error':a.target_error,'max_depth':a.max_depth,'max_nodes':a.max_nodes}}
    (a.out_dir/('n'+str(a.n)+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k]for k in ('n','area_lower','area_upper','area_lower_decimal','area_upper_decimal','target_reached','processed_cells','leaf_counts','wall_seconds')}))


if __name__=='__main__':main()
