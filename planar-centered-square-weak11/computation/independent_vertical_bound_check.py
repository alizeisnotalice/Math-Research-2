"""Independently verify saved vertical-family B&B leaves; never search.

Rebuild masses by x-column/y-row accumulation from the exact base input.
Use independent active-interval slab union for every saved leaf envelope.
Keep all potentially active boxes across the full parameter interval.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import argparse,bisect,hashlib,json,math,platform,sys,time
from independent_periodic_halo import setup_grid,rounded_radii,slab_union_area,ceil,floor


def rectangles(boxes,grid,px,py,a,b,inner=False):
    assert not inner or a==b
    Q=math.lcm(a.denominator,b.denominator);L,H=int(a*Q),int(b*Q)
    X,Y=int(px*grid['S']),int(py*grid['S'])
    out=set()
    for left,right,bottom,top,rlo,rhi in boxes:
        r=rlo if inner else rhi
        rr=(max(0,right-r)*Q,min(X,left+r)*Q,
            max(0,min(L*top,H*top)-r*Q),
            min(H*Y,max(L*bottom,H*bottom)+r*Q))
        if rr[0]<rr[1] and rr[2]<rr[3]:out.add(rr)
    return out,grid['S']*Q


def build(data,bmin,bmax,bits,cap):
    px,py,lam=(F(data[k]) for k in ('Px','Py','lam'))
    base=defaultdict(F)
    for x,y,m in data['pts']:
        x,y,m=F(x)%px,F(y)%py,F(m)
        assert m>=0
        if m:base[x,y]+=m
    mass=sum(base.values(),F(0))
    tail=mass/(px*bmin*py)*(1+px/(2*cap))*(1+bmin*py/(2*cap))
    assert tail<lam
    lifted=[];margin=cap/bmin
    for (x,y),m in base.items():
        for k in range(ceil((-cap-x)/px),floor((px+cap-x)/px)+1):
            for l in range(ceil((-margin-y)/py),floor((py+margin-y)/py)+1):
                lifted.append((x+k*px,y+l*py,m))
    grid=setup_grid(lifted,(px,py),bits,cap,lam)
    S=grid['S'];xs=sorted({p[0] for p in grid['points']});ys=sorted({p[1] for p in grid['points']})
    yi={y:i for i,y in enumerate(ys)};columns=defaultdict(list)
    for x,y,m in grid['points']:columns[x].append((yi[y],m))
    capi=grid['cap'];cache={};boxes=set();count=Counter()
    globalQ=math.lcm(bmin.denominator,bmax.denominator)
    L,H=int(bmin*globalQ),int(bmax*globalQ)
    X,Y=int(px*S),int(py*S)
    ystops=[bisect.bisect_right(ys,y+floor(2*capi/bmin)) for y in ys]
    for li,left in enumerate(xs):
        rows=[0]*len(ys)
        for right in xs[li:bisect.bisect_right(xs,left+2*capi)]:
            for j,m in columns[right]:rows[j]+=m
            for bottomi,bottom in enumerate(ys):
                m=0
                for topi in range(bottomi,ystops[bottomi]):
                    count['all_span_boxes']+=1;m+=rows[topi]
                    if not m:continue
                    count['positive_mass_boxes']+=1
                    radii=cache.get(m)
                    if radii is None:radii=rounded_radii(m,grid);cache[m]=radii
                    rlo,rhi=radii;top=ys[topi]
                    # A whole-family empty outer receiver cannot be needed.
                    if right-left>=2*rhi or bmin*(top-bottom)>=2*rhi:continue
                    if max(0,right-rhi)>=min(X,left+rhi):continue
                    if max(0,min(L*top,H*top)-rhi*globalQ)>=min(H*Y,max(L*bottom,H*bottom)+rhi*globalQ):continue
                    boxes.add((left,right,bottom,top,rlo,rhi))
    count.update({'lifted_atoms':len(lifted),'x_coordinates':len(xs),'y_coordinates':len(ys),
                  'unique_masses':len(cache),'retained_potentially_active_boxes':len(boxes)})
    return boxes,grid,px,py,lam,mass,tail,dict(count)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--certificate',type=Path,required=True)
    ap.add_argument('--round004',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    start=time.perf_counter();raw=args.input.read_bytes();cert_raw=args.certificate.read_bytes()
    data=json.loads(raw);cert=json.loads(cert_raw)
    assert hashlib.sha256(raw).hexdigest()==cert['input_sha256']
    bmin,bmax=F(cert['beta_min']),F(cert['beta_max']);bits=cert['bits']
    leaves=sorted((F(r['lo']),F(r['hi']),F(r['upper']),r['status']) for r in cert['leaves'])
    assert leaves and leaves[0][0]==bmin and leaves[-1][1]==bmax
    assert all(a<b and status in ('active','pruned') for a,b,u,status in leaves)
    assert all(leaves[i][1]==leaves[i+1][0] for i in range(len(leaves)-1))
    assert sum((b-a for a,b,u,status in leaves),F(0))==bmax-bmin
    assert len(leaves)==cert['divisions']+1
    boxes,grid,px,py,lam,mass,tail,counts=build(data,bmin,bmax,bits,F(5))
    assert tail==F(cert['tail_bound'])
    factor=lam/mass
    print(json.dumps({'phase':'independent family rebuilt','counts':counts,'elapsed_seconds':time.perf_counter()-start}),flush=True)
    incumbent=F(cert['incumbent_lower']);best=F(cert['best_beta'])
    rr,scale=rectangles(boxes,grid,px,py,best,best,True)
    own_incumbent=F(slab_union_area(rr)[0],scale*scale)*factor
    assert own_incumbent==incumbent
    results=[];equal=0
    for index,(a,b,saved,status) in enumerate(leaves):
        rr,scale=rectangles(boxes,grid,px,py,a,b)
        upper=F(slab_union_area(rr)[0],scale*scale)*factor
        assert upper<=saved,(index,a,b,upper,saved)
        equal+=upper==saved
        if status=='pruned':assert saved<incumbent
        results.append({'lo':str(a),'hi':str(b),'saved_upper':str(saved),
                        'independent_upper':str(upper),'equal':upper==saved,'status':status})
        if (index+1)%50==0:print(json.dumps({'leaves_verified':index+1,'total':len(leaves),'elapsed_seconds':time.perf_counter()-start}),flush=True)
    upper=max([incumbent]+[u for a,b,u,status in leaves])
    assert upper==F(cert['global_upper'])
    assert upper-incumbent==F(cert['certified_gap'])<=F(cert['target_gap'])
    assert cert['status']=='target_gap'
    samples={F(beta):(F(v['lower']),F(v['upper'])) for beta,v in cert['samples'].items()}
    assert all(bmin<=beta<=bmax and lo<=hi for beta,(lo,hi) in samples.items())
    assert max(lo for lo,hi in samples.values())==incumbent==samples[best][0]
    for a,b,u,status in leaves:
        for beta in (a,b):
            if beta in samples:assert samples[beta][1]<=u
    r4raw=args.round004.read_bytes();r4=json.loads(r4raw)
    separation=F(r4['weak_ratio_lower']['fraction'])-upper
    assert separation>0
    out={'input_sha256':hashlib.sha256(raw).hexdigest(),
         'leaf_certificate_sha256':hashlib.sha256(cert_raw).hexdigest(),
         'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':sys.version,'platform':platform.platform(),'verified_utc':datetime.now(timezone.utc).isoformat(),
         'algorithm':'independent x-column accumulation, complete uniform window, slab interval merging',
         'source_family_or_primary_code_imported':False,'branch_search_rerun':False,
         'counts':counts,'tail_bound':str(tail),'leaf_count':len(leaves),
         'leaf_coverage_exact':True,'leaf_bounds_equal_count':equal,
         'pruned_leaves':sum(s=='pruned' for a,b,u,s in leaves),
         'active_leaves':sum(s=='active' for a,b,u,s in leaves),
         'best_beta':str(best),'independently_verified_incumbent_lower':str(incumbent),
         'global_family_upper':str(upper),'certified_gap':str(upper-incumbent),
         'round004_certificate_sha256':hashlib.sha256(r4raw).hexdigest(),
         'round004_lower_minus_restricted_family_upper':str(separation),
         'round004_separation_decimal_approximation':float(separation),
         'scope':'only the fixed R3 vertical family on the recorded beta interval; not a general C2 upper bound',
         'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv],
         'leaves':results}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))

if __name__=='__main__':main()
