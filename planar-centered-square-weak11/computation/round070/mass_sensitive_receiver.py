#!/usr/bin/env python3
"""S-only mixed-receiver envelope, analytic mass elimination, exact cell bounds.

For U captured from S, a=m(U)<L, r=max_U distance, c=max_S distance,
minimize sqrt(L+m-a)-sqrt(m)/2 on max(a,4r*r)<=m<=L.
The minimizer is q=max(a,4r*r,(L-a)/3).  Strict boundaries may be added
to an outer enclosure; lower boxes require strict inequalities throughout.
No grid in m, floating geometry, optimization, or random samples is used.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import heapq
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2)+'\n')


def dist2_bounds(cell,bbox):
    """Return twice the minimum/maximum farthest L-infinity distance."""
    x0,x1,y0,y1=cell;a,b,c,d=bbox
    lo=max(b-a,2*(b-x1),2*(x0-a),d-c,2*(d-y1),2*(y0-c))
    hi=2*max(b-x0,x1-a,d-y0,y1-c)
    assert 0<=lo<=hi
    return lo,hi


def above_minimum(C2,R2,a,L):
    """Exact c>min_m f_a(m), where C2=2c, R2=2r, units fixed.

    q3=3q; z4=4(L-a)+q3-C2².  For z4>=0, the inequality
    c*sqrt(q)>z4/4 is equivalent to 4*C2²*q3>3*z4².
    """
    q3=max(3*a,3*R2*R2,L-a)
    if q3>3*L:return False
    z4=4*(L-a)+q3-C2*C2
    return z4<0 or 4*C2*C2*q3>3*z4*z4


def classify(cell,candidates,labels,Sbox,L):
    cmin,cmax=dist2_bounds(cell,Sbox);possible=[]
    for j in candidates:
        lab=labels[j];rmin,rmax=dist2_bounds(cell,lab['bbox'])
        if rmax*rmax<L and above_minimum(cmin,rmax,lab['mass'],L):
            return 1,()  # closed cell lies strictly inside one envelope label
        if rmin*rmin<=L and above_minimum(cmax,rmin,lab['mass'],L):
            possible.append(j)
    return (0,tuple(possible)) if possible else (-1,())


def area(cell):return (cell[1]-cell[0])*(cell[3]-cell[2])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--job',required=True)
    ap.add_argument('--out-dir',type=Path,required=True);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data') and __debug__
    start=time.perf_counter();sp=args.spec if args.spec.is_absolute() else args.archive_root/args.spec
    spec=json.loads(sp.read_text());job=next(j for j in spec['jobs'] if j['id']==args.job)
    ip=args.archive_root/job['input'];assert sha(ip)==job['sha256'];data=json.loads(ip.read_text())
    gp=args.archive_root/spec['geometry'];assert sha(gp)==spec['geometry_sha256']
    loader=importlib.util.spec_from_file_location('frozen_r064_geometry',gp)
    mod=importlib.util.module_from_spec(loader);loader.loader.exec_module(mod)
    pts=[tuple(map(F,p)) for p in data['points']];mass=list(map(F,data['initial_masses']))
    S=data['S_mask'];n=len(pts);ids=[i for i in range(n) if S>>i&1];rem=[i for i in range(n) if not S>>i&1]
    assert min(mass)>0 and len(pts)==len(set(pts))==len(mass) and ids
    S_only=data.get('mode')=='S_only';assert bool(rem)!=S_only
    base=math.lcm(*(v.denominator for p in pts for v in p),*(m.denominator for m in mass))
    bits=spec['bits'];D=base*(1<<bits);den=D*D
    ps=[tuple(v*base for v in p) for p in pts];ms=[m*base*base for m in mass]
    ipts=[ps[i] for i in ids];ims=[ms[i] for i in ids]
    rpts=[ps[i] for i in rem];rms=[ms[i] for i in rem]
    full=mod.closures(ps);inside=mod.closures(ipts);res=mod.closures(rpts)
    L=sum(mass[i] for i in ids);Lu=int(L*den);assert Lu==L*den
    eligible=[b for b in full if sum(ms[i] for i in b['ids'])>b['span2']]
    maximum=max(sum(ms[i] for i in b['ids']) for b in eligible)
    maxima=[b['mask'] for b in eligible if sum(ms[i] for i in b['ids'])==maximum]
    assert maximum==L*base*base and S in maxima
    labels=[]
    for b in inside:  # all bbox closures, not merely internally eligible ones
        aa=sum(ims[i] for i in b['ids']);a=int(aa*(1<<(2*bits)))
        if a>=Lu:continue
        bb=tuple(int(v*(1<<bits)) for v in b['bbox'])
        labels.append({'mask':b['mask'],'mass':a,'bbox':bb,
                       'internal_eligible':aa>b['span2']})
    assert len(labels)==len(inside)-1
    Sbox=(min(int(pts[i][0]*D) for i in ids),max(int(pts[i][0]*D) for i in ids),
          min(int(pts[i][1]*D) for i in ids),max(int(pts[i][1]*D) for i in ids))
    rad=math.isqrt(Lu//4)
    if 4*rad*rad<Lu:rad+=1
    domain=(Sbox[0]-rad,Sbox[1]+rad,Sbox[2]-rad,Sbox[3]+rad)
    pending=[(-area(domain),0,0,domain,tuple(range(len(labels))))]
    leaves=[];serial=0;processed=0;budget=spec['max_nodes'];depthmax=spec['max_depth']
    while pending:
        _,_,depth,cell,candidates=heapq.heappop(pending)
        if processed>=budget:
            leaves.append((0,depth,cell));continue
        state,poss=classify(cell,candidates,labels,Sbox,Lu);processed+=1
        if state or depth>=depthmax:
            leaves.append((state,depth,cell));continue
        x0,x1,y0,y1=cell;xm=(x0+x1)//2;ym=(y0+y1)//2
        assert x0<xm<x1 and y0<ym<y1
        for child in [(x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1)]:
            serial+=1;heapq.heappush(pending,(-area(child),serial,depth+1,child,poss))
    assert sum(area(c) for _,_,c in leaves)==area(domain)
    mi=[c+(0,) for state,_,c in leaves if state==1]
    mo=[c+(0,) for state,_,c in leaves if state!=-1]
    ii,io=mod.exact_rectangles(ipts,ims,inside,bits)
    ri,ro=mod.exact_rectangles(rpts,rms,res,bits)
    DA=mod.difference_area
    def enc(lo,hi):
        assert 0<=lo<=hi
        return {'lower':str(F(lo,den)),'upper':str(F(hi,den)),
                'normalized_lower':str(F(lo,den)/L),'normalized_upper':str(F(hi,den)/L)}
    ilo,ihi=DA(ii,[]),DA(io,[])
    mlo,mhi=DA(mi,[]),DA(mo,[])
    m_ilo,m_ihi=DA(mi,io),DA(mo,ii)
    m_irlo,m_irhi=DA(mi,io+ro),DA(mo,ii+ri)
    creditlo=DA(ii,[])-DA(ii,ri);credithi=DA(io,[])-DA(io,ro)
    regions={'M_S_envelope':enc(mlo,mhi),'I':enc(ilo,ihi),
      'M_minus_I_S_only':enc(m_ilo,m_ihi),
      'I_union_M_S_only_budget':enc(ilo+m_ilo,ihi+m_ihi),
      'M_minus_I_union_R_diagnostic':enc(m_irlo,m_irhi),
      'internal_old_credit':enc(creditlo,credithi)}
    baseline=None
    if 'baseline' in data:
        bp=args.archive_root/data['baseline']['path'];assert sha(bp)==data['baseline']['sha256']
        baseline=json.loads(bp.read_text());assert baseline['S_mask']==S and F(baseline['L'])==L
        assert baseline['input_sha256']==data['provenance']['sha256']
        for key,old in [('true_P','P_pure_mixed'),('true_Delta','Delta')]:
            regions[key]=baseline['regions'][old]['global']
        assert m_irhi>=F(regions['true_P']['lower'])*den
    elif not S_only:
        fi,fo=mod.exact_rectangles(ps,ms,full,bits)
        plo,phi=DA(fi,io+ro),DA(fo,ii+ri);dlo,dhi=DA(fi,ro),DA(fo,ri)
        regions.update({'true_P':enc(plo,phi),'true_Delta':enc(dlo,dhi)})
        assert m_irhi>=plo  # interval compatibility, not a standalone inclusion proof
    cert={'schema':'r070-exact-envelope-cell-certificate-v1','id':args.job,'scale':D,
          'mass_unit_denominator':den,'L_units':Lu,'S_bbox':Sbox,'domain':domain,
          'labels':labels,'leaves':[{'state':s,'depth':dep,'cell':c} for s,dep,c in leaves],
          'legend':{'1':'strictly inside envelope','-1':'envelope excluded','0':'undecided; included only in outer area'}}
    args.out_dir.mkdir(parents=True,exist_ok=True);cp=args.out_dir/(args.job+'.cells.json.gz')
    with cp.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as g:
            g.write((json.dumps(cert,separators=(',',':'))+'\n').encode())
    unknown=sum(area(c) for s,_,c in leaves if s==0)
    output={'schema':'r070-analytic-mass-envelope-v1','status':'PASS','id':args.job,
      'python':sys.version,'source_sha256':sha(Path(__file__)),'spec_sha256':sha(sp),
      'input_sha256':sha(ip),'geometry_sha256':sha(gp),'source_provenance':data['provenance'],
      'source_count':n,'S_mask':S,'S_rank':len(ids),'L':str(L),'all_maximum_masks':maxima,
      'ambient_A_supplied':not S_only,'historical_true_geometry_reused':data.get('baseline'),
      'full_closed_bboxes':len(full),'internal_all_closed_bboxes':len(inside),
      'proper_envelope_labels':len(labels),'internally_ineligible_envelope_labels':sum(not b['internal_eligible'] for b in labels),
      'regions':regions,'unresolved_cell_area':enc(unknown,unknown),
      'cell_counts':{str(s):sum(ss==s for ss,_,_ in leaves) for s in [-1,0,1]},
      'processed_nodes':processed,'max_nodes':budget,'max_depth':depthmax,
      'certificate_file':cp.name,'certificate_sha256':sha(cp),
      'S_only_budget_at_most_3point5_certified':F(ihi+m_ihi,den)<=F(7,2)*L,
      'S_only_budget_strictly_above_3point5_certified':F(ilo+m_ilo,den)>F(7,2)*L,
      'wall_seconds':time.perf_counter()-start,
      'scope':'Three selected inputs only; S-only envelope uses all proper internal bbox closures and exact analytic elimination of m. Undecided cells retained outward. Uniform4 has no ambient A; historical true P/Delta reused from pinned R69 output, and no previous W is recalculated. No general theorem follows from these finite cases.'}
    write(args.out_dir/(args.job+'.json'),output)
    print(json.dumps({'status':'PASS','id':args.job,'S_only_budget':regions['I_union_M_S_only_budget'],
                      'true_P':regions.get('true_P'),'seconds':output['wall_seconds']}))


if __name__=='__main__':main()
