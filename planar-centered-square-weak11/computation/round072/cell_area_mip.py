#!/usr/bin/env python3
"""R72 finite shared-measure whole-cell area MIP, with residual separation.

Source/cell geometry uses units 1/200. Every proposed selected cell is
certified against ALL outside bounding-box closures using Fraction masses.
MIP objectives/bounds have only finite support/cell/menu/margin scope.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,sys,time,warnings
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
from scipy.optimize import milp,Bounds,LinearConstraint,linprog
from scipy.io import savemat

D=40000;EPS=F(1,1000000)
S=[(48*i,48*j) for i in range(5) for j in range(5)]
O=[(48*i,48*j) for i in range(-4,9) for j in range(-4,9) if not(0<=i<=4 and 0<=j<=4)]+[(-60,47)]

def indices(mask):
    a=[]
    while mask:
      b=mask&-mask;a.append(b.bit_length()-1);mask-=b
    return a

def geometry(step=32):
    n=len(O);P=O+S;allout=(1<<n)-1;xs=sorted(set(x for x,y in P));ys=sorted(set(y for x,y in P))
    xm=[];ym=[]
    for i,a in enumerate(xs):
      for b in xs[i:]:xm.append((a,b,sum(1<<k for k,(x,y) in enumerate(P) if a<=x<=b)))
    for i,a in enumerate(ys):
      for b in ys[i:]:ym.append((a,b,sum(1<<k for k,(x,y) in enumerate(P) if a<=y<=b)))
    caps={}
    for a,b,mx in xm:
      for c,d,my in ym:
        full=mx&my;out=full&allout
        if not out:continue
        rhs=max(D,max(b-a,d-c)**2)-1600*(full>>n).bit_count()
        if out not in caps or rhs<caps[out][0]:caps[out]=(rhs,(a,b,c,d))
    masks=list(caps);rr=[indices(m) for m in masks];boxes=[]
    for ii in rr:
      pp=[O[k] for k in ii];boxes.append((min(p[0] for p in pp),max(p[0] for p in pp),min(p[1] for p in pp),max(p[1] for p in pp)))
    ri=[];ci=[]
    for i,row in enumerate(rr):ri.extend([i]*len(row));ci.extend(row)
    C=sparse.csr_matrix((np.ones(len(ri)),(ri,ci)),shape=(len(rr),n));capnum=np.array([caps[m][0] for m in masks],dtype=np.int64)
    cells=[(x,x+step,y,y+step) for x in range(-96,288,step) for y in range(-96,288,step)]
    if 384%step:raise ValueError("step must divide384")
    labels=[];bycell=[];pure=[]
    for cellid,(u,v,s,t) in enumerate(cells):
      dist=[max(x-u,v-x,y-s,t-y) for x,y in P]
      pp=[];always=False
      for r in sorted(set(d for d in dist if d<100)):
        internal=sum(d<=r for d in dist[n:]);hits=tuple(k for k,d in enumerate(dist[:n]) if d<=r)
        threshold=F(4*r*r-1600*internal,D)+EPS
        if threshold<=0:always=True;break
        if not internal or not hits:continue
        pp.append((r,hits,threshold))
      pure.append(always);bycell.append([])
      if not always:
        for r,hit,threshold in pp:bycell[-1].append(len(labels));labels.append((cellid,r,hit,threshold))
    boxes=np.array(boxes,dtype=np.int64)
    delta=[]
    for u,v,s,t in cells:
      a,b,c,d=boxes.T
      # 2*delta is integer, allowing half-unit bbox centers exactly.
      dd=np.maximum.reduce([b-a,d-c,2*(b-v),2*(u-a),2*(d-t),2*(s-c),np.zeros(len(a),dtype=np.int64)])
      delta.append(dd*dd) # 4*delta^2 in units /D
    return {"n":n,"masks":masks,"rows":rr,"boxes":boxes,"capnum":capnum,"C":C,"cells":cells,"labels":labels,"bycell":bycell,"pure":pure,"resnum":np.array(delta,dtype=np.int64),"cap_geometries":[caps[m][1] for m in masks]}

def matrix(g,cuts):
    n=g['n'];nc=len(g['cells']);nh=len(g['labels']);nv=n+nc+nh
    rows=[];vals=[];rhs=[];rowlabels=[]
    for i,(rr,b) in enumerate(zip(g['rows'],g['capnum'])):rows.append(rr);vals.append([1]*len(rr));rhs.append(float(F(int(b),D)));rowlabels.append(['maxcap',i])
    for h,(cell,r,hit,thr) in enumerate(g['labels']):rows.append(list(hit)+[n+nc+h]);vals.append([-1]*len(hit)+[float(thr)]);rhs.append(0.);rowlabels.append(['witness',h])
    for cell,hs in enumerate(g['bycell']):
      if g['pure'][cell]:continue
      rows.append([n+cell]+[n+nc+h for h in hs]);vals.append([1]+[-1]*len(hs));rhs.append(0.);rowlabels.append(['select',cell])
    for cell,k in sorted(cuts):
      cap=int(g['capnum'][k]);res=int(g['resnum'][cell,k]);M=F(cap-res,D)
      if M<=0:continue
      rows.append(g['rows'][k]+[n+cell]);vals.append([1]*len(g['rows'][k])+[float(M)]);rhs.append(float(F(cap,D)));rowlabels.append(['residual',cell,k])
    ri=[];ci=[];vv=[]
    for i,(rr,vs) in enumerate(zip(rows,vals)):ri.extend([i]*len(rr));ci.extend(rr);vv.extend(vs)
    A=sparse.csc_matrix((vv,(ri,ci)),shape=(len(rows),nv));b=np.array(rhs);c=np.zeros(nv);c[n:n+nc]=-1
    integrality=np.r_[np.zeros(n,dtype=int),np.ones(nc+nh,dtype=int)]
    return A,b,c,integrality,rowlabels

def certify(g,x):
    n=g['n'];nc=len(g['cells']);sel=[k for k in range(nc) if x[n+k]>.5];q=[F(max(0.,float(a))).limit_denominator(10**9) for a in x[:n]]
    cm=[sum(q[j] for j in rr) for rr in g['rows']];scale=F(1)
    for z,b in zip(cm,g['capnum']):
      if z>F(int(b),D):scale=min(scale,F(int(b),D)/z)
    for cell in sel:
      for z,b in zip(cm,g['resnum'][cell]):
        if z>F(int(b),D):scale=min(scale,F(int(b),D)/z)
    q=[a*scale for a in q];cert=[]
    for cell in sel:
      u,v,s,t=g['cells'][cell];dist=[max(p[0]-u,v-p[0],p[1]-s,t-p[1]) for p in O+S]
      candidates=[]
      for r in sorted(set(a for a in dist if a<100)):
        mu=sum(q[k] for k,d in enumerate(dist[:n]) if d<=r)+F(sum(d<=r for d in dist[n:]),25)
        margin=mu-F(4*r*r,D)
        if margin>0:candidates.append((margin,r,mu))
      if not candidates:return {"success":False,"reason":"rationalized full witness lost","failed_cell":cell,"scale":str(scale)}
      margin,r,mu=max(candidates)
      cert.append({'cell_id':cell,'bounds':[str(F(a,200)) for a in g['cells'][cell]],'closed_full_label_radius':str(F(r,200)),'captured_mass':str(mu),'strict_full_margin':str(margin),'label_outside_indices':[k for k,d in enumerate(dist[:n]) if d<=r],'label_S_indices':[k for k,d in enumerate(dist[n:]) if d<=r]})
    # Exact final all-row audit after global scaling.
    cm=[sum(q[j] for j in rr) for rr in g['rows']]
    assert all(z<=F(int(b),D) for z,b in zip(cm,g['capnum']))
    assert all(z<=F(int(b),D) for cell in sel for z,b in zip(cm,g['resnum'][cell]))
    area=sum(F((g['cells'][k][1]-g['cells'][k][0])*(g['cells'][k][3]-g['cells'][k][2]),D) for k in sel)
    return {'success':True,'selected_cells':cert,'selected_count':len(sel),'certified_area':str(area),'certified_area_float':float(area),'outside_weights':[str(a) for a in q],'rationalization_scale':str(scale),'all_caps_exact_checked':len(g['rows']),'all_selected_cell_residual_rows_exact_checked':len(g['rows'])*len(sel),'smallest_strict_full_margin':str(min((F(a['strict_full_margin']) for a in cert),default=F(0)))}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=int,default=32);ap.add_argument('--out',required=True);ap.add_argument('--threads',type=int,default=2);args=ap.parse_args();t=time.monotonic();g=geometry(args.step)
    n=g['n'];nc=len(g['cells']);cuts=set();masklookup={m:i for i,m in enumerate(g['masks'])}
    # Exact zero-distance cell rows are cheap and prevent placing mass inside
    # a selected residual-free cell before any MIP separation round.
    for k,(u,v,s,tt) in enumerate(g['cells']):
      m=sum(1<<j for j,(x,y) in enumerate(O) if u<=x<=v and s<=y<=tt)
      if m:cuts.add((k,masklookup[m]))
    print(json.dumps({'event':'built','caps':len(g['rows']),'cells':nc,'pureS_full_covered_cells':sum(g['pure']),'full_labels':len(g['labels']),'binaries':nc+len(g['labels']),'initial_cuts':len(cuts),'seconds':time.monotonic()-t}),flush=True)
    ledger=[];certificate=None;iteration=0
    while True:
      A,b,c,integrality,rowlabels=matrix(g,cuts);it=time.monotonic()
      with warnings.catch_warnings(record=True) as ww:
        res=milp(c,integrality=integrality,bounds=Bounds(np.zeros(len(c)),np.ones(len(c))),constraints=LinearConstraint(A,-np.inf,b),options={'disp':False,'presolve':True,'mip_rel_gap':0.,'threads':args.threads,'random_seed':0,'mip_feasibility_tolerance':1e-9})
      if res.x is None:raise RuntimeError(res.message)
      masses=g['C']@res.x[:n];selected=np.flatnonzero(res.x[n:n+nc]>.5).tolist();new=set();worst=0.
      for k in selected:
        violation=masses-g['resnum'][k]/D;worst=max(worst,float(violation.max()))
        bad=np.flatnonzero(violation>1e-8)
        if len(bad):
          order=bad[np.argsort(violation[bad])[-5:]]
          new.update((k,int(j)) for j in order if (k,int(j)) not in cuts)
      item={'iteration':iteration,'status':res.status,'message':res.message,'selected':len(selected),'objective_cellcount':float(-res.fun),'mip_gap':float(res.mip_gap),'mip_nodes':int(res.mip_node_count),'rows':A.shape[0],'columns':A.shape[1],'nnz':A.nnz,'cuts':len(cuts),'new_cuts':len(new),'worst_full_residual_violation':worst,'seconds':time.monotonic()-it,'warnings':[str(a.message) for a in ww]};ledger.append(item);print(json.dumps(item),flush=True)
      if not new:
        certificate=certify(g,res.x)
        if not certificate['success']:raise RuntimeError(certificate)
        break
      cuts.update(new);iteration+=1
    out={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parameters':{'source_grid_n':5,'S_side':'24/25','candidate_outside_count':145,'cell_step_units_1over200':args.step,'strict_mip_full_margin':str(EPS),'threads':args.threads,'seed':0,'mip_feasibility_tolerance':1e-9,'mip_relative_gap':0.,'radius_menu':'all full-source max-distance critical labels R_C<1/2, with at least one S point'},'scope':'Exact positive shared primal cell-area lower bound on one fixed S/O/cell/menu family. Floating MIP optimum/upper not general bound. Closed receiver cells, closed full source labels with strict positive margin imply slightly enlarged open full witnesses.','ledger':ledger,'certificate':certificate,'seconds':time.monotonic()-t}
    p=Path(args.out);p.write_text(json.dumps(out,indent=2)+'\n');stem=str(p.with_suffix(''))
    q=list(map(F,certificate['outside_weights']));positive=[(O[k],a) for k,a in enumerate(q) if a>0]
    inp={'threshold':'1','background':'Lebesgue R2','geometry':'centered axis-parallel common-radius squares','S_indices':list(range(25)),'positions':[[str(F(x,200)),str(F(y,200))] for x,y in S]+[[str(F(x,200)),str(F(y,200))] for (x,y),a in positive],'masses':['1/25']*25+[str(a) for pp,a in positive]}
    Path(stem+'-actual-input.json').write_text(json.dumps(inp,indent=2)+'\n')
    savemat(stem+'.mat',{'A':A,'b':b[:,None],'f':c[:,None],'lb':np.zeros((len(c),1)),'ub':np.ones((len(c),1)),'intcon':(np.flatnonzero(integrality)+1)[:,None],'scipy_x':res.x[:,None],'scipy_cellcount':-res.fun},do_compression=True)
    Path(stem+'-model-exact.json').write_text(json.dumps({'units_denominator':200,'S_positions_units':S,'outside_positions_units':O,'cap_rows':g['rows'],'cap_rhs_numerators':g['capnum'].tolist(),'rhs_denominator':D,'residual_bbox_units':g['boxes'].tolist(),'cells_units':g['cells'],'selected_residual_cut_indices':sorted(cuts),'strict_full_margin':str(EPS),'labels':[{'cell':a,'radius_units':r,'outside_indices':hit,'threshold':str(thr)} for a,r,hit,thr in g['labels']],'row_labels':rowlabels},separators=(',',':'))+'\n')
    print(json.dumps({'event':'finished','output':str(p),'selected_count':certificate['selected_count'],'certified_area':certificate['certified_area'],'positive_outside':len(positive),'seconds':out['seconds']}),flush=True)

if __name__=='__main__':main()
