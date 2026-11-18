"""Certify directional contact errors and adaptive two-sided compatible-angle brackets."""
import argparse
import contextlib
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

from flint import arb, ctx
from certify_ah_margin_envelope import ANCHORED_CUTS, added_area, replace_line
from certify_ah_three_line_cut import dyadic, fit_archived_witnesses, format_area, specialized_source


def make_geometry(archive, extra):
    path = archive/'certify_ah_two_line_cap.py'
    code = specialized_source(path,F('.404'),F('.00046'))
    cuts = ANCHORED_CUTS[:-1] + ''.join(f', (F("{p}"), F("{d}"))' for p,d in extra) + ']'
    full = replace_line(code, ANCHORED_CUTS, cuts)
    ns = {'__file__':str(path),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(full.split('angular_cells=32')[0],str(path),'exec'), ns)
    return ns, full


def verify(ns, bins=None):
    dep = ('A','F_at','centered_points','sub','dot','cross','points','orders','sector','rows','n0','n1',
           'trig','rotate','placement','box','r3','root_seeds','low_root_seeds','deficit_by_bin','D')
    (A,F_at,vertices_for,sub,dot,cross,points,orders,sector,rows,n0,n1,trig,rotate,placement,box,r3,
     high,low,deficits,D) = [ns[k] for k in dep]
    # In either labeling: DC/EF contacts, reflected/ordinary cap points, DE contact.
    tol = list(map(A, ('.70','.70','.84','1','.84')))
    zero, one = A(0), A(1)
    side_errors=[]
    for j,(k,arc) in enumerate(sector):
        e, f = orders[j], orders[(j+1)%10]
        assert (e[1] if arc else e[0]) == k == (f[1] if arc else f[0])
        adj = (e[0], f[0]) if arc else (e[1], f[1])
        side_errors.append((k,arc,adj))
    def error(edges,u):
        cr=[cross(e,u) for e in edges]
        intervals=[]
        for j,(k,arc,(c1,c2)) in enumerate(side_errors):
            nxt=(j+1)%10
            if cr[j]<0 or cr[nxt]>0:
                continue
            det=cross(edges[j],edges[nxt]); assert det>0
            left=-cr[nxt]/det; right=cr[j]/det
            upper=min((max(tol[c1],tol[c2])*A('1.07')).upper(),
                      (left*tol[c1]+right*tol[c2]).upper())
            intervals.append((-upper,tol[k].upper()) if arc else (-tol[k].upper(),upper))
        assert intervals
        return min(r[0] for r in intervals).union(max(r[1] for r in intervals))
    result=[]
    orient_normals = {}
    for sgn in (-1,1):
        orient_normals[sgn] = [(name,(u[0],sgn*u[1]),u[0]-u[1]/r3,2*u[1]/r3) for name,u,c in rows]
    for j,E in enumerate(deficits):
        if E.is_zero() or (bins is not None and j not in bins):
            continue
        if j<9: turn,sign,cells=300,1,40
        elif j<30: turn,sign,cells=120,-1,40
        elif j<35: turn,sign,cells=300,1,75
        else: turn,sign,cells=0,1,90
        seeds=high if turn else low
        orientation=1 if turn else -1
        best_root=best_pose=one
        worst_root=worst_pose=None
        widths=[]
        for i in range(cells):
            a,b=(F(j,100)+F(i,100*cells), F(j,100)+F(i+1,100*cells))
            v=vertices_for(a,b); edges=[sub(v[k],v[l]) for k,l in orders]
            for k,e in enumerate(edges):
                f=edges[(k+1)%10]
                assert ((e[1]>0) if k<5 else (e[1]<0)) and cross(e,f)>0 and dot(e,f)>A('.75')
            for k in range(5):
                for l in range(k):
                    if (k,l) not in orders and (l,k) not in orders:
                        d=sub(v[k],v[l]); assert dot(d,d)<1
            _,gl,gn=points(box(a,b)); assert dot(gn,D)/dot(gl,D)<A('-3.5')
            t=F(2*i+1,2*cells); seed=((1-t)*seeds[j]+t*seeds[j+1])/1000
            endpoints=[]
            def sign_margin(phi,s):
                err=[error(edges,rotate(n0,trig(A(phi+60*k)))) for k in range(3)]
                return (s*(F_at(v,edges,phi)+E*(err[0]-err[1]+err[2]))).lower()
            for side in (-1,1):
                # Predict minimal useful width from neighboring cells, then verify each end.
                h=widths[-1][(side+1)//2] if widths else F(1)
                s=-side*orientation
                while sign_margin(seed+side*h,s) <= A('.000001'):
                    h+=F(1,20)
                    assert h < 5, (j,i,side,h)
                while h>=F(1,10) and sign_margin(seed+side*(h-F(1,20)),s)>A('.000001'):
                    h-=F(1,20)
                phi=seed+side*h; val=sign_margin(phi,s)
                if val<best_root: best_root,worst_root=val,[i,side,str(phi)]
                endpoints.append((phi,h))
            (al,h1),(bl,h2)=endpoints; widths.append((h1,h2))
            # Fixed 12 angular cells first. Recursively bisect only inconclusive comparisons.
            def check(left,right,depth):
                nonlocal best_pose,worst_pose
                ms=placement(v,edges,left,right,turn,sign)
                sn=trig(box(left-turn,right-turn)); n1p=(n1[0],sign*n1[1]); er0=error(edges,rotate(n0,sn)); er1=error(edges,rotate(n1p,sn))
                for (name,un,d0,d1),m in zip(orient_normals[sign],ms,strict=True):
                    cheap=(m-A('3.21')*E).lower()
                    # The rough bound is already rigorous. Refine near limiting edges.
                    val=cheap if cheap>A('.00012') else (m+E*(d0*er0+d1*er1-error(edges,rotate(un,sn)))).lower()
                    if val<=0:
                        assert depth < 4, (j,i,name,str(val),float(left),float(right))
                        mid=(left+right)/2
                        return check(left,mid,depth+1)+check(mid,right,depth+1)
                    if val<best_pose: best_pose,worst_pose=val,[i,name,str(left),str(right)]
                return 1
            if i==0: angle_cells=0
            for k in range(12):
                angle_cells+=check(al+(bl-al)*F(k,12),al+(bl-al)*F(k+1,12),0)
        row={'bin':j,'parameter_cells':cells,'angle_cells':angle_cells,'E_upper':str(E),'root_margin_lower':str(best_root),'root_index':worst_root,'placement_margin_lower':str(best_pose),'placement_index':worst_pose,'max_half_widths_degrees':[str(max(x[k] for x in widths)) for k in (0,1)]}
        result.append(row); print(j,str(best_root),str(best_pose), flush=True,file=sys.stderr)
    return result


def run(path, inp, out, bins):
    extra=[tuple(map(F,r)) for r in json.loads(inp.read_bytes())]
    ar,full=make_geometry(path,extra)
    assert all(0<p<F('.41') and 0<d<F('.00060') for p,d in extra)
    verified=verify(ar,bins)
    if bins is not None: return
    ctx.prec=160
    old={**ar,'cuts':ar['cuts'][:3],'tangents':ar['tangents'][:3]}
    shape,gain,inner=added_area(old,ar)
    prior_path=path/'ah-three-line-cut-certificate.json'; prior=json.loads(prior_path.read_bytes()); a,b=map(F,prior['area_exact_interval']); lo,hi=a-dyadic(gain.upper()),b-dyadic(gain.lower())
    assert lo<hi<a and hi < F('.844093395218375410388648671351126643458')
    obj={'lines':json.loads(inp.read_bytes()),'geometry_sha256':sha256(full.encode()).hexdigest(),'R167_area_sha256':sha256(prior_path.read_bytes()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ar['flint_version'],'geometry_and_placement_bits':128,'area_bits':160,'bins':verified,'shape':shape,'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),'R150_R156_witnesses_at_identity':fit_archived_witnesses(path,ar,inner)}
    out.write_text(json.dumps(obj,indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path); parser.add_argument('lines',type=Path); parser.add_argument('output',type=Path); parser.add_argument('--bins',nargs='+',type=int)
    opt=parser.parse_args(); sys.path.insert(0,str(opt.archive.resolve())); run(opt.archive.resolve(),opt.lines,opt.output,opt.bins)
