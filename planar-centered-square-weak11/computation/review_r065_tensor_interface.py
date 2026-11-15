#!/usr/bin/env python3
"""Bounded independent nonuniform-grid/D4/anisotropic/source-oracle audit."""
from collections import Counter,defaultdict
from fractions import Fraction as F
import ast,hashlib,importlib.util,json,os
from pathlib import Path
import sys,time

HERE=Path(__file__).resolve().parent
OWNER=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11')
CHECK=OWNER/'computation/r065_tensor_source_verify.py'
LP=OWNER/'computation/cauchy_occupation_tensor_jump_lp_v2.m'
GRID=OWNER/'computation/inputs/r065-tensor-occupation/central-fine-outer-coarse.json'
PINS={CHECK:'5124b452b656fb093a30f40b619bfbc31d7d3bb095385aa6b10f8e0fccc23e38',LP:'3d696e00cbbd4edb0c0a4392efd820c763bb54ac3d8c73f991d3d530418bbfe2'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def independent_grid(tokens):
    e=list(map(F,tokens));assert all(x<y for x,y in zip(e,e[1:]))and e==[-x for x in reversed(e)]
    c=[(x+y)/2 for x,y in zip(e,e[1:])];w=[y-x for x,y in zip(e,e[1:])]
    groups=defaultdict(set)
    for i,x in enumerate(c):
        for j,y in enumerate(c):groups[tuple(sorted((abs(x),abs(y)),reverse=True))].add((i,j))
    return e,c,w,[(key,groups[key])for key in sorted(groups)]


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    for p,h in PINS.items():assert sha(p)==h
    m=load('r065_checked',CHECK);ref=load('independent_ref',HERE/'review_r062_coupled_source.py')
    v3=OWNER/'computation/r062_coupled_source_verify_v3.py'
    def functions(p):return {n.name:ast.dump(n,include_attributes=False)for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
    old,new=functions(v3),functions(CHECK)
    for name in ('ceildiv','iv','plus','minus','times','scale','divide','square','power','atan_small','atan_pos','atan_signed','boxprob_point','cauchyjet','kprime','boxjet','split'):assert old[name]==new[name]
    tokens=json.loads(GRID.read_text())['axis_edges'];e,c,w,groups=independent_grid(tokens)
    ae,ac,aw,ao,ai=m.tensor_geometry(tokens)
    assert (ae,ac,aw)==(e,c,w)and len(groups)==len(ao)==276 and len(ai)==2025
    assert len(w)==45 and Counter(w)=={F(1,32):17,F(1,8):28}
    assert sum(w)==F(129,32)and sum(len(ids)*aw[next(i for i,x in enumerate(c)if x==key[0])]*aw[next(j for j,y in enumerate(c)if y==key[1])]for key,ids in groups)==F(16641,1024)
    for (key,ids),orb in zip(groups,ao):
        assert set(orb['indices'])==ids and orb['center']==key
        for i,j in ids:
            x,y=c[i],c[j];wx,wy=w[i],w[j]
            if (j,i)in ids:assert w[j]*w[i]==wx*wy
            assert wx*wy==orb['width'][0]*orb['width'][1]
    # Unequal-cell toy shows a common short-side error may falsely certify.
    corner=F(1,2);true_min=F(-1,2)
    correct=corner-F(2)*F(2)**2/8;wrong=corner-F(2)*F(1)**2/8
    assert correct==true_min<0<wrong
    cases=[('constant_nonuniform3',['-2','-1/2','1/2','2'],F(4),[F(1,10000)]*3),
           ('nonmonotone_nonuniform3',['-2','-1/2','1/2','2'],F(4),[F(1,10000),F(3,10000),F(2,10000)]),
           ('nonuniform5',['-3','-2','-1/4','1/4','2','3'],F(8),[F(i+1,100000)for i in range(6)])]
    rows=[];samples=0;boundary=0
    for name,tokens,h,coeff in cases:
        data=dict(axis_edges_exact=tokens,h=str(h),source_origin_mass='1',density_coefficients=list(map(str,coeff)))
        model=m.Model(data,F(0));ee,cc,ww,gg=independent_grid(tokens)
        grid={ij:a for a,(_,ids)in zip(coeff,gg)for ij in ids}
        mass=F(1)+sum(a*ww[i]*ww[j]for (i,j),a in grid.items());assert model.mass==mass
        # Independently sum per-cell left/right endpoint contributions.
        xedges=defaultdict(F);yedges=defaultdict(F)
        for (i,j),a in grid.items():
            xedges[ee[i],cc[j],ww[j]]+=a;xedges[ee[i+1],cc[j],ww[j]]-=a
            yedges[ee[j],cc[i],ww[i]]+=a;yedges[ee[j+1],cc[i],ww[i]]-=a
        assert {k:v for k,v in xedges.items()if v}=={(x,y,d):v for x,y,d,v in model.xjumps}
        assert {k:v for k,v in yedges.items()if v}=={(x,y,d):v for x,y,d,v in model.yjumps}
        for j in range(model.m):
            assert xedges[ee[0],cc[j],ww[j]]==grid[0,j]
            assert xedges[ee[-1],cc[j],ww[j]]==-grid[model.m-1,j];boundary+=2
        if name.startswith('constant'):assert len(model.xjumps)==len(model.yjumps)==2*model.m
        areas=F(0);checked=[]
        for target,B,_ in m.initial_boxes(model,1):
            o=model.oo[target];areas+=len(o['indices'])*(B[1]-B[0])*(B[3]-B[2])
            lo,err,req,neg,curv=model.box_bound(target,B);assert lo>0 and req==0 and not neg
            checked.append(dict(cell=target,box=list(map(str,B)),lower_fixed=str(lo),error_fixed=str(err)))
            for x in (B[0],(B[0]+B[1])/2,B[1]):
                for y in (B[2],(B[2]+B[3])/2,B[3]):
                    rx=ref.multiply(ref.reference_jet(x,h)[2],ref.reference_jet(y,h)[0]);ry=ref.multiply(ref.reference_jet(x,h)[0],ref.reference_jet(y,h)[2])
                    for (i,j),a in grid.items():
                        tx=ref.reference_prob(x-cc[i],ww[i],h);ty=ref.reference_prob(y-cc[j],ww[j],h)
                        ddx=ref.minus(ref.reference_jet(x-ee[i],h)[1],ref.reference_jet(x-ee[i+1],h)[1])
                        ddy=ref.minus(ref.reference_jet(y-ee[j],h)[1],ref.reference_jet(y-ee[j+1],h)[1])
                        rx=ref.plus(rx,ref.times(ref.multiply(ddx,ty),a));ry=ref.plus(ry,ref.times(ref.multiply(tx,ddy),a))
                    ref.contains(curv[0],rx,m.S);ref.contains(curv[1],ry,m.S);samples+=1
        assert areas==(ee[-1]-ee[0])**2
        rows.append(dict(name=name,axis=tokens,physical_cells=model.m**2,orbits=len(gg),mass_exact=str(mass),vertical_jumps=len(model.xjumps),horizontal_jumps=len(model.yjumps),checks=checked))
    for p,h in PINS.items():assert sha(p)==h
    out=dict(status='PASS_R65_NONUNIFORM_INTERFACE_BOUNDED_AUDIT',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,
             prepared_grid=dict(axis_intervals=45,physical_cells=2025,D4_orbits=276,fine_intervals=17,coarse_intervals=28,total_area='16641/1024'),
             scalar_interval_functions_AST_identical_to_reviewed_R62V3=True,anisotropic_error_counterexample=dict(correct_lower=str(correct),wrong_shared_short_side_lower=str(wrong)),
             toy_cases=rows,fullcell_derivative_reference_points=samples,outer_boundary_checks=boundary,
             source_pins=[dict(path=str(p),sha256=sha(p))for p in [CHECK,LP,GRID,v3,HERE/'review_r062_coupled_source.py']],
             scope='Independent exact prepared-grid/D4/mass/jump/anisotropic oracle checks plus3tiny nonuniform source cases. No MATLAB/LP run, no4416row exact certificate, no optimized source/fullD/generalC2 claim.')
    dest=HERE/'results/round065-tensor-interface-independent.json';dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('toy_cases','source_pins')},indent=2))


if __name__=='__main__':main()
