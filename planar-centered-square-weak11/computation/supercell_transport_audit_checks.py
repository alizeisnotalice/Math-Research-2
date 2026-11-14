"""Independent exact transport-lift and tiny geometry checks; no search."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,copy,hashlib,importlib.util,json,math,platform,re,sys,time
from independent_periodic_halo import setup_grid,enumerate_boxes,slab_union_area,ceil,floor

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copies(base,nx,ny):
    px,py=F(base['Px']),F(base['Py']);d=copy.deepcopy(base)
    d['Px'],d['Py']=str(nx*px),str(ny*py)
    d['pts']=[[str(F(x)+a*px),str(F(y)+b*py),w] for a in range(nx) for b in range(ny) for x,y,w in base['pts']]
    return d

def independent_area(data,bits):
    px,py,lam=map(F,(data['Px'],data['Py'],data['lam']));R=F(5)
    raw=[(F(x)%px,F(y)%py,F(m)) for x,y,m in data['pts'] if F(m)>0]
    mass=sum(p[2] for p in raw)
    assert mass/(px*py)*(1+px/10)*(1+py/10)<lam
    lifted=[]
    for x,y,m in raw:
        for k in range(ceil((-R-x)/px),floor((px+R-x)/px)+1):
            for l in range(ceil((-R-y)/py),floor((py+R-y)/py)+1):lifted.append((x+k*px,y+l*py,m))
    grid=setup_grid(lifted,(px,py),bits,R,lam);inn,out,count=enumerate_boxes(grid)
    return tuple(F(slab_union_area(rr)[0],grid['S']**2) for rr in (inn,out))

def square_mass(data,cx,cy,r):
    px,py=F(data['Px']),F(data['Py']);total=F(0)
    for x,y,m in data['pts']:
        x,y,m=F(x),F(y),F(m)
        nx=max(0,floor((cx+r-x)/px)-ceil((cx-r-x)/px)+1)
        ny=max(0,floor((cy+r-y)/py)-ceil((cy-r-y)/py)+1)
        total+=m*nx*ny
    return total

def lifts(data,base,factors,labels):
    px,py=F(data['Px']),F(data['Py']);bx,by=F(base['Px']),F(base['Py']);answer=[]
    for (x,y,m),(a,b,i) in zip(data['pts'],labels):
        sx,sy,w=map(F,base['pts'][i]);x,y=F(x)%px,F(y)%py
        dx,dy=x-(sx%bx+a*bx),y-(sy%by+b*by)
        hx=min((dx-floor(dx/px)*px,dx-ceil(dx/px)*px),key=abs)
        hy=min((dy-floor(dy/py)*py,dy-ceil(dy/py)*py),key=abs)
        assert F(m)==w
        answer.append({'label':list((a,b,i)),'hx':str(hx),'hy':str(hy),'linf':str(max(abs(hx),abs(hy)))})
    return answer

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.perf_counter()
    spec=importlib.util.spec_from_file_location('transport_under_audit',args.module)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    calc=module.calculate_transport;bits=20;results={}
    base={'Px':'3','Py':'5','lam':'1','pts':[['299/100','499/100','1/100']]}
    data=copies(base,2,2);last=data['pts'][-1]
    last[0]=str((F(last[0])+F(1,50))%6);last[1]=str((F(last[1])+F(1,50))%10)
    result=calc(data,base,(2,2),F(1,50),bits=bits)
    own=independent_area(data,bits)
    assert own==(F(result['inner_area']),F(result['outer_area']))
    labels=[(a,b,0) for a in range(2) for b in range(2)]
    signed=lifts(data,base,(2,2),labels)
    assert max(F(z['linf']) for z in signed)==F(1,50)
    checks=0
    for cx in (F(0),F(3),F(6),F(1,100)):
        for cy in (F(0),F(5),F(10),F(1,100)):
            for r in (F(1,100),F(1,50),F(1,10),F(1),F(5)):
                assert square_mass(data,cx,cy,r)<=square_mass(base,cx,cy,r+F(1,50));checks+=1
    results['two_axis_supercell_wrap']={'data':data,'signed_label_lifts':signed,
       'independent_inner_area':str(own[0]),'independent_outer_area':str(own[1]),
       'closed_square_transport_inequality_checks':checks}
    def reject(name,data,base,factors,delta,labels=None):
        try:calc(data,base,factors,delta,bits=bits,labels=labels)
        except ValueError as e:results[name]={'rejected':True,'message':str(e)};return
        raise AssertionError(name+' accepted')
    reject('delta_just_below_true_wrap',data,base,(2,2),F(19,1000))
    b={'Px':'2','Py':'3','lam':'3/2','pts':[['1/10','1/5','1/100']]}
    d=copies(b,2,1);d['pts'].reverse()
    reject('equal_mass_parity_swap_without_labels',d,b,(2,1),F(0))
    result=calc(d,b,(2,1),F(0),bits=bits,labels=[(1,0,0),(0,0,0)])
    own=independent_area(d,bits)
    assert own==(F(result['inner_area']),F(result['outer_area']))
    results['explicit_parity_labels_and_nonunit_lambda']={'inner_area':str(own[0]),'outer_area':str(own[1]),'lambda':'3/2'}
    b={'Px':'3','Py':'5','lam':'1','pts':[['1','1','1/100'],['2','2','1/50']]}
    d=copies(b,2,1);d['pts'][0][2],d['pts'][1][2]=d['pts'][1][2],d['pts'][0][2]
    reject('mass_swap_preserves_total_but_breaks_label_mass',d,b,(2,1),F(0))
    # Bind the supplied diff to the actual two source versions, ignoring dates.
    original=args.module.with_name('periodic_halo.py')
    provided=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/transport-tail/original_to_transport.diff')
    oldlines=original.read_text().splitlines(); patch=provided.read_text().splitlines()
    rebuilt=[];position=0;index=2
    while index<len(patch):
        header=re.fullmatch(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@.*',patch[index])
        assert header,patch[index]
        oldstart,oldcount,newstart,newcount=header.groups()
        oldstart=int(oldstart)-1;oldcount=int(oldcount or 1);newcount=int(newcount or 1)
        assert oldstart>=position
        rebuilt.extend(oldlines[position:oldstart]);position=oldstart
        seen_old=seen_new=0;index+=1
        while index<len(patch) and not patch[index].startswith('@@'):
            line=patch[index];assert line and line[0] in ' +-'
            sign,content=line[0],line[1:]
            if sign in ' -':
                assert oldlines[position]==content
                position+=1;seen_old+=1
            if sign in ' +':rebuilt.append(content);seen_new+=1
            index+=1
        assert seen_old==oldcount and seen_new==newcount
    rebuilt.extend(oldlines[position:])
    assert rebuilt==args.module.read_text().splitlines()
    p4=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/r004_final_0504.json')
    p8=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r008_final_0611.json')
    d4,d8=json.loads(p4.read_text()),json.loads(p8.read_text())
    inv=lambda z:(F(z['Px']),F(z['Py']),F(z['lam']),sum(F(p[2]) for p in z['pts']))
    assert inv(d4)==inv(d8)
    px,py,lam,m=inv(d8);rho=m/(px*py);R=F(5)
    bound=lambda delta:rho*(1+(px+2*delta)/(2*R))*(1+(py+2*delta)/(2*R))
    safe,rejected=bound(F(1,32)),bound(F(1,31));assert safe<lam<rejected
    A,B,C=2*R+px,2*R+py,4*R*R/rho
    discr=(A-B)**2+4*C;D=1<<80;rt=math.isqrt(discr.numerator*D*D//discr.denominator)
    assert rt*rt*discr.denominator<=discr.numerator*D*D<(rt+1)**2*discr.denominator
    dslo,dshi=(-(A+B)+F(rt,D))/4,(-(A+B)+F(rt+1,D))/4
    assert F(1,32)<dslo<dshi<F(1,31)
    assert F(1,2)<1 and F(1,2)*(1+F(1)/(2*F(1,10)))**2==18
    probe_root=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/transport-tail')
    probe=json.loads((probe_root/'probe_transport_r5.json').read_text())
    probe_input=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round010/supercell-probe/candidate_0015.json')
    pd=json.loads(probe_input.read_text());matching=probe['transport_matching']
    assert sha(probe_input)==probe['input_sha256']
    assert matching['base_input']==d8 and matching['copy_factors']==[1,2]
    canonical=lambda z:hashlib.sha256(json.dumps(z,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert canonical(pd)==matching['data_canonical_sha256'] and canonical(d8)==matching['base_canonical_sha256']
    assert set(map(tuple,matching['labels']))=={(0,b,i) for b in (0,1) for i in range(51)}
    signed_probe=lifts(pd,d8,(1,2),matching['labels'])
    probe_delta=F(matching['delta']);assert max(F(z['linf']) for z in signed_probe)==F(matching['maximum_torus_distance'])==probe_delta<=F(1,32)
    probe_bound=bound(probe_delta)
    assert probe_bound==F(probe['tail_bound'])==F(matching['transport_tail_bound'])<lam
    assert F(pd['Px'])==px and F(pd['Py'])==2*py
    assert sum(F(z[2]) for z in pd['pts'])==2*m
    pp=json.loads((probe_root/'probe_primary_r8.json').read_text())
    pi=json.loads((probe_root/'probe_independent_r8.json').read_text())
    assert pp['input_sha256']==pi['input_sha256']==probe['input_sha256']
    assert F(pp['tail_radius'])==F(pi['radius_cap'])==8
    assert F(probe['inner_area'])==F(pp['inner_area'])==F(pi['area_lower']['fraction'])
    assert F(probe['outer_area'])==F(pp['outer_area'])==F(pi['area_upper']['fraction'])
    probe_record={'input_path':str(probe_input),'input_sha256':sha(probe_input),
      'copy_factors':[1,2],'labels_independently_lifted':len(signed_probe),
      'maximum_linf_lift':str(probe_delta),'tail_bound_recomputed':str(probe_bound),
      'R5_new_vs_R8_original_and_independent_areas_equal':True,
      'large_geometry_rerun':False,'signed_label_lifts':signed_probe}
    checker_path=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-supercell-proof/planar-centered-square-weak11/computation/check_supercell_transport_tail.py')
    checker_spec=importlib.util.spec_from_file_location('tail_diagnostic_under_audit',checker_path)
    checker=importlib.util.module_from_spec(checker_spec);checker_spec.loader.exec_module(checker)
    regression=checker.check_nonunit_lambda_regression()
    assert F(regression['lambda'])==3 and F(regression['delta'])==F(1,5)
    assert F(regression['C'])==12 and F(regression['tail'])==F(289,100)<3
    out={'python':sys.version,'platform':platform.platform(),'verified_utc':datetime.now(timezone.utc).isoformat(),
      'module_sha256':sha(args.module),'diff_sha256':sha(provided),'diff_matches_actual_sources':True,
      'audit_script_sha256':sha(Path(__file__)),'small_exact_cases':results,
      'R4_input_sha256':sha(p4),'R8_input_sha256':sha(p8),'R4_R8_period_mass_level_equal':True,
      'density':str(rho),'G_1_over_32_at_5':str(safe),'safe_margin':str(lam-safe),
      'G_1_over_31_at_5':str(rejected),'reject_excess':str(rejected-lam),
      'delta_star_rational_interval':[str(dslo),str(dshi)],
      'supercell_q':ceil(R/min(2*px,py))+1,
      'fixed_R_positive_root_caveat':'rho<1 alone does not imply G_0(R)<1; Px=Py=1,m=1/2,R=1/10 gives G_0=18',
      'saved_102_row_probe_independent_matching':probe_record,
      'diagnostic_checker_sha256':sha(checker_path),'nonunit_lambda_root_regression':regression,
      'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('small_exact_cases','saved_102_row_probe_independent_matching')},indent=2))

if __name__=='__main__':main()
