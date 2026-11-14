"""Independent small boundary checks; no rerun of the saved five large points."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,platform,sys,time
from datetime import datetime,timezone
from independent_periodic_halo import setup_grid,enumerate_boxes,slab_union_area,ceil,floor


def normalized(data):
    px,py=F(data['Px']),F(data['Py'])
    pts={}
    for x,y,m in data['pts']:
        key=F(x)%px,F(y)%py
        pts[key]=pts.get(key,F(0))+F(m)
    return [(x,y,m) for (x,y),m in pts.items() if m]


def independent_at(data,beta,bits):
    px,py,lam=F(data['Px']),F(data['Py'])*beta,F(data['lam'])
    cap=F(5);points=[]
    for x,y,m in normalized(data):
        y*=beta
        for k in range(ceil((-cap-x)/px),floor((px+cap-x)/px)+1):
            for l in range(ceil((-cap-y)/py),floor((py+cap-y)/py)+1):
                points.append((x+k*px,y+l*py,m))
    g=setup_grid(points,(px,py),bits,cap,lam)
    inner,outer,counts=enumerate_boxes(g)
    return tuple(F(slab_union_area(rects)[0],g['S']**2) for rects in (inner,outer))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    start=time.perf_counter();sys.path.insert(0,str(args.module.parent))
    spec=importlib.util.spec_from_file_location('vertical_family_under_audit',args.module)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    VF=module.VerticalFamily;bits=16;bmin,bmax=F(1,2),F(3,2)
    examples=[
      ('new_activation',{'lam':'1','Px':'20','Py':'20','pts':[['10','8','4'],['10','12','4']]},[F(1,2),F(2,3),F(1),F(3,2)]),
      ('periodic_boundary_negative_lift',{'lam':'1','Px':'10','Py':'10','pts':[['0','19/2','4']]},[F(1,2),F(3,4),F(1),F(5,4),F(3,2)]),
      ('duplicate_zero_and_periodic_normalization',{'lam':'1','Px':'10','Py':'10','pts':[['0','19/2','1'],['10','-1/2','3'],['1','1','0']]},[F(1,2),F(3,4),F(1),F(5,4),F(3,2)]),
      ('interior_value_exceeds_both_endpoints',{'lam':'1','Px':'20','Py':'20','pts':[['10','8','1'],['10','9','1']]},[F(1,2),F(1),F(6,5),F(3,2)]),
    ]
    records=[];area_records={}
    for name,data,betas in examples:
        family=VF(data,bmin,bmax,bits,F(5))
        whole=family.interval_upper(bmin,bmax)
        env,scale=family._rectangles(bmin,bmax)
        assert whole==F(slab_union_area(env)[0],scale*scale)
        results={}
        for beta in betas:
            point=family.evaluate(beta);ours=independent_at(data,beta,bits)
            assert point==ours,(name,beta,point,ours)
            assert point[1]<=whole
            if name in ('periodic_boundary_negative_lift','duplicate_zero_and_periodic_normalization'):
                assert point[0]<4<point[1]
            results[str(beta)]={'inner_area':str(point[0]),'outer_area':str(point[1]),'independent_complete_area_equal':True}
        if name=='new_activation':
            S=family.scale
            chosen=[box for box in family.boxes if box[:4]==(10*S,10*S,8*S,12*S)]
            assert chosen
            assert all(4*S>=2*box[5] and F(1,2)*4*S<2*box[4] for box in chosen)
        if name=='interior_value_exceeds_both_endpoints':
            assert F(results['6/5']['inner_area'])>max(F(results['1/2']['outer_area']),F(results['3/2']['outer_area']))
        records.append({'name':name,'data':data,'enumeration':family.enumeration,
                        'evaluations':results,'whole_interval_upper':str(whole)})
        area_records[name]=results
    assert area_records['periodic_boundary_negative_lift']==area_records['duplicate_zero_and_periodic_normalization']
    # Direct endpoint-envelope checks isolate negative slopes and moving clipping.
    models=[]
    for name,c,d,r,expected in [('negative_source_y',-100,-100,200,(F(0),F(3),F(0),F(3,2))),
                               ('moving_cell_top',900,900,200,(F(0),F(3),F(5,2),F(15)))]:
        family=object.__new__(VF);family.bmin=bmin;family.bmax=bmax
        family.px=family.py=F(10);family.scale=100;family.boxes=[(100,100,c,d,r,r)]
        rects,s=family._rectangles(bmin,bmax)
        assert len(rects)==1
        envelope=tuple(F(v,s) for v in next(iter(rects)))
        assert envelope==expected
        for i in range(21):
            beta=bmin+F(i,20)
            point,ss=family._rectangles(beta,beta)
            for rr in point:
                actual=tuple(F(v,ss) for v in rr)
                assert envelope[0]<=actual[0]<actual[1]<=envelope[1]
                assert envelope[2]<=actual[2]<actual[3]<=envelope[3]
        models.append({'name':name,'exact_interval_envelope':list(map(str,envelope)),'all_21_rational_subchecks':'contained'})
    out={'python':sys.version,'platform':platform.platform(),'run_at_utc':datetime.now(timezone.utc).isoformat(),
         'module_sha256':hashlib.sha256(args.module.read_bytes()).hexdigest(),
         'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'bits':bits,'beta_interval':[str(bmin),str(bmax)],'independent_small_cases':records,
         'direct_affine_endpoint_models':models,'saved_five_large_points_rerun':False,
         'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
