"""Independent bounded exact source/bbox oracle for position_slab_enumerator.
No large target geometry enumeration is performed.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
import argparse,hashlib,importlib.util,json,random,sys,time
PIN='acede2c710655d0c179e7bea3ee9483271caaf37968227e7e9fce6442907056b'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();assert sha(a.module)==PIN;spec=importlib.util.spec_from_file_location('reviewed_enumerator',a.module);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);rng=random.Random(390151);cases=[];total=0;catalog_checks=0
    inputs=[]
    for k in range(8):
        data={'lam':['2/3','5/7'][k%2],'pts':[[str(F(rng.randrange(-16,17),7)),str(F(rng.randrange(-15,16),11)),str(F(rng.randrange(1,8),13))]for _ in range(5)]};inputs.append((data,[rng.randrange(2)for _ in range(5)],F(k+1,1009),F(1,1000),F(3),[F(-5,3),F(7,4),F(-9,5),F(11,6)],False))
    for k in range(3):
        data={'lam':'2/3','Px':'13/2','Py':'17/2','pts':[['0','0','1/3'],['11/2','6','2/5'],['92/7','3','0']]};inputs.append((data,[1,k%2,1],F(-1,37)+F(k,101),F(1,10),F(3),[F(0),F(13,2),F(0),F(17,2)],True))
    for ix,(data,vel,s0,H,R,domain,periodic)in enumerate(inputs):
        out=mod.enumerate_slab(data,vel,s0,H,R,domain=domain,periodic=periodic);S=out['source_catalog'];xs=[(F(x),int(v))for x,v in out['x_classes']];ys=list(map(F,out['y_classes']));dt=F(out['delta_empty']);h=F(out['H']);lam=F(data['lam']);fl,fr,fb,ft=domain
        assert 0<dt<=h<=H
        expected=set()
        for j,((x,y,w),v)in enumerate(zip(data['pts'],vel)):
            x,y,w=map(F,(x,y,w))
            if not w:continue
            shifts=product(range(-20,21),repeat=2)if periodic else[(0,0)]
            for k,l in shifts:
                xx=x+v*s0+(k*F(data['Px'])if periodic else 0);yy=y+(l*F(data['Py'])if periodic else 0)
                if not periodic or (xx+v*H>=fl-R and xx-v*H<=fr+R and fb-R<=yy<=ft+R):expected.add((j,k,l,xx,yy,v,w))
        actual={(p['label'],*p['period_shift'],F(p['x']),F(p['y']),p['v'],F(p['w']))for p in S};assert expected==actual;catalog_checks+=1
        recs={tuple(z['bbox']):z for z in out['box_records']};assert len(recs)==len(out['boxes'])
        for t in(-dt,F(0),dt):
            current=[(F(p['x'])+p['v']*t,F(p['y']),F(p['w']))for p in S];active=set()
            for i,k in product(range(len(xs)),repeat=2):
                if i>k:continue
                lo,hi=xs[i][0]+xs[i][1]*t,xs[k][0]+xs[k][1]*t
                for j,l in product(range(len(ys)),repeat=2):
                    if j>l:continue
                    mass=sum(w for x,y,w in current if lo<=x<=hi and ys[j]<=y<=ys[l]);ids=i,k,j,l
                    if not mass:continue
                    r2=min(R*R,mass/(4*lam));positive=(hi-lo)**2<4*r2 and (ys[l]-ys[j])**2<4*r2 and max(F(0),fl-lo)**2<r2 and max(F(0),hi-fr)**2<r2 and max(F(0),fb-ys[j])**2<r2 and max(F(0),ys[l]-ft)**2<r2
                    if positive:
                        active.add(ids);assert ids in recs and F(recs[ids]['mass'])==mass and F(recs[ids]['radius_squared'])==r2
                    total+=1
            assert active<=set(recs)
            if not t:assert active==set(recs)
        margins=0
        for w in out['rejection_witnesses'].values():
            if 'radius'in w:
                q,L=map(F,(w['witness']['q'],w['witness']['L']));assert q>0 and L>0 and F(w['radius'])==min(h,q/(2*L));margins+=1
        cases.append({'case':ix,'periodic':periodic,'source_count':len(S),'active_boxes':len(recs),'H':str(h),'delta_empty':str(dt),'margin_witnesses_checked':margins})
    failures=[]
    for data,vel,kind,domain,R in [({'lam':'1','pts':[['0','0','1'],['0','0','1']]},[0,1],'source_order',[-10,11,-4,4],3),({'lam':'1','pts':[['0','0','4'],['3','0','5']]},[0,1],'x_mass_span',[-10,11,-4,4],3),({'lam':'1','pts':[['0','0','1'],['6','0','1']]},[0,1],'span_2R',[-10,11,-4,4],3),({'lam':'1','pts':[['14','0','1']]},[1],'reach_right',[0,11,-4,4],3),({'lam':'1','pts':[['-3','0','1']]},[1],'reach_left',[0,11,-4,4],3),({'lam':'1','pts':[['12','0','4']]},[1],'x_clip_right',[0,11,-4,4],3),({'lam':'1','pts':[['-1','0','4']]},[1],'x_clip_left',[0,11,-4,4],3)]:
        try:mod.enumerate_slab(data,vel,0,F(1,10),R,domain=domain,periodic=False)
        except mod.WallAtCenter as e:assert e.kind==kind;failures.append(kind)
        else:raise AssertionError('expected wall rejection '+kind)
    assert sha(a.module)==PIN
    out={'status':'passed','module_sha256':PIN,'cases':cases,'catalogues_independently_reconstructed':catalog_checks,'nonzero_bbox_radius_tests':total,'center_active_sets_exact':True,'no_rejected_box_activates_at_delta_endpoints':True,'walls_fail_closed':failures,'scope':'eleven finite/periodic bounded exact cases only; no real-target large enumeration replay','seed':390151,'script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
