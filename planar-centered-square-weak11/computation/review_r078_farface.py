"""Independent far-face input: exact new cap rows, raw-prefix horizontal areas."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,os,platform,sys,time
start=time.monotonic();ip=Path(sys.argv[1]);pp=Path(sys.argv[2]);ir=ip.read_bytes();pr=pp.read_bytes();raw=json.loads(ir);claim=json.loads(pr)
assert hashlib.sha256(ir).hexdigest()=='e343e448b9582a731ddaccca37452adcce1c6f28b4153aa520a0ec2b4a3ef574'
points=[tuple(map(F,p))for p in raw['positions']];masses=list(map(F,raw['masses']));old=raw['S_indices']
assert old==list(range(400)) and points[:400]==[(F(21*i,400),F(21*j,400))for i in range(20)for j in range(20)]
assert masses[:400]==[F(1,400)]*400 and points[400:]==[(-F(1,2),F(0)),(F(0),-F(1,2))] and masses[400:]==[F(37,100)]*2
assert F(raw['threshold'])==1
h=F(21,400);e=F(1,400);B=F(37,50);prefix_caps=[];jointcaps=[]
for k in range(1,21):
    rhs=max(F(1),(F(1,2)+(k-1)*h)**2)-F(k,20)
    for tag in (1,2):assert rhs>F(37,100);prefix_caps.append((tag,k,20,rhs-F(37,100)))
for k in range(1,21):
    for l in range(1,21):
        rhs=max(F(1),max(F(1,2)+(k-1)*h,F(1,2)+(l-1)*h)**2)-F(k*l,400)
        assert rhs>B;prefix_caps.append((3,k,l,rhs-B));jointcaps.append(rhs)
assert len(prefix_caps)==440 and min(jointcaps)==F(1197,1600)
assert min(t[3]for t in prefix_caps)==F(13,1600)
assert all(F(37,100)*j<1 for j in (1,1,2))
nonterminal=(F(599,400))**2-F(35,33)*F(87,50);assert nonterminal==F(698811,1760000)>0
local=[i for i in old if all(-1<t<F(1,2)-e for t in points[i])];xx=sorted({points[i][0]for i in local});yy=sorted({points[i][1]for i in local})
assert len(local)==100 and len(xx)==len(yy)==10
L=400*(1<<72);unit=L//400;assert claim['integer_scale']==L and claim['bits']==72
def clip(r):
    z=[max(r[0],-L//2),min(r[1],-unit),max(r[2],-L//2),min(r[3],-unit)]
    return z if z[0]<z[1]and z[2]<z[3]else None
def make(k,l,tag):
    ids=[i for i in local if points[i][0]<=xx[k-1]and points[i][1]<=yy[l-1]]if k else[]
    assert len(ids)==k*l
    ids+=[400+j for j in range(2)if tag>>j&1]
    a,b=min(points[i][0]for i in ids),max(points[i][0]for i in ids)
    c,d=min(points[i][1]for i in ids),max(points[i][1]for i in ids)
    mass=sum((masses[i]for i in ids),F(0));span=max(b-a,d-c);delta=max((b-a)/2,(d-c)/2,b+e,d+e,-F(1,2)-a,-F(1,2)-c,F(0))
    row={'id':f'old{k}x{l}-outside{tag}','k':k,'l':l,'outside_mask':tag,'old_count':k*l,'mass':str(mass),'bbox_units_400':[int(v*400)for v in (a,b,c,d)],'eligible':mass>span*span,'Gamma_strict_active':mass>4*delta*delta,'twice_min_distance_units_400':int(800*delta)}
    assert mass<=1
    if row['eligible']:
        z=mass*L*L/4;low=math.isqrt(z.numerator//z.denominator);assert low*low<=z<(low+1)*(low+1)
        high=low if low*low==z else low+1
        assert z<=high*high and high<=L//2
        aa,bb,cc,dd=[int(v*L)for v in (a,b,c,d)]
        row.update(radius_lower_units=low,radius_upper_units=high,inner=clip((bb-low,aa+low,dd-low,cc+low)),outer=clip((bb-high,aa+high,dd-high,cc+high)))
    else:assert not row['Gamma_strict_active']
    return row
nodes=[make(k,l,t)for k in range(1,11)for l in range(1,11)for t in range(4)];res=[make(0,0,t)for t in (1,2,3)];nodes+=res
claimed={r['id']:r for r in claim['menu']};assert len(claimed)==len(nodes)==403
fieldchecks=0
for row in nodes:
    for key,val in row.items():assert val==claimed[row['id']].get(key),(row['id'],key);fieldchecks+=1
def rects(rows,key):return [r[key]for r in rows if r.get(key)is not None]
ri,ro=rects(res,'inner'),rects(res,'outer')
def merge(iv):
    out=[]
    for a,b in sorted(iv):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out
slabs=0
def horizontal(full,subtract):
    global slabs
    yy=sorted({v for r in full+subtract for v in r[2:]});answer=0
    for a,b in zip(yy,yy[1:]):
        f=merge([(x0,x1)for x0,x1,y0,y1 in full if 2*y0<a+b<2*y1]);r=merge([(x0,x1)for x0,x1,y0,y1 in subtract if 2*y0<a+b<2*y1])
        length=sum(y-x for x,y in f)
        for x0,x1 in f:
            for y0,y1 in r:length-=max(0,min(x1,y1)-max(x0,y0))
        assert length>=0;answer+=(b-a)*length;slabs+=1
    return answer
def enc(source,expected,si=None,so=None):
    low=F(horizontal(rects(source,'inner'),ro if so is None else so),L*L)
    high=F(horizontal(rects(source,'outer'),ri if si is None else si),L*L)
    assert low==F(expected['lower'])and high==F(expected['upper'])and low<=high
    return {'lower':str(low),'upper':str(high),'display':float((low+high)/2)}
total=enc(nodes,claim['actual_Gamma_deletion_area']);multi=[r for r in nodes if r['old_count']and r['outside_mask']==3 and r['k']>=2 and r['l']>=2];kept=[r for r in nodes if r not in multi]
extra=claim['farface_checks'];assert len(multi)==81 and len(kept)==322
multi_area=enc(multi,extra['multirow_multicolumn_after_complete_residual'])
kept_area=enc(kept,extra['area_after_deleting_all_multirow_multicolumn_joint_labels'])
loss=enc(multi,extra['net_area_lost_by_deleting_those_labels'],rects(kept,'inner')+ri,rects(kept,'outer')+ro)
assert F(loss['lower'])>F(1,250000) and F(total['upper'])<F(1,8)
individual=[];byjoint={(r['k'],r['l']):r for r in claim['joint_label_classification']}
for r in nodes:
    if r['old_count']and r['outside_mask']==3 and r['Gamma_strict_active']:
        z=enc([r],byjoint[r['k'],r['l']]['area_outside_complete_residual']);assert F(z['lower'])>0
        individual.append((r['k'],r['l']))
assert len(individual)==76 and sum(k>=2 and l>=2 for k,l in individual)==61
E=(-F(61,1000),-F(59,1000),-F(61,1000),-F(59,1000));exclusion=[];esaved={r['id']:r for r in extra['E_exclusion_rows']}
for r in kept:
    a,b,c,d=map(lambda v:F(v,400),r['bbox_units_400']);delta=max((b-a)/2,(d-c)/2,b-E[1],E[0]-a,d-E[3],E[2]-c,F(0))
    slack=4*delta*delta-F(r['mass']);assert slack>0
    assert str(delta)==esaved[r['id']]['minimum_farthest_distance']and str(slack)==esaved[r['id']]['four_distance_squared_minus_mass']
    exclusion.append(slack)
assert len(exclusion)==322 and min(exclusion)==F(2721,250000)
radius=F(89,200)
def far(z):return max(abs(z[0]-E[0]),abs(z[0]-E[1]),abs(z[1]-E[2]),abs(z[1]-E[3]))
def near(z):return max(E[0]-z[0],z[0]-E[1],E[2]-z[1],z[1]-E[3],F(0))
always=[i for i,z in enumerate(points)if far(z)<radius];possible=[i for i,z in enumerate(points)if near(z)<radius]
assert always==possible and len(always)==66 and sum(i<400 for i in always)==64
assert sum((masses[i]for i in always),F(0))==F(9,10)
assert max(far(points[i])for i in always)==F(441,1000)
assert F(9,10)-4*radius*radius==F(1079,10000)
active={str(t):sum(r['Gamma_strict_active']for r in nodes if r['old_count']and r['outside_mask']==t)for t in range(4)}
assert active==claim['active_counts_by_outside_mask']
out={'status':'PASS','method':'new440cap rows via frozen all-n theorem; raw402atom prefix reconstruction; independent horizontal yslab integer union/difference; rawwholecell witness and322labelminimax exclusions',
 'new_prefix_cap_checks':440,'pure_outside_cap_checks':3,'minimum_joint_capacity':'1197/1600','minimum_outside_prefix_slack':'13/1600','unique_maximum_by_strict_rows':True,'connected_by_origin_edges':True,'nonterminal_gap':str(nonterminal),
 'catalogue_labels':403,'catalogue_field_comparisons':fieldchecks,'active_old_label_counts':active,
 'actual_Gamma_deletion':total,'multirow_joint_after_residual':multi_area,'reduced_menu_Gamma_deletion':kept_area,'direct_missing_menu_area':loss,
 'positive_joint_label_count':len(individual),'positive_multirow_joint_label_count':61,'horizontal_slabs_executed':slabs,
 'E_certificate':{'cell':list(map(str,E)),'area':'1/250000','common_open_radius':'89/200','old_count':64,'captured_total_mass':'9/10','open_distance_margin':'1/250','mass_margin':'1079/10000','kept_labels_excluded':322,'minimum_exact_exclusion_slack':'2721/250000'},
 'scope':'One new specified full S20 input;440newparameterrows actually checked, no fullbboxcap re-enumeration. Exact regional areas show dropping all higher joint labels loses positive area, not a payment/C2 counterexample. Individual surviving labels are not summed.',
 'input_sha256':hashlib.sha256(ir).hexdigest(),'producer_result_sha256':hashlib.sha256(pr).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[3]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
