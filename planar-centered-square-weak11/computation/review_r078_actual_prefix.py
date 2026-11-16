"""Raw-atom exact Gamma catalogue/area audit; no producer geometry imports."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import hashlib,json,math,os,platform,sys,time
start=time.monotonic();ip=Path(sys.argv[1]);pp=Path(sys.argv[2]);qp=Path(sys.argv[3]);ir=ip.read_bytes();pr=pp.read_bytes();qr=qp.read_bytes()
assert hashlib.sha256(ir).hexdigest()=='d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a'
assert hashlib.sha256(qr).hexdigest()=='8c45f7adcf9fe0b3ba1ae6eb99a417f6923714e4dae7a43808eae4c010e7121a'
raw=json.loads(ir);claimed=json.loads(pr);qualification=json.loads(qr)
assert qualification['status']=='PASS' and qualification['actual402_input_sha256']==hashlib.sha256(ir).hexdigest()
points=[tuple(map(F,p))for p in raw['positions']];mass=list(map(F,raw['masses']));old=raw['S_indices']
assert old==list(range(400)) and len(points)==len(mass)==402 and F(raw['threshold'])==1
assert points[:400]==[(F(21*i,400),F(21*j,400))for i in range(20)for j in range(20)] and mass[:400]==[F(1,400)]*400
outside=[i for i in range(402)if i not in set(old)];assert outside==[400,401]
e=F(1,400);edge=F(1,2)-e
local=[i for i in old if all(-1<t<edge for t in points[i])]
xx=sorted({points[i][0]for i in local});yy=sorted({points[i][1]for i in local});assert len(local)==100 and len(xx)==len(yy)==10
L=400*(1<<72);unit=L//400;Gamma=(-L//2,-unit,-L//2,-unit)
assert claimed['bits']==72 and claimed['integer_scale']==L and claimed['input_sha256']==hashlib.sha256(ir).hexdigest()
def clip(rr):
    a,b,c,d=rr;z=(max(a,Gamma[0]),min(b,Gamma[1]),max(c,Gamma[2]),min(d,Gamma[3]))
    return list(z)if z[0]<z[1]and z[2]<z[3]else None
def construct(k,l,tag):
    ids=[i for i in local if points[i][0]<=xx[k-1]and points[i][1]<=yy[l-1]]if k else[]
    assert len(ids)==k*l
    ids +=[outside[j]for j in range(2)if tag>>j&1]
    a=min(points[i][0]for i in ids);b=max(points[i][0]for i in ids);c=min(points[i][1]for i in ids);d=max(points[i][1]for i in ids)
    total=sum((mass[i]for i in ids),F(0));span=max(b-a,d-c)
    delta=max((b-a)/2,(d-c)/2,b+e,d+e,-F(1,2)-a,-F(1,2)-c,F(0))
    active=total>4*delta*delta;eligible=total>span*span;assert total<=1
    bbox=[int(t*400)for t in (a,b,c,d)];assert all(F(q,400)==t for q,t in zip(bbox,(a,b,c,d)))
    node={'k':k,'l':l,'outside_mask':tag,'old_count':k*l,'mass':str(total),'bbox_units_400':bbox,'eligible':eligible,'Gamma_strict_active':active,'twice_min_distance_units_400':int(delta*800)}
    assert F(node['twice_min_distance_units_400'],800)==delta
    if eligible:
        radius_squared=total*L*L/4;low=math.isqrt(radius_squared.numerator//radius_squared.denominator)
        assert low*low<=radius_squared<(low+1)*(low+1)
        high=low if low*low==radius_squared else low+1
        assert radius_squared<=high*high and high<=L//2
        aa,bb,cc,dd=(int(t*L)for t in (a,b,c,d))
        node.update(radius_lower_units=low,radius_upper_units=high,inner=clip((bb-low,aa+low,dd-low,cc+low)),outer=clip((bb-high,aa+high,dd-high,cc+high)))
    else:assert not active
    return node
own=[construct(k,l,tag)for k in range(1,11)for l in range(1,11)for tag in range(4)]
res=[construct(0,0,tag)for tag in (1,2,3)];own+=res
assert len(own)==403
reported={(r['k'],r['l'],r['outside_mask']):r for r in claimed['menu']};assert len(reported)==403
field_checks=0
for node in own:
    other=reported[node['k'],node['l'],node['outside_mask']]
    for key,value in node.items():assert value==other.get(key),(key,node,other);field_checks+=1
def rects(nodes,key):return [n[key]for n in nodes if n.get(key)is not None]
fi,fo,ri,ro=(rects(nodes,key)for nodes,key in ((own,'inner'),(own,'outer'),(res,'inner'),(res,'outer')))
for key,value in [('full_inner_rectangles',fi),('full_outer_rectangles',fo),('residual_inner_rectangles',ri),('residual_outer_rectangles',ro)]:assert value==claimed[key]
def merge(intervals):
    out=[]
    for a,b in sorted(intervals):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out
slab_count=0
def horizontal(full,residual):
    global slab_count
    ycuts=sorted({q for rr in full+residual for q in rr[2:]});answer=0
    for a,b in zip(ycuts,ycuts[1:]):
        ff=merge([(x0,x1)for x0,x1,y0,y1 in full if 2*y0<a+b<2*y1])
        rr=merge([(x0,x1)for x0,x1,y0,y1 in residual if 2*y0<a+b<2*y1])
        length=sum(x1-x0 for x0,x1 in ff)
        for x0,x1 in ff:
            for z0,z1 in rr:length-=max(0,min(x1,z1)-max(x0,z0))
        assert length>=0;answer+=(b-a)*length;slab_count+=1
    return answer
def enclosure(source,expected):
    low=F(horizontal(rects(source,'inner'),ro),L*L);high=F(horizontal(rects(source,'outer'),ri),L*L)
    assert low<=high and low==F(expected['lower'])and high==F(expected['upper'])
    return {'lower':str(low),'upper':str(high),'display':float((low+high)/2),'width':str(high-low)}
total=enclosure(own,claimed['actual_Gamma_deletion_area']);assert F(total['upper'])<F(1,125)<F(1,8)
groups={}
for name,source in [('pure_old',[n for n in own if n['old_count']and n['outside_mask']==0]),('old_plus_one_outside',[n for n in own if n['old_count']and n['outside_mask']in(1,2)]),('old_plus_both_outside',[n for n in own if n['old_count']and n['outside_mask']==3])]:
    groups[name]=enclosure(source,claimed['groups_individually_after_complete_residual'][name])
assert F(groups['pure_old']['lower'])==F(groups['pure_old']['upper'])==F(81,160000)
reported_joint={(r['k'],r['l']):r for r in claimed['joint_label_classification']};individual=[]
for node in own:
    if node['old_count']and node['outside_mask']==3 and node['Gamma_strict_active']:
        other=reported_joint[node['k'],node['l']]
        interval=enclosure([node],other['area_outside_complete_residual']);assert F(interval['lower'])>0
        assert other['classification']=='positive_surviving_area'
        individual.append({'k':node['k'],'l':node['l'],'surviving_area':interval})
assert len(individual)==len(reported_joint)==9
active={str(t):sum(n['Gamma_strict_active']for n in own if n['old_count']and n['outside_mask']==t)for t in range(4)}
assert active==claimed['active_counts_by_outside_mask']=={'0':1,'1':22,'2':0,'3':9}
out={'status':'PASS','method':'actual raw402atom old-prefix extraction, own label bbox/mass/radius construction, exact integer horizontal y-slabs and merged xinterval differences; no producer imports',
 'complete_menu_reconstructed':403,'old_prefix_times_outside_labels':400,'residual_labels':3,'catalogue_field_comparisons':field_checks,
 'actual_deletion_Gamma':total,'group_unions_individually_after_complete_residual':groups,'joint_individual_survivors':individual,
 'active_old_label_counts':active,'all_nine_joint_labels_strictly_survive_residual':True,'horizontal_slabs_executed':slab_count,
 'all_producer_area_endpoints_matched':True,'no_recomputed_capacity_catalogue':True,
 'scope':'One unchanged qualified S20 input, actual Gamma only, no collapsed old mass. R77 caps reused by exact input/result pins. Nine positive individual joint survivors are not summed. No general area inequality or C2 bound.',
 'input_sha256':hashlib.sha256(ir).hexdigest(),'producer_area_source':str(pp),'producer_area_sha256':hashlib.sha256(pr).hexdigest(),
 'inherited_qualification_sha256':hashlib.sha256(qr).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[4]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('joint_individual_survivors',)}))
