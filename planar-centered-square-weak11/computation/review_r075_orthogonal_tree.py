"""Independent exact replay of a saved continuum cover, not a new B&B run."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import gzip,hashlib,json,os,platform,sys,time

start=time.monotonic(); cp=Path(sys.argv[1]); raw=cp.read_bytes(); data=json.loads(gzip.decompress(raw))
params=data['parameters']; assert params['n']==3 and params['bits']==72
assert list(map(F,params['root_box']))==[F(1,9),F(1),F(0),F(4,9)]
target=F(params['target']); assert target==F(1,18)
D=1<<72; assert data['coordinate_root_scale']==D
theory=Path(sys.argv[3]); theoryraw=theory.read_bytes(); theorysha=hashlib.sha256(theoryraw).hexdigest()
assert theorysha==data['theory_sha256']==params['theory_sha256']=='7a0e6c14b12f1ddb610e7c63f061400254f38863b3c314ee1528fad7e2553cf3'
def cap(a):
    return min(max(F(1),(a+F(4*(k-1),9))**2)-F(k*k,9) for k in range(1,4))/2
def clip(rr):
    out=[]
    for a,b,c,d in rr:
        z=(max(a,-F(1,2)),min(b,-F(1,9)),max(c,-F(1,2)),min(d,-F(1,9)))
        if z[0]<z[1] and z[2]<z[3]:out.append(z)
    return out
def merged(intervals):
    result=[]
    for a,b in sorted(intervals):
        if not result or a>result[-1][1]:result.append([a,b])
        else:result[-1][1]=max(result[-1][1],b)
    return result
def horizontal_area(full,res):
    yy=sorted({v for r in full+res for v in r[2:]}); total=F(0); count=0
    for y0,y1 in zip(yy,yy[1:]):
        y=(y0+y1)/2
        ff=merged([(a,b)for a,b,c,d in full if c<y<d])
        rr=merged([(a,b)for a,b,c,d in res if c<y<d])
        length=sum((b-a for a,b in ff),F(0))
        for a,b in ff:
            for c,d in rr:length-=max(F(0),min(b,d)-max(a,c))
        assert length>=0
        total+=(y1-y0)*length;count+=1
    return total,count
nodes={r['id']:r for r in data['nodes']};assert len(nodes)==len(data['nodes'])
root=data['root_id'];assert nodes[root]['parent'] is None and nodes[root]['depth']==0
assert list(map(F,nodes[root]['original_box']))==list(map(F,params['root_box']))
seen=set();stack=[root];counts=Counter();strips=roots=slabs=0;leaves=[];recorded_areas=[]
while stack:
    nid=stack.pop();assert nid not in seen;seen.add(nid);node=nodes[nid];counts[node['status']]+=1
    box=tuple(map(F,node['original_box']));aL,aU,bL,bU=box
    assert F(1,9)<=aL<=aU<=1 and 0<=bL<=bU<=F(4,9)
    bound=cap(aU);assert F(node['cap_at_aU'])==bound
    if bL>bound:
        assert node['status']=='cap_infeasible' and node['kept_box'] is None
        assert node['removed_cap_strip'] is None
        leaves.append(nid);continue
    kept=(aL,aU,bL,min(bU,bound)); assert tuple(map(F,node['kept_box']))==kept
    strip=node['removed_cap_strip']
    if bU>bound:
        assert strip is not None and list(map(F,strip['a_interval']))==[aL,aU]
        assert list(map(F,strip['b_interval']))==[bound,bU]
        assert strip['b_lower_open'] is True and strip['b_upper_closed'] is True
        strips+=1
    else:assert strip is None
    aL,aU,bL,bU=kept
    masses={'uU':F(1,9)+bU,'vU':F(1,9)+2*bU,'sL':bL,'tL':2*bL}
    radius={}
    for key,mass in masses.items():
        rec=node['root_bounds'][key];assert F(rec['mass'])==mass
        lo,hi=rec['lower_integer'],rec['upper_integer'];assert type(lo) is int and type(hi) is int and lo>=0
        scaled=mass*D*D
        assert 4*lo*lo<=scaled<4*(lo+1)*(lo+1)
        assert hi==(lo if 4*lo*lo==scaled else lo+1)
        assert scaled<=4*hi*hi
        radius[key]=F(hi if key in ('uU','vU') else lo,D);roots+=1
    u,v,s,t=(radius[k]for k in ('uU','vU','sL','tL'))
    full=clip([(-F(1,6),F(1,6),-F(1,6),F(1,6)),(-u,u-aL,-u,u),(-u,u,-u,u-aL),(-v,v-aL,-v,v-aL)])
    res=clip([(-aL-s,-aU+s,-s,s),(-s,s,-aL-s,-aU+s),(-t,t-aU,-t,t-aU)])
    assert full==[tuple(map(F,r))for r in node['full_outer_clipped']]
    assert res==[tuple(map(F,r))for r in node['residual_inner_clipped']]
    area,ns=horizontal_area(full,res);slabs+=ns
    assert area==F(node['area_upper']) and 0<=area<=F(49,324)
    recorded_areas.append(area)
    if node['status']=='split':
        dim={'a':0,'b':1}[node['split_dimension']];mid=F(node['split_midpoint'])
        assert kept[2*dim]<mid<kept[2*dim+1] and mid==(kept[2*dim]+kept[2*dim+1])/2
        left=list(kept);right=list(kept);left[2*dim+1]=mid;right[2*dim]=mid
        children=node['children'];assert len(children)==2 and children[0]!=children[1]
        for child,expected in zip(children,(left,right)):
            assert nodes[child]['parent']==nid and nodes[child]['depth']==node['depth']+1
            assert list(map(F,nodes[child]['original_box']))==expected
        stack.extend(children)
    else:
        assert node['status']=='certified' and area<=target
        leaves.append(nid)
assert seen==set(nodes) and data['node_count']==len(nodes) and data['leaf_count']==len(leaves)
assert data['unresolved_leaf_ids']==[] and data['all_feasible_parameter_domain_covered_by_certified_leaves'] is True
terminal=[F(nodes[i]['area_upper'])for i in leaves if nodes[i]['status']=='certified'];worst=max(terminal)
assert worst<target
assert data['a_below_1over9_degenerate_branch']=={'a_interval':'0<a<1/9','admissible_b':'0','exact_Gamma_area':'1/324'}
assert F(80,27)+F(1,18)+3*F(1,324)==F(109,36)<F(7,2)
out={'status':'PASS','method':'saved-tree traversal; exact cap rows and root squaring; reconstructed rectangles with horizontal y-slabs/merged xintervals, no producer imports or tree generation',
 'node_count':len(nodes),'leaf_count':len(leaves),'status_counts':dict(counts),'cap_clip_strips':strips,'root_bounds_checked':roots,
 'all_node_area_bounds_recomputed':len(recorded_areas),'horizontal_slabs':slabs,'maximum_depth':max(n['depth']for n in nodes.values()),
 'target':'1/18','maximum_certified_leaf_upper':str(worst),'minimum_target_margin':str(target-worst),
 'outside_root_domain':'0<a<1/9 forces b=0 and pureold Gamma area1/324; a=1/9 is included in root',
 'complete_Delta_upper':'109/36','complete_Delta_upper_scope':'full S3 equal orthogonal family,0<a<=1,0<=b<=B3(a); conditional only on frozen cap theorem and verified cover, no connected/nonterminal hypothesis',
 'certificate_source':str(cp),'certificate_sha256':hashlib.sha256(raw).hexdigest(),'theory_sha256':theorysha,
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
