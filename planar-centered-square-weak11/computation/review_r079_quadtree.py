"""Complete saved-tree audit; no adaptive search and no producer imports."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import gzip,hashlib,json,math,os,platform,sys,time
start=time.monotonic();root=Path(sys.argv[1]);outpath=Path(sys.argv[2]);prod=root/'computation/results/round079-quad/production';parampath=root/'computation/inputs/r079-quad/parameters.json'
paramsraw=parampath.read_bytes();params=json.loads(paramsraw);summaryraw=(prod/'summary.json').read_bytes();summary=json.loads(summaryraw)
SOURCE='428c3c90516d9d51a062c602fefb85c8da722ae9db7abf331c1f0c5eb7645aee';PARAM='380edaf5c4fa12b9fa4cb58caab5fc3e8e5bc8d8fd0c7d7dc2fd42e0d5cfaf73'
assert hashlib.sha256(paramsraw).hexdigest()==PARAM
assert hashlib.sha256((root/'computation/round079-quad/quadratic_envelope.py').read_bytes()).hexdigest()==SOURCE
assert summary['source_sha256']==SOURCE and summary['parameters_sha256']==PARAM
VARIANTS=['exact_pureold_plus_three_mixed','four_soft_branches'];H=F(1,2);e=F(1,400);h=F(21,400);c=F(400,441);A=4-c;tol=F(params['target_area_gap'])
assert params['tiles_per_axis']==2 and params['residual_bits']==72 and tol==F(1,10000)
expected_inputs={'corrected':'d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a','farface':'e343e448b9582a731ddaccca37452adcce1c6f28b4153aa520a0ec2b4a3ef574'}
expected_certs={'corrected':'fb00faf83ebce3009e49ae0160d1c70669fcbcd4444060a72c9a532d31be8d78','farface':'f1568f9afa41688a2a3f68705024d25874c4f196e8c4b52f541db30fbe90cc50'}
auditbase=Path(__file__).resolve().parents[1]
priorfiles={'corrected':('actual-prefix-replay.json','41313173bf5345ea714da183414df417a09617e9879d27abde4816980ab79718'),'farface':('farface-replay.json','17746df9254b3dab70b05880581b6ee55e2b86cd89cd7718b169825b20ea4a5e')}
counts=Counter();case_results=[];task_pids=set()
def area(r):return (r[1]-r[0])*(r[3]-r[2])if r else F(0)
def intersect(r,s):
    if r is None or s is None:return None
    z=(max(r[0],s[0]),min(r[1],s[1]),max(r[2],s[2]),min(r[3],s[3]))
    return z if z[0]<z[1]and z[2]<z[3]else None
def remaining(region,rects):
    """Horizontal strips of clipped residual rectangles, no inclusion-exclusion."""
    counts['residual_area_calls']+=1
    if region is None:return F(0)
    rr=[z for r in rects if (z:=intersect(region,r))is not None]
    if not rr:return area(region)
    yy=sorted({region[2],region[3]}|{v for r in rr for v in r[2:]});ans=F(0)
    for low,high in zip(yy,yy[1:]):
        y=(low+high)/2;iv=sorted((a,b)for a,b,c0,d0 in rr if c0<y<d0);merged=[]
        for a,b in iv:
            if not merged or a>merged[-1][1]:merged.append([a,b])
            else:merged[-1][1]=max(merged[-1][1],b)
        free=region[1]-region[0]-sum((b-a for a,b in merged),F(0));assert free>=0
        ans+=(high-low)*free;counts['residual_horizontal_slabs']+=1
    assert 0<=ans<=area(region);return ans
def bbox(points):return(min(x for x,y in points),max(x for x,y in points),min(y for x,y in points),max(y for x,y in points))
def radius_range(cell,box):
    xl,xu,yl,yu=cell;a,b,d,f=box
    lo=max((b-a)/2,(f-d)/2,b-xu,xl-a,f-yu,yl-d,F(0))
    hi=max(xu-a,b-xl,yu-d,f-yl,F(0))
    assert max(xl,yl)<=lo<=hi and hi>=max(xu,yu)
    return lo,hi
def minimum(x,y,lo):
    """Use derivative signs instead of the producer's vertex clamp."""
    assert max(x,y)<=lo<=H
    beta=c*(x+y-2*h);constant=-c*(x-h)*(y-h)
    if 2*A*lo+beta>=0:result=A*lo*lo+beta*lo+constant
    elif A+beta<=0:result=A/4+beta/2+constant
    else:result=constant-beta*beta/(4*A)
    counts['quadratic_minima']+=1;return result
def branch(cell,box,mass):
    lo,hi=radius_range(cell,box);counts['branch_tests']+=1
    if lo>=H:return 'out'
    v=minimum(cell[0],cell[2],lo)
    if v==mass:counts['lower_polynomial_equalities']+=1
    if v>=mass:return 'out'
    if hi<H and minimum(cell[1],cell[3],hi)<mass:return 'in'
    return 'unknown'
def union(states):return 'in'if'in'in states else'out'if all(t=='out'for t in states)else'unknown'
def expected_leaf(cell,branches,ri,ro,oldrect,path):
    ai=remaining(cell,ri);ao=remaining(cell,ro);assert ao<=ai
    if ai==0:return {'path':path,'branch_states':['skipped_residual_cover']*4,'mixed_state':'residual_covered','soft_state':'residual_covered','remaining_residual_inner':ai,'remaining_residual_outer':ao,'known_pureold_remaining_lower':F(0),'known_pureold_remaining_upper':F(0),'lower':[F(0)]*2,'upper':[F(0)]*2,'unknown_area':[F(0)]*2}
    states=[branch(cell,b,m)for b,m in branches];mixed=union(states[1:]);soft=union(states)
    old=intersect(cell,oldrect);oldlo=remaining(old,ro);oldhi=remaining(old,ri)
    if mixed=='in':exactlo,exacthi=ao,ai
    elif mixed=='out':exactlo,exacthi=oldlo,oldhi
    else:exactlo,exacthi=oldlo,ai
    lower=[exactlo,ao if soft=='in'else F(0)];upper=[exacthi,F(0)if soft=='out'else ai]
    assert all(0<=a<=b<=area(cell)for a,b in zip(lower,upper))
    return {'path':path,'branch_states':states,'mixed_state':mixed,'soft_state':soft,'remaining_residual_inner':ai,'remaining_residual_outer':ao,'known_pureold_remaining_lower':oldlo,'known_pureold_remaining_upper':oldhi,'lower':lower,'upper':upper,'unknown_area':[area(cell)if mixed=='unknown'else F(0),area(cell)if soft=='unknown'else F(0)]}
def compare(value,record):
    counts['leaf_fields_checked']+=1
    if isinstance(value,F):assert value==F(record)
    elif isinstance(value,list):
        assert len(value)==len(record)
        for a,b in zip(value,record):compare(a,b)
    else:assert value==record
def children(box):
    xl,xu,yl,yu=box;xm=(xl+xu)/2;ym=(yl+yu)/2
    return[(xl,xm,yl,ym),(xm,xu,yl,ym),(xl,xm,ym,yu),(xm,xu,ym,yu)]
for cp in params['cases']:
    cid=cp['id'];ip=parampath.parent/cp['input'];ir=ip.read_bytes();assert hashlib.sha256(ir).hexdigest()==expected_inputs[cid]
    inp=json.loads(ir);pts=[tuple(map(F,p))for p in inp['positions']];masses=list(map(F,inp['masses']))
    assert inp['S_indices']==list(range(400))and pts[:400]==[(i*h,j*h)for i in range(20)for j in range(20)]and masses[:400]==[e]*400
    assert len(pts)==len(masses)==402 and F(inp['threshold'])==1
    pf,psha=priorfiles[cid];previousraw=(auditbase/'computation/results/round078-independent'/pf).read_bytes();assert hashlib.sha256(previousraw).hexdigest()==psha
    previous=json.loads(previousraw);assert previous['status']=='PASS'and previous['input_sha256']==expected_inputs[cid]
    reflected=[(-x,-y)for x,y in pts[400:]];branches=[]
    for tag in range(4):
        jj=[j for j in range(2)if tag>>j&1];branches.append((bbox([(F(0),F(0))]+[reflected[j]for j in jj]),sum((masses[400+j]for j in jj),F(0))))
    sm=next(z for z in summary['cases']if z['case_id']==cid);certpath=prod/sm['certificate'];certbytes=certpath.read_bytes();assert hashlib.sha256(certbytes).hexdigest()==sm['certificate_sha256']==expected_certs[cid]
    cert=json.loads(gzip.decompress(certbytes));assert cert['case_id']==cid and cert['variants']==VARIANTS and cert['input_sha256']==expected_inputs[cid] and cert['source_sha256']==SOURCE and cert['parameters_sha256']==PARAM
    assert cert['n']==20 and F(cert['h'])==h and F(cert['c'])==c and cert['residual_bits']==params['residual_bits']
    assert list(map(F,cert['positive_Gamma']))==[e,H,e,H]
    assert cert['children_order']==['lowerX_lowerY','upperX_lowerY','lowerX_upperY','upperX_upperY']
    oldrect=(e,F(1,40),e,F(1,40));assert list(map(F,cert['pureold_rectangle']))==list(oldrect)
    denom=math.lcm(400,*[v.denominator for p in reflected for v in p],*[m.denominator for m in masses[400:]]);scale=denom*(1<<cert['residual_bits']);ri=[];ro=[]
    cat={r['outside_mask']:r for r in cert['residual_catalog']};assert set(cat)=={1,2,3}
    for tag in (1,2,3):
        ids=[j for j in range(2)if tag>>j&1];mass=sum((masses[400+j]for j in ids),F(0));box=bbox([reflected[j]for j in ids]);span=max(box[1]-box[0],box[3]-box[2]);saved=cat[tag]
        assert F(saved['mass'])==mass and list(map(F,saved['bbox']))==list(box)and saved['eligible']==(mass>span*span)
        if mass<=span*span:continue
        square=mass*scale*scale/4;lo=math.isqrt(square.numerator//square.denominator);hi=lo if lo*lo==square else lo+1
        assert lo*lo<=square<(lo+1)*(lo+1)and square<=hi*hi
        low,high=F(lo,scale),F(hi,scale);assert high<=H
        assert F(saved['radius_lower'])==low and F(saved['radius_upper'])==high
        a,b,d,f=box;inner=(b-low,a+low,f-low,d+low);outer=(b-high,a+high,f-high,d+high)
        assert list(map(F,saved['inner']))==list(inner)and list(map(F,saved['outer']))==list(outer)
        ri.append(inner);ro.append(outer);counts['residual_root_bounds']+=1
    assert len(cert['tiles'])==4
    expected_roots={f'{i}-{j}':(e+(H-e)*i/2,e+(H-e)*(i+1)/2,e+(H-e)*j/2,e+(H-e)*(j+1)/2)for i in range(2)for j in range(2)}
    assert {t['tile_id']for t in cert['tiles']}==set(expected_roots)
    tiles=[];case_low=[F(0),F(0)];case_high=[F(0),F(0)];case_unknown=[F(0),F(0)];case_states=Counter()
    for tile in cert['tiles']:
        rootbox=expected_roots[tile['tile_id']];assert tuple(map(F,tile['root_box']))==rootbox and tile['case_id']==cid
        trie={};records=tile['leaf_records'];paths=set()
        for rec in records:
            path=rec['path'];assert path not in paths and all(ch in'0123'for ch in path);paths.add(path);d=trie
            for ch in path:
                assert '_leaf'not in d;d=d.setdefault(ch,{})
            assert not d;d['_leaf']=rec
        stack=[('',rootbox,trie)];low=[F(0),F(0)];high=[F(0),F(0)];unknown=[F(0),F(0)];area_sum=F(0);leaf_count=split_count=depth=0
        while stack:
            path,cell,tr=stack.pop();depth=max(depth,len(path))
            if '_leaf'in tr:
                assert len(tr)==1;rec=tr['_leaf'];want=expected_leaf(cell,branches,ri,ro,oldrect,path)
                assert set(rec)==set(want)
                for key,value in want.items():compare(value,rec[key])
                for i in range(2):low[i]+=want['lower'][i];high[i]+=want['upper'][i];unknown[i]+=want['unknown_area'][i]
                for k,state in enumerate(want['branch_states']):case_states[f'branch{k}_{state}']+=1
                area_sum+=area(cell);leaf_count+=1;counts['all_leaves']+=1
            else:
                assert set(tr)==set('0123');split_count+=1
                for j,sub in enumerate(children(cell)):stack.append((path+str(j),sub,tr[str(j)]))
        assert area_sum==area(rootbox)and leaf_count==len(records)==tile['leaf_count']==1+3*split_count
        assert tile['split_count']==split_count and tile['created_node_count']==1+4*split_count and tile['maximum_depth']==depth<=params['maximum_depth']
        assert 1+4*split_count<=params['maximum_nodes_per_tile']
        assert low==list(map(F,tile['lower']))and high==list(map(F,tile['upper']))and unknown==list(map(F,tile['unresolved_leaf_area']))
        gaps=[b-a for a,b in zip(low,high)];assert gaps==list(map(F,tile['gaps']))and F(tile['target_per_variant'])==tol/4
        if tile['termination']=='target_met':assert max(gaps)<=tol/4
        elif tile['termination']=='maximum_nodes':assert max(gaps)>tol/4 and tile['created_node_count']+4>params['maximum_nodes_per_tile']
        elif tile['termination']=='maximum_depth':assert max(gaps)>tol/4 and depth==params['maximum_depth']
        else:assert tile['termination']=='no_refinable_leaf'and max(gaps)>tol/4
        for i in range(2):case_low[i]+=low[i];case_high[i]+=high[i];case_unknown[i]+=unknown[i]
        task_pids.add(tile['pid']);tiles.append({'tile_id':tile['tile_id'],'leaf_count':leaf_count,'splits':split_count,'depth':depth,'termination':tile['termination'],'variant_gaps':list(map(str,gaps))})
        print(json.dumps({'progress':cid,'tile':tile['tile_id'],'verified_leaves':leaf_count,'seconds':time.monotonic()-start}),flush=True)
    assert sum(t['leaf_count']for t in tiles)==sm['leaf_count'] and sum(1+4*t['splits']for t in tiles)==sm['created_node_count']
    assert sm['all_tiles_target_met']==all(t['termination']=='target_met'for t in tiles)
    vs={}
    for i,name in enumerate(VARIANTS):
        saved=sm['variants'][name];gap=case_high[i]-case_low[i]
        assert F(saved['envelope_minus_true_residual_lower'])==case_low[i]and F(saved['envelope_minus_true_residual_upper'])==case_high[i]and F(saved['gap'])==gap
        assert F(saved['unresolved_leaf_area'])==case_unknown[i]and saved['target_met']==(gap<=tol)
        assert saved['actual_deletion_upper_valid']is True and saved['this_lower_is_not_actual_deletion_lower']is True
        vs[name]={'state_lower':str(case_low[i]),'state_upper':str(case_high[i]),'gap':str(gap),'case_target_met':gap<=tol,'unresolved_leaf_area':str(case_unknown[i]),'actual_deletion_lower_claimed':False}
    case_results.append({'case_id':cid,'leaf_count':sm['leaf_count'],'tiles':tiles,'branch_states':dict(case_states),'variants':vs,'input_sha256':expected_inputs[cid],'certificate_sha256':expected_certs[cid],'inherited_qualification_result_sha256':psha})
assert counts['all_leaves']==68912 and len(task_pids)==5 and sorted(task_pids)==summary['worker_pids']and summary['workers']==8
out={'status':'PASS','method':'complete prefix-trie coverage; rawinput/reflectedbbox/rational derivative-minimum cell tests; residual roots squared and horizontal residual-union subtraction; all leaf fields and all variant area sums rebuilt without producer imports or search',
 'cases':case_results,'counts':dict(counts),'configured_production_workers':8,'recorded_production_task_workers':5,
 'scope':'Two fixed uniform-grid inputs, continuous quadratic OUTER-state areas minus identical actual residual. Lower endpoints are NOT actual halo/deletion lower bounds. Unknown/node-budget leaves remain fully in upper; no tree regenerated, old caps/actual halo geometry not rerun.',
 'producer_source_sha256':SOURCE,'parameters_sha256':PARAM,'producer_summary_sha256':hashlib.sha256(summaryraw).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','counts':dict(counts),'cases':[{k:v for k,v in z.items()if k not in('branch_states','tiles')}for z in case_results],'runtime':out['runtime']}),flush=True)
