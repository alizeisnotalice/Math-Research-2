"""Read-only exact contract/tree check of one atomic R24 B&B snapshot.
No per-node large-family envelope or halo is recomputed.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
start=time.process_time();snapshot_path,contract_path,parallel_audit_path,main_audit_path,refinput_path,refcert_path,out=map(Path,sys.argv[1:])
raw=snapshot_path.read_bytes();record=json.loads(raw);contract=json.loads(contract_path.read_text());parallel_audit=json.loads(parallel_audit_path.read_text());main_audit=json.loads(main_audit_path.read_text());reference=json.loads(refinput_path.read_text());certificate=json.loads(refcert_path.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert record['input_sha256']==contract['input_sha256']
assert record['driver_sha256']==parallel_audit['source_sha256']
expected_envelope=next(v for k,v in main_audit['sources'].items() if k.endswith('/mass_family_envelope.py'));assert record['envelope_sha256']==expected_envelope
assert record['reference_input_sha256']==sha(refinput_path)==certificate['input_sha256']
assert record['reference_certificate_sha256']==sha(refcert_path)
assert certificate['code_sha256']=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'
assert F(record['lower'])==F(certificate['inner_ratio'])
assert F(certificate['mass'])==F(contract['actual_constant_mass'])
assert F(certificate['inner_ratio'])==F(certificate['inner_area'])*F(reference['lam'])/F(certificate['mass'])
a0,a1=map(F,contract['bounds']['a']);b0,b1=map(F,contract['bounds']['b']);a,b=map(F,record['reference_point']);assert a0<=a<=a1 and b0<=b<=b1
assert len(reference['pts'])==len(contract['all_label_records'])
for row,label in zip(reference['pts'],contract['all_label_records']):
 x,y,w=map(F,row);assert [str(x),str(y)]==label['coordinate']
 assert w==F(label['intercept'])+a*F(label['u'])+b*F(label['v'])
assert record['family']['all_corner_labels_strictly_positive']
for key,source_key in [('mass','actual_constant_mass'),('minimum_corner_label_mass','minimum_corner_mass'),('uniform_kappa','uniform_kappa'),('delta','delta'),('uniform_tail_bound','uniform_tail')]:assert F(record['family'][key])==F(contract[source_key])
nodes=record['nodes'];assert nodes and nodes[0]['parent'] is None and list(map(F,nodes[0]['cell']))==[a0,a1,b0,b1]
all_children=[];splits=0
for i,n in enumerate(nodes):
 assert n['id']==i;al,ah,bl,bh=map(F,n['cell']);assert a0<=al<ah<=a1 and b0<=bl<bh<=b1
 assert n['status'] in ['split','active','pruned']
 if i:
  assert type(n['parent'])is int and 0<=n['parent']<i
  assert F(n['upper'])<=F(nodes[n['parent']]['upper'])
 if n['status']=='split':
  splits+=1;assert len(n['children'])==2 and len(set(n['children']))==2
  left,right=[nodes[j] for j in n['children']];assert left['parent']==right['parent']==i
  l,r=tuple(map(F,left['cell'])),tuple(map(F,right['cell']));am,bm=(al+ah)/2,(bl+bh)/2
  assert (l,r) in [((al,am,bl,bh),(am,ah,bl,bh)),((al,ah,bl,bm),(al,ah,bm,bh))]
  all_children+=n['children']
 else:assert n['children']==[]
assert sorted(all_children)==list(range(1,len(nodes))) and len(nodes)==2*splits+1
leaves=[n for n in nodes if n['status']!='split'];assert len(leaves)==splits+1==record['leaf_count'] and splits==record['splits']
assert len(leaves)<=record['max_leaves']
assert sum((F(n['cell'][1])-F(n['cell'][0]))*(F(n['cell'][3])-F(n['cell'][2])) for n in leaves)==(a1-a0)*(b1-b0)
# Exact midpoint child unions plus an acyclic complete tree prove full coverage
# and interior disjointness; area equality is a separate arithmetic check.
lower=F(record['lower']);assert all(F(n['upper'])<lower for n in leaves if n['status']=='pruned')
assert all(F(n['upper'])>=lower for n in leaves if n['status']=='active')
upper=max([lower]+[F(n['upper']) for n in leaves]);assert upper==F(record['upper']) and upper-lower==F(record['gap'])
if record['status']=='target_gap':assert upper-lower<=F(record['target_gap'])
result={'scope':'one saved running/final checkpoint contract and complete tree; no per-node large geometry replay',
 'snapshot_sha256':hashlib.sha256(raw).hexdigest(),'snapshot_status':record['status'],'driver_sha256':record['driver_sha256'],'envelope_sha256':record['envelope_sha256'],
 'reference_point':record['reference_point'],'reference_all_labels_exactly_in_family':True,'reference_certificate_identity_and_mass_ratio_checked':True,
 'nodes':len(nodes),'splits':splits,'leaves':len(leaves),'all_binary_splits_and_backlinks_checked':True,'full_domain_cover_and_interior_disjointness_proved_by_tree':True,'sum_leaf_areas':str((a1-a0)*(b1-b0)),
 'lower':str(lower),'upper':str(upper),'gap':str(upper-lower),'cpu_seconds':time.process_time()-start,'script_sha256':sha(Path(__file__)),
 'source_paths':{str(p):sha(p) for p in [contract_path,parallel_audit_path,main_audit_path,refinput_path,refcert_path]}}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
