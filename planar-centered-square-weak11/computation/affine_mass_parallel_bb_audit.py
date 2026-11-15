"""Independent tiny control-flow audit for the parallel R24 B&B.
Does not evaluate any real-family halo. It compares a 2-worker four-leaf run
against the previously independent full-box replay and checks every tree edge.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,subprocess,sys,time
start=time.monotonic();topic=Path(sys.argv[1]);baseline_path=Path(sys.argv[2]);temp=Path(sys.argv[3]);out=Path(sys.argv[4]);comp=topic/'computation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
baseline=json.loads(baseline_path.read_text());input_path=temp/'positive-family.json';data=json.loads(input_path.read_text())
ref=data['base'];refinput=temp/'parallel-reference-input.json';refinput.write_text(json.dumps(ref))
old=load('parallel_audit_old',comp/'periodic_halo.py');certificate=old.calculate(ref,bits=5,tail=F(2));certificate.update(input_sha256=sha(refinput),code_sha256=sha(comp/'periodic_halo.py'))
assert F(certificate['inner_ratio'])==F(baseline['tiny_branch_and_bound']['lower'])
refcert=temp/'parallel-reference-certificate.json';refcert.write_text(json.dumps(certificate))
output=temp/'parallel-two-worker-four-leaf.json';cmd=[sys.executable,'-B',str(comp/'parallel_branch_mass_family.py'),str(input_path),'--output',str(output),'--reference-input',str(refinput),'--reference-certificate',str(refcert),'--reference-point','0','0','--workers','2','--bits','5','--radius','2','--gap','1/10000000000','--max-leaves','4','--seconds','30','--checkpoint-seconds','0.001']
run=subprocess.run(cmd,capture_output=True,text=True);assert run.returncode==0,run.stderr
result=json.loads(output.read_text());nodes=result['nodes'];assert len(nodes)==7 and result['leaf_count']==4 and result['splits']==3
for i,n in enumerate(nodes):
 assert n['id']==i
 cell=tuple(map(F,n['cell']));a,b,c,d=cell;assert -1<=a<b<=1 and -1<=c<d<=1
 if n['status']=='split':
  assert len(n['children'])==2
  left,right=[nodes[j] for j in n['children']];assert left['parent']==right['parent']==i
  lcell,rcell=tuple(map(F,left['cell'])),tuple(map(F,right['cell']))
  am,bm=(a+b)/2,(c+d)/2
  assert (lcell,rcell) in [((a,am,c,d),(am,b,c,d)),((a,b,c,bm),(a,b,bm,d))]
 else:assert n['children']==[] and n['status'] in ['active','pruned']
 if n['parent'] is not None:assert F(n['upper'])<=F(nodes[n['parent']]['upper'])
leaves=[n for n in nodes if n['status']!='split'];assert sum((F(n['cell'][1])-F(n['cell'][0]))*(F(n['cell'][3])-F(n['cell'][2])) for n in leaves)==4
lookup={tuple(n['cell']):F(n['upper']) for n in leaves};expected={tuple(n['cell']):F(n['upper']) for n in baseline['tiny_branch_and_bound']['leaves']};assert lookup==expected
lower=F(result['lower']);assert all(F(n['upper'])<lower for n in leaves if n['status']=='pruned')
assert F(result['upper'])==max([lower]+[F(n['upper']) for n in leaves])
assert F(result['gap'])==F(result['upper'])-lower
if result['status']=='target_gap':assert F(result['gap'])<=F(result['target_gap'])
# Invalid reference point must fail before worker scheduling.
badcmd=cmd.copy();i=badcmd.index('--reference-point');badcmd[i+1:i+3]=['2','0'];badcmd[badcmd.index('--output')+1]=str(temp/'outside-reference-rejected.json')
bad=subprocess.run(badcmd,capture_output=True,text=True);assert bad.returncode!=0 and 'outside family' in bad.stderr
record={'scope':'parallel-driver tiny exact tree audit only; no real-family envelope evaluation',
 'python':sys.version,'elapsed_seconds':time.monotonic()-start,'script_sha256':sha(Path(__file__)),
 'source_sha256':sha(comp/'parallel_branch_mass_family.py'),'reference_point_exactly_bound':True,'outside_reference_point_rejected':True,
 'all_seven_node_edges_and_splits_checked':True,'all_four_leaves_match_independent_single_worker_certificate':True,
 'tree':result,'baseline_sha256':sha(baseline_path),'actual_command':cmd}
out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in ['tree','actual_command']},indent=2))
