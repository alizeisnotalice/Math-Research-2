"""Audit archived .83745 leaf receipts and transfer margins. Run with -I and assertions enabled."""
import argparse,collections,hashlib,io,json,pathlib,sys,tarfile
from fractions import Fraction as F
import numpy as np
P=argparse.ArgumentParser();P.add_argument('package',type=pathlib.Path);P.add_argument('comp',type=pathlib.Path);P.add_argument('tree',type=pathlib.Path);P.add_argument('output',type=pathlib.Path);a=P.parse_args();assert __debug__
sys.path.insert(0,str(a.comp.resolve()));sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
from certify_tomar_six_high import SourceTree,FINE
from verify_five_selected_083740_traces import audit
sha=lambda x:hashlib.sha256(x).hexdigest()
def bundle(name):
 with tarfile.open(a.package/name,'r:gz') as t:
  rows=json.load(t.extractfile('MANIFEST.json'));D={}
  assert len(rows)+1==len(t.getnames())
  for n,z,h in rows:
   raw=t.extractfile(n).read();assert len(raw)==z and sha(raw)==h and n not in D;D[n]=raw
 return D
ordinary=bundle('native45-receipts.tar.gz');D={}
for name in ['native45-screen-1.tar.gz','native45-screen-2.tar.gz','native45-screen-3.tar.gz','native45-screen-4.tar.gz','native45-pure-a.tar.gz','native45-pure-b.tar.gz']:
 q=bundle(name);assert not D.keys()&q.keys();D.update(q)
j=lambda data:json.loads(data);scope=j(ordinary['inputs/ordinary/scope.json']); idx=dict(j(ordinary['inputs/ordinary/stage4-parts.json']));assert len(idx)==1024 and sha(ordinary['inputs/ordinary/stage4-parts.json'])==scope['parent_parts_sha256']
O={tuple(x) for x in scope['leaves']}; S={tuple(x) for x in scope['deferred']};assert len(O)==3071 and len(S)==202 and not O&S
T=SourceTree(a.tree); densekeys={(176,35174),(176,35426),(176,35428),(180,23933),(180,23934),(180,24225)};assert densekeys<=S
remaining=set(); layers=[];goal=120689264254125501
for layer in range(1,5):
 c=collections.Counter()
 for part in idx:
  prefix=f'screen/stage{layer}/part{part:04}';raw=D[prefix+'.json'];r=j(raw); cfg=r['config']; arr=D[prefix+'.npz'];assert cfg['part']==part and cfg['layer']==layer and cfg['goal']==goal and cfg['radial_goal']=='16749/20000' and sha(arr)==r['data_sha256']
  if layer>1:
   prev=D[f'screen/stage{layer-1}/part{part:04}.npz'];assert cfg['source_sha256']==sha(prev)
  q=np.load(io.BytesIO(arr),allow_pickle=False);A=list(map(int,q['accepted']));R=list(map(int,q['remaining']));assert len(A)==r['accepted'] and len(R)==r['remaining'] and len(A)+len(R)==r['tested'] and len(set(A))==len(A) and len(set(R))==len(R) and not set(A)&set(R) and all(int(z)>=goal for z in q['full_floor'])
  if layer>1:assert set(np.load(io.BytesIO(prev),allow_pickle=False)['remaining'])==set(A)|set(R)
  if layer==4:assert sha(raw)==idx[part];remaining.update((part,pos) for pos in R)
  for k in ('tested','accepted','remaining'):c[k]+=r[k]
 layers.append(dict(c))
assert remaining==O|S
seen=set();failed=set();ordinary_visits=0
for name in sorted(ordinary):
 if '/rows/' not in name:continue
 r=j(ordinary[name]);key=(r['part'],r['pos']);assert key not in seen;seen.add(key);t=r['stat'];ordinary_visits+=t[1]
 if t[0]:assert t[5]==0 and t[1]==2*(t[2]+t[3])-1
 else:failed.add(key)
assert seen==O and len(failed)==18
group_union=set()
for name in sorted(n for n in ordinary if n.startswith('ordinary/') and n.endswith('/result.json')):
 group=name.split('/')[1]; result=j(ordinary[name]); number=group.rsplit('-',1)[1]; selection_raw=ordinary[f'inputs/ordinary/group-{number}.json']; selection=j(selection_raw); keys={tuple(z) for z in selection['leaves']}; assert not group_union&keys;group_union|=keys
 rows=[j(v) for n,v in ordinary.items() if n.startswith('ordinary/'+group+'/rows/')]; assert {(r['part'],r['pos']) for r in rows}==keys and result['trees']==len(keys) and result['selection_sha256']==sha(selection_raw) and {tuple(x) for x in result['failed']}=={(r['part'],r['pos']) for r in rows if not r['stat'][0]}
assert group_union==O
purekeys=set();pure_visits=0
for name in sorted(D):
 if not name.startswith('pure/actual/0.837455/') or not name.endswith('.json'):continue
 r=j(D[name]);part,pos=r['part'],r['pos'];key=(part,pos);assert key not in purekeys;purekeys.add(key);assert all(sha(D['pure/'+k])==v for k,v in r['sources'].items()) and F(r['threshold'])==F('.837455') and not r['result']['remaining'] and r['original_header_sha256']==T.pins[f'part{part:04}.json'] and r['original_tree_sha256']==T.pins[f'part{part:04}.npy']
 box=T.leaf(part,pos);root=[int(z*2**44) for k in (0,1) for z in (v[k] for v in box)];assert root==r['root'];trace=D[name[:-5]+'.bin'];assert sha(trace)==r['trace_sha256'];out=audit(root,trace);assert out=={k:r['result'][k] for k in out};pure_visits+=out['vis']
assert purekeys==S-densekeys and len(purekeys)==196
all_dense=[]
for name,raw in D.items():
 if name.startswith('dense/proof-') and name.endswith('/result.json'):all_dense.append(('original',name,raw,D))
for name,raw in ordinary.items():
 if name.startswith('dense/') and name.endswith('/result.json'):all_dense.append(('rescue',name,raw,ordinary))
assert len(all_dense)==24
seen_dense=set();dense_visits=0
for category,name,raw,source in all_dense:
 x=j(raw);key=(x['part'],x['pos']); selection_name=('dense/' if category=='original' else 'inputs/dense/')+str(key[0])+'-'+str(key[1])+'.json';selection_raw=source[selection_name];selection=j(selection_raw);assert sha(selection_raw)==x['selection_sha256'] and selection['paths']==[''] and selection['terminal_guard']=='mesh64-dense';assert key not in seen_dense;seen_dense.add(key);assert x['all_closed'] and x['complete'] and x['threshold_integer']==goal and F(x['target'])==F('.83745') and len(x['rows'])==1
 box=T.leaf(*key);lo,hi=[[int(z*FINE) for z in row] for row in zip(*box)];assert (lo,hi)==(x['original_lo'],x['original_hi']);r=x['rows'][0];s=r['stat'];assert r['bits']=='' and r['terminal_guard']=='mesh64-dense' and r['lo']==lo and r['hi']==hi and s[0]==1 and s[5]==0 and s[1]==2*(s[2]+s[3])-1;dense_visits+=s[1]
assert seen_dense==densekeys|failed
# R208 wall holds conservatively over every original leaf.
mx=max(sum(max(abs(lo),abs(hi))**2 for lo,hi in T.leaf(*key)[:2]) for key in O|S);assert mx<F('.02816')**2
G=F(goal,2**57);scale=F('.999999006');rad=F(1,2);R=F('.703');h=F('.0125');m=F('.000000860')*F('.528619')/(1-scale);v=(1+scale)*rad*rad+(1-scale)*m*m;loss=F(22,7)*(2*R*F('.6')*(1-scale)+(F('.6')*(1-scale))**2)
shallow=scale*scale*F('.837295')+scale*rad*h*h/F('.495')-G;cap=4*(F('.47184')-m)*scale*scale-((1+scale)*h+max(F(0),(v-2*m*(rad+h))/(2*rad+h)));pure=F('.837455')-loss-G;factor=F(500000000377835707,499999570000000000);bound=factor*factor*G
assert min(shallow,cap,pure)>0 and bound>F('.83745')
out={'complete':True,'meaning':'receipt, branch and rational transfer checks; does not rerun area kernels or certify source-tree outside exclusion','layers':layers,'ordinary_closed':len(O)-len(failed),'ordinary_rescued':len(failed),'five_only':len(purekeys),'dense_original':len(densekeys),'ordinary_nodes':ordinary_visits,'pure_nodes':pure_visits,'dense_nodes':dense_visits,'max_translation_square':str(mx),'shallow_margin':str(shallow),'cap_margin':str(cap),'pure_margin':str(pure),'conditional_bound':str(bound)}
a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
