from pathlib import Path
import json, argparse
p=argparse.ArgumentParser();p.add_argument("--out",type=Path);args=p.parse_args()
S={0,1,2,3,4,5,8}
base={'I':({0,1,2,3,4,5},{0}),'II':({0,1,2},{0,3}),'III':({0,1},{0,2,4})}
out={}
for name,(aa,bb) in base.items():
 for t in (0,1):
  apos=aa|({8}if t else set());bpos=bb|({8}if not t else set())
  au={r for r in S if all((r+j)%9 in S for j in bpos)};bu={r for r in S if all((r+i)%9 in S for i in apos)}
  assert apos<=au and bpos<=bu
  rows={str(k):[(i,j,(i+j)//9)for i in sorted(au)for j in sorted(bu)if(i+j)%9==k]for k in range(9)}
  out[name+str(t)]={'known_positive_A':sorted(apos),'known_positive_B':sorted(bpos),'upper_support_bound_A':sorted(au),'upper_support_bound_B':sorted(bu),'rows':rows,'warning':'Upper bounds are not assumptions of positivity; unknown rows must remain.'}
assert out['II0']['upper_support_bound_B']==[0,1,2,3,8]
assert out['II1']['upper_support_bound_B']==[0,1,2,3]
if args.out: args.out.write_text(json.dumps(out,indent=2)+'\n')
print({k:(v['upper_support_bound_A'],v['upper_support_bound_B'])for k,v in out.items()})
