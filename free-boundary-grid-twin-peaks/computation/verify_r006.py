"""Exact certificate arithmetic and independent GPU witness checks."""
import json
from pathlib import Path
from fractions import Fraction
from decimal import Decimal
B=Path(__file__).parent
R=B/'results'
def read(p):return json.loads(p.read_text())
a=read(R/'n7-exact-r005.json');b=read(R/'n7-top2-r006.json')
for key in a:assert a[key]==b[key]
assert b['A1_from_top2_adj']==a['A1']
top=int(b['A2_top2_nonadj']);assert 0<top<int(a['A2'])
s=read(R/'top2-h-small-r006.json')
for n,row in s['results'].items():assert 2*sum(map(int,row['by_distance'].values()))==int(row['ordered_top2_sum'])
assert int(s['n5_selected_h']['opposite_corners_distance_8'])<int(s['n5_selected_h']['middle_straight_distance_2'])
c=read(R/'n8-census-r006.json');assert sum(c['histogram'].values())==c['independent_sets']
k,v=max(c['histogram'].items(),key=lambda x:x[1]);assert int(k)==15 and v==111494429212
assert 32*((v+7)//8)==445977716864
folder=B/'certificates/r006-gpu';moves=0;wall=Decimal(0);count=0;summary=[]
for f in sorted(folder.glob('gpu-search*.jsonl')):
 for line in f.read_text().splitlines():
  z=json.loads(line,parse_float=Decimal)
  if z.get('summary'):
   assert z['moves']==z['chains']*z['rounds']*z['steps']
   moves+=int(z['moves']);wall+=z['gpu_seconds'];summary.append(z['n']);continue
  n=z['n'];labels=z['labels'];assert sorted(labels)==list(range(1,n*n+1))
  def neighbors(v):
   x,y=v%n,v//n
   return [w for w,ok in ((v-1,x>0),(v+1,x<n-1),(v-n,y>0),(v+n,y<n-1)) if ok]
  peaks=[v for v in range(n*n) if all(labels[v]>labels[w] for w in neighbors(v))]
  assert len(peaks)==2
  u,v=peaks;assert abs(u%n-v%n)+abs(u//n-v//n)==z['distance']
  count+=1
assert moves==2348810240 and count==209 and len(summary)==8
out={'verified':True,'n7_top2_probability':str(Fraction(top,int(a['A2']))),'n8_optimistic_32byte_layer_lower_bound_bytes':445977716864,'gpu_proposals':moves,'gpu_saved_labelings_cpu_verified':count,'submission_wait_wall_seconds':str(wall),'hardware_gpu_seconds_measured_in_r006':False}
print(json.dumps(out,indent=2))
