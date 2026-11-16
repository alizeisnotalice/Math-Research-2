#!/usr/bin/env python3
"""One bounded exact full-parameter upper certificate, rational area only."""
from fractions import Fraction as F
from collections import defaultdict,Counter
from pathlib import Path
import argparse,gzip,hashlib,heapq,json,math,os,platform,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def B(a):return min(max(F(1),(a+(k-1)*F(4,9))**2)-F(k*k,9) for k in (1,2,3))/2
def encbox(box):return list(map(str,box))
def merge(iv):
 out=[]
 for a,b in sorted(iv):
  if not out or a>out[-1][1]:out.append([a,b])
  elif b>out[-1][1]:out[-1][1]=b
 return out
def intersection(a,b):
 i=j=0;z=F(0)
 while i<len(a) and j<len(b):
  z+=max(F(0),min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<=b[j][1]:i+=1
  else:j+=1
 return z
def area(full,res):
 events=defaultdict(list)
 for ch,rr in enumerate((full,res)):
  for j,(a,b,c,d) in enumerate(rr):events[a].append((ch,j,(c,d)));events[b].append((ch,j,None))
 if not events:return F(0)
 active=[{},{}];last=min(events);out=F(0)
 for x,ee in sorted(events.items()):
  f=merge(active[0].values());r=merge(active[1].values());out+=(x-last)*(sum((b-a for a,b in f),F(0))-intersection(f,r))
  for ch,j,value in ee:
   if value is None:del active[ch][j]
   else:active[ch][j]=value
  last=x
 assert not active[0] and not active[1];return out
def clip(rects):
 out=[]
 for a,b,c,d in rects:
  z=(max(F(-1,2),a),min(F(-1,9),b),max(F(-1,2),c),min(F(-1,9),d))
  if z[0]<z[1] and z[2]<z[3]:out.append(z)
 return out
def radius(m,D):
 assert m>=0;num=m.numerator*D*D;den=4*m.denominator;k=math.isqrt(num//den);assert den*k*k<=num<den*(k+1)*(k+1);h=k if den*k*k==num else k+1
 return F(k,D),F(h,D),{'mass':str(m),'lower_integer':k,'upper_integer':h}
def bound(box,D):
 aL,aU,bL,bU=box;m=F(1,9);_,u,ur=radius(m+bU,D);_,v,vr=radius(m+2*bU,D);s,_,sr=radius(bL,D);t,_,tr=radius(2*bL,D)
 full=clip([(F(-1,6),F(1,6),F(-1,6),F(1,6)),(-u,u-aL,-u,u),(-u,u,-u,u-aL),(-v,v-aL,-v,v-aL)])
 res=clip([(-aL-s,-aU+s,-s,s),(-s,s,-aL-s,-aU+s),(-t,t-aU,-t,t-aU)])
 upper=area(full,res);assert 0<=upper<=F(49,324)
 return upper,{'uU':ur,'vU':vr,'sL':sr,'tL':tr},[encbox(r) for r in full],[encbox(r) for r in res]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--theory',required=True);ap.add_argument('--out',required=True);ap.add_argument('--certificate',required=True);a=ap.parse_args();start=time.monotonic();ip=Path(a.input);params=json.loads(ip.read_text());tp=Path(a.theory);assert sha(tp)==params['theory_sha256'];assert params['n']==3;D=1<<params['bits'];target=F(params['target']);nodes=[];heap=[]
 def create(box,parent,depth):
  i=len(nodes);aL,aU,bL,bU=box;assert F(1,9)<=aL<=aU<=1 and 0<=bL<=bU<=F(4,9);cap=B(aU);rec={'id':i,'parent':parent,'depth':depth,'original_box':encbox(box),'cap_at_aU':str(cap),'removed_cap_strip':None};nodes.append(rec)
  if bL>cap:rec.update(status='cap_infeasible',kept_box=None);return i
  if bU>cap:
   rec['removed_cap_strip']={'a_interval':[str(aL),str(aU)],'b_interval':[str(cap),str(bU)],'b_lower_open':True,'b_upper_closed':True,'reason':'b>B3(aU)>=B3(a)'};bU=cap
  kept=(aL,aU,bL,bU);upper,roots,full,res=bound(kept,D);rec.update(kept_box=encbox(kept),root_bounds=roots,full_outer_clipped=full,residual_inner_clipped=res,area_upper=str(upper))
  if upper<=target:rec['status']='certified'
  else:rec['status']='pending';heapq.heappush(heap,(-float(upper),i))
  return i
 rootbox=tuple(map(F,params['root_box']));create(rootbox,None,0);budget_stop=False;split_count=0
 while heap:
  if len(nodes)+2>params['max_nodes']:budget_stop=True;break
  _,i=heapq.heappop(heap);r=nodes[i]
  if r['status']!='pending':continue
  box=tuple(map(F,r['kept_box']));aL,aU,bL,bU=box
  if r['depth']>=params['max_depth']:r['status']='unresolved_depth';continue
  if aL==aU and bL==bU:r['status']='unresolved_point';continue
  dim=0 if aU-aL>=bU-bL else 1;mid=(box[2*dim]+box[2*dim+1])/2;left=list(box);right=list(box);left[2*dim+1]=mid;right[2*dim]=mid;r.update(status='split',split_dimension='a' if dim==0 else'b',split_midpoint=str(mid));c1=create(tuple(left),i,r['depth']+1);c2=create(tuple(right),i,r['depth']+1);r['children']=[c1,c2];split_count+=1
  if split_count%500==0:print(json.dumps({'event':'progress','nodes':len(nodes),'pending':len(heap),'splits':split_count,'seconds':time.monotonic()-start}),flush=True)
 # Internal structural checks. Audit independently replays them and leaf bounds.
 seen=set();stack=[0]
 while stack:
  i=stack.pop();assert i not in seen;seen.add(i);r=nodes[i];box=tuple(map(F,r['original_box']));cap=F(r['cap_at_aU']);assert cap==B(box[1])
  if r['status']=='cap_infeasible':assert box[2]>cap;continue
  kept=tuple(map(F,r['kept_box']));assert kept==(box[0],box[1],box[2],min(box[3],cap))
  if r['status']=='split':
   dim=0 if r['split_dimension']=='a' else 1;mid=F(r['split_midpoint']);left=list(kept);right=list(kept);left[2*dim+1]=mid;right[2*dim]=mid;c1,c2=r['children'];assert nodes[c1]['parent']==nodes[c2]['parent']==i;assert nodes[c1]['original_box']==encbox(left) and nodes[c2]['original_box']==encbox(right);stack.extend([c1,c2])
  elif r['status']=='certified':assert F(r['area_upper'])<=target
 assert len(seen)==len(nodes)
 counts=Counter(r['status'] for r in nodes);leaves=[r for r in nodes if r['status']!='split'];unresolved=[r['id'] for r in leaves if r['status'] not in ('certified','cap_infeasible')];proved=not unresolved;worst=max((F(r['area_upper']) for r in leaves if r['status']!='cap_infeasible'),default=F(0));cert={'schema':'r075-orthogonal-n3-full-parameter-tree-v1','parameters':params,'root_id':0,'coordinate_root_scale':D,'theory_path':str(tp),'theory_sha256':sha(tp),'nodes':nodes,'node_count':len(nodes),'leaf_count':len(leaves),'unresolved_leaf_ids':unresolved,'all_feasible_parameter_domain_covered_by_certified_leaves':proved,'a_below_1over9_degenerate_branch':{'a_interval':'0<a<1/9','admissible_b':'0','exact_Gamma_area':'1/324'},'scope':'Exact parameter-tree arithmetic. Closed children cover each retained parent box; each removed upper-b strip is infeasible by monotone B3. Uniform fullouter/resinner envelopes from frozen theory; no claim beyond specifiedn3 orthogonal equalmass family.'};cp=Path(a.certificate);cp.write_bytes(gzip.compress((json.dumps(cert,separators=(',',':'))+'\n').encode(),mtime=0));out={'schema':'r075-orthogonal-n3-box-summary-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'input_path':str(ip),'input_sha256':sha(ip),'theory_path':str(tp),'theory_sha256':sha(tp),'parameters':params,'certificate_path':str(cp),'certificate_sha256':sha(cp),'node_count':len(nodes),'leaf_count':len(leaves),'status_counts':dict(counts),'removed_cap_strip_count':sum(r['removed_cap_strip'] is not None for r in nodes),'maximum_depth':max(r['depth'] for r in nodes),'target_proved':proved,'unresolved_leaf_count':len(unresolved),'budget_stop':budget_stop,'maximum_terminal_area_upper':str(worst),'maximum_terminal_area_upper_decimal':float(worst),'target':'1/18','minimum_target_margin':str(target-worst),'strict_target_margin':proved and worst<target,'seconds':time.monotonic()-start,'scope':'One bounded parameter-box run only. All root/area/cap/coverage decisions exact; floating heap priority changes order only. No maximum-area optimization or generaloutside-family claim. Independent leaf and coverage audit recorded separately.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
