#!/usr/bin/env python3
"""Exact prefix caps and independent complete bbox qualification, no halo."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prefix(a,d,b,c):
 rows=[]
 def add(kind,k,l,m,s):
  rhs=max(160000,s*s);rows.append({'kind':kind,'k':k,'l':l,'mass_numerator':m,'side_units_400':s,'cap_numerator':rhs,'slack_numerator':rhs-m})
 for k in range(1,21):
  add('left',k,20,b+8000*k,a+21*(k-1))
  add('bottom',20,k,c+8000*k,d+21*(k-1))
 for k in range(1,21):
  for l in range(1,21):add('both',k,l,b+c+400*k*l,max(a+21*(k-1),d+21*(l-1)))
 add('outside_left',0,0,b,0);add('outside_bottom',0,0,c,0);add('outside_pair',0,0,b+c,max(a,d))
 assert len(rows)==443 and min(r['slack_numerator'] for r in rows)>=0
 return {'denominator':160000,'rows':rows,'minimum_slack_numerator':min(r['slack_numerator'] for r in rows),'all_pass':True}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();t=time.monotonic();ip=Path(a.input);data=json.loads(ip.read_text());points=[tuple(map(F,p)) for p in data['positions']];m=list(map(F,data['masses']));assert len(points)==402
 p=[tuple(int(x*400) for x in q) for q in points];assert all(F(z,400)==x for q,r in zip(points,p) for x,z in zip(q,r));mass=[int(x*160000) for x in m];assert all(F(z,160000)==x for x,z in zip(m,mass))
 old=(1<<400)-1;q1=1<<400;q2=1<<401;xx=sorted({x for x,y in p});yy=sorted({y for x,y in p})
 def groups(axis,coords):return [sum(1<<i for i,q in enumerate(p) if q[axis]==v) for v in coords]
 xg,yg=groups(0,xx),groups(1,yy)
 def intervals(g):
  out=[]
  for i in range(len(g)):
   z=0
   for j in range(i,len(g)):z|=g[j];out.append(z)
  return out
 closures={x&y for x in intervals(xg) for y in intervals(yg)}-{0};maxmass=0;maxmasks=[];eligible=0;mincap=10**100;flags=set();capsamples=[]
 for mask in sorted(closures):
  xmin=next(x for x,g in zip(xx,xg) if mask&g);xmax=next(x for x,g in zip(reversed(xx),reversed(xg)) if mask&g)
  ymin=next(y for y,g in zip(yy,yg) if mask&g);ymax=next(y for y,g in zip(reversed(yy),reversed(yg)) if mask&g)
  num=400*(mask&old).bit_count()+40000*bool(mask&q1)+80000*bool(mask&q2)
  side=max(xmax-xmin,ymax-ymin);slack=max(160000,side*side)-num;assert slack>=0;mincap=min(mincap,slack)
  if num>side*side:
   eligible+=1;f=int(bool(mask&old))+2*bool(mask&q1)+4*bool(mask&q2);flags.add(f)
   if num>maxmass:maxmass=num;maxmasks=[mask]
   elif num==maxmass:maxmasks.append(mask)
 assert maxmass==160000 and maxmasks==[old]
 # Old atoms are connected by the eligible S block. Two explicit eligible
 # edges attach both outside atoms, so all 402 sources are connected.
 edges=[]
 for i,j in [(0,400),(400,401)]:
  span=max(abs(p[i][0]-p[j][0]),abs(p[i][1]-p[j][1]));num=mass[i]+mass[j];assert num>span*span
  edges.append({'indices':[i,j],'mass_numerator':num,'span_squared_numerator':span*span,'strict_margin_numerator':num-span*span})
 fullspan=max(max(xx)-min(xx),max(yy)-min(yy));total=sum(m);gap=F(fullspan*fullspan,160000)-F(35,33)*total
 assert gap>0
 menu=[]
 for k in range(1,11):
  for l in range(1,11):
   xmax=21*(k-1);ymax=21*(l-1);twodelta=max(xmax+142,ymax+342,2*(xmax+1),2*(ymax+1),2*(-200+142),2*(-200+342),0)
   fullnum=120000+400*k*l;credit=twodelta*twodelta;menu.append({'k':k,'l':l,'old_mass':str(F(k*l,400)),'twice_delta_units_400':twodelta,'full_mass_numerator':fullnum,'four_delta_squared_numerator':credit,'strict_activity':fullnum>credit,'activity_margin_numerator':fullnum-credit})
 # Separately verify exact root bounds behind the short hand formula.
 assert F(433,250)**2<3 and F(283,200)**2>2
 handlower=F(43,200)*F(433,250)-F(101,1600)*F(283,200)-F(24641,160000);assert handlower==F(206483,1600000)>F(1,8)
 # Full synthetic measure after collapsing P (not merely the 3-source proxy)
 # violates a unit-square cap; record the position-loss mechanism exactly.
 window=(F(-71,200),F(129,200),F(-171,200),F(29,200));captured=[i for i,(x,y) in enumerate(points[:400]) if window[0]<=x<=window[1] and window[2]<=y<=window[3]]
 Pin=[i for i in captured if points[i][0]<F(199,400) and points[i][1]<F(199,400)]
 assert len(captured)==39 and len(Pin)==30
 out={'schema':'r077-n20-qualification-v1','pid':os.getpid(),'python':sys.version,'input_sha256':sha(ip),'script_sha256':sha(Path(__file__)),'integer_coordinate_denominator':400,'integer_mass_cap_denominator':160000,'prefix_original':prefix(124,332,20000,80000),'prefix_corrected':prefix(142,342,40000,80000),'actual_402_complete_bbox':{'axis_coordinate_counts':[len(xx),len(yy)],'closure_count':len(closures),'eligible_count':eligible,'all_caps_pass':True,'minimum_cap_slack_numerator':mincap,'maximum_eligible_mass':'1','maximum_masks':[str(z) for z in maxmasks],'unique_S20':True,'eligible_connectivity_flags':sorted(flags),'connected':True,'connection_certificate':{'old_S_mass':'1','old_S_span_squared':'159201/160000','outside_edges':edges},'total_mass':str(total),'full_span_squared':str(F(fullspan*fullspan,160000)),'nonterminal_gap_span2_minus_35over33total':str(gap),'nonterminal':True},'actual_joint_prefix_menu':{'count':100,'strict_active_count':sum(r['strict_activity'] for r in menu),'null_count':sum(not r['strict_activity'] for r in menu),'rows':menu,'scope':'Full actual-old-prefix plus both-outside witness activity only; no residual subtraction or actual grid corner area claim.'},'proxy_hand_formula':'43sqrt(3)/200-101sqrt(2)/1600-24641/160000','proxy_simple_strict_lower':str(handlower),'collapse_cap_loss':{'unit_closed_square':list(map(str,window)),'original_old_atoms_captured':len(captured),'captured_old_in_P':len(Pin),'original_full_mass':str(F(39,400)+F(3,4)),'synthetic_full_mass_after_collapse':str(F(9,400)+F(1,4)+F(3,4))},'scope':'Original S20+two outside is legal, unique maximum S20, connected, nonterminal. Mixed collapse proxy violation is not actual S20 corner violation. No full 400-source halo evaluated.','seconds':time.monotonic()-t}
 Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['pid','actual_402_complete_bbox','proxy_simple_strict_lower','collapse_cap_loss','seconds']},indent=2));print('joint active',out['actual_joint_prefix_menu']['strict_active_count'])
if __name__=='__main__':main()
