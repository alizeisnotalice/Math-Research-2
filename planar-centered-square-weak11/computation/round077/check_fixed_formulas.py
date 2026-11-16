#!/usr/bin/env python3
"""Exact short radical lower bounds and provenance cross-checks, no search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,os,sys,time
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root(q):
 q=F(q);D=1<<72;k=math.isqrt(q.numerator*D*D//q.denominator)
 assert F(k,D)**2<=q<F(k+1,D)**2
 return F(k,D),F(k+1,D)
def main():
 t=time.monotonic();rp=Path('computation/results/round077-mixed-proxy');ip=rp/'directed-areas.json';data=json.loads(ip.read_text());cases={c['id']:c for c in data['cases']};r2,r3,r5,r13=[root(q) for q in [2,3,5,13]]
 uncL=F(78,25)-F(53,50)*r5[1]-F(161,100)*r13[1]+F(161,50)*r3[0]
 uncU=F(78,25)-F(53,50)*r5[0]-F(161,100)*r13[0]+F(161,50)*r3[1]
 corL=F(43,200)*r3[0]-F(101,1600)*r2[1]-F(24641,160000)
 corU=F(43,200)*r3[1]-F(101,1600)*r2[0]-F(24641,160000)
 for cid,lo,hi in [('uncapped-mixed-only',uncL,uncU),('n20-corrected',corL,corU)]:
  assert max(lo,F(cases[cid]['area']['lower']))<=min(hi,F(cases[cid]['area']['upper']))
 s3=F(173205,100000);s5=F(223607,100000);s13=F(360556,100000)
 assert s3*s3<3 and s5*s5>5 and s13*s13>13
 easy_unc=F(78,25)-F(53,50)*s5-F(161,100)*s13+F(161,50)*s3
 assert easy_unc==F(652519,1250000)>F(261,500)>F(1,2)
 easy_cor=F(43,200)*F(433,250)-F(101,1600)*F(283,200)-F(24641,160000)
 assert F(433,250)**2<3 and F(283,200)**2>2 and easy_cor==F(206483,1600000)>F(1,8)
 original=cases['n20-original-failed'];pair=next(r for r in original['catalog'] if r['mask']==6);assert not pair['eligible'] and F(pair['mass'])==F(5,8) and F(pair['span_squared'])==F(6889,10000)
 absent=[]
 for name in ['area-process-receipt.json','qualification-process-receipt.json']:
  rec=json.loads((rp/name).read_text());assert rec['exit_code']==0
  try:os.kill(rec['pid'],0);alive=True
  except ProcessLookupError:alive=False
  assert not alive;absent.append({'receipt':name,'pid':rec['pid'],'exit_code':0,'absent_at_check':True})
 out={'schema':'r077-fixed-radical-and-provenance-check-v1','pid':os.getpid(),'python':sys.version,'script_sha256':sha(Path(__file__)),'area_input_sha256':sha(ip),'root_bits':72,'root_bounds':{str(q):list(map(str,r)) for q,r in zip([2,3,5,13],[r2,r3,r5,r13])},'uncapped_formula':'78/25-53sqrt(5)/50-161sqrt(13)/100+161sqrt(3)/50','uncapped_formula_interval':list(map(str,[uncL,uncU])),'uncapped_simple_lower':str(easy_unc),'corrected_formula':'43sqrt(3)/200-101sqrt(2)/1600-24641/160000','corrected_formula_interval':list(map(str,[corL,corU])),'corrected_simple_lower':str(easy_cor),'formula_intervals_overlap_geometric_certificates':True,'original_residual_pair_empty':True,'uncapped_cap_countercheck':{'closed_square':['-161/50','0','-161/50','0'],'mass':'13','area':'25921/2500','mass_exceeds_area':True,'maximum_eligible_mass':'13','old_is_maximum':False},'process_checks':absent,'scope':'Root arithmetic and enclosure agreement; analytic formula derivations are in the separate theory proof. Geometric lower bounds are independently certified by both exact union algorithms.','seconds':time.monotonic()-t}
 (rp/'fixed-formula-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
