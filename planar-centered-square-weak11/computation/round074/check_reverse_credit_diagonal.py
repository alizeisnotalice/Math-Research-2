#!/usr/bin/env python3
"""One deterministic exact diagonal check; no1200-case area replay."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,platform,random,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--probe-root',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();root=Path(a.probe_root);sp=root/'probe_reverse_credit.py';op=root/'reverse-credit-probe.json';rp=root/'reverse-credit-probe-receipt.json';d=json.loads(op.read_text());receipt=json.loads(rp.read_text());assert sha(sp)==d['script_sha256']==receipt['script_sha256'];assert sha(op)==receipt['output_sha256'];argv=receipt['argv'];gp=Path(argv[argv.index('--geometry')+1]);assert sha(gp)==d['geometry_sha256'];assert len(d['all_results'])==d['count']==1200
 top=d['top20'][0];assert top['count']==2 and top['seed']==742090 and top['mode']==1
 # Reconstruct ONLY this selected seed's positions/rawweights independently.
 rng=random.Random(top['seed']);points=set()
 while len(points)<2:
  x=rng.choice((-10,10));y=rng.randint(-10,10)
  if rng.randrange(2):x,y=y,x
  if (x,y)!=(0,0):points.add((x,y))
 points=sorted(points);weights=[rng.randint(1,20) for _ in points];groups={}
 for p,m in zip(points,weights):t=max(abs(p[0]),abs(p[1]));groups[t]=groups.get(t,0)+m
 cumulative=0;scale=None
 for t,m in sorted(groups.items()):
  cumulative+=m;b=F(4*t*t,400*cumulative);scale=b if scale is None else min(scale,b)
 masses=[scale*m for m in weights];positions=[(F(x,20),F(y,20)) for x,y in points];assert [list(map(str,p)) for p in positions]==top['positions'];assert list(map(str,masses))==top['masses'];assert positions==[(F(-1,2),F(-1,2)),(F(1,2),F(1,2))];assert masses==[F(8,15),F(7,15)]
 prefixes=[]
 for t in sorted(set(max(map(abs,p)) for p in positions)):
  mass=sum(m for p,m in zip(positions,masses) if max(map(abs,p))<=t);assert mass<=4*t*t;prefixes.append({'radius':str(t),'mass':str(mass),'capacity':str(4*t*t),'slack':str(4*t*t-mass)})
 # Complete nonempty source-subset catalog for exactly2points:1,2,3.
 subsets=[]
 for mask in [1,2,3]:
  pp=[p for i,p in enumerate(positions) if mask>>i&1];mass=sum(m for i,m in enumerate(masses) if mask>>i&1);span=max(max(p[j] for p in pp)-min(p[j] for p in pp) for j in (0,1));subsets.append({'mask':mask,'mass':str(mass),'span_squared':str(span*span),'eligible':mass>span*span})
 assert [x['eligible'] for x in subsets]==[True,True,False];assert top['closure_count']==3
 radius_squared=[m/4 for m in masses];assert all(0<r<F(1,4) for r in radius_squared)
 # Each singleton radius is<1/2. Its square clips to a quarter-square
 # insideQ; the opposite quarters are disjoint since both radii<1/2.
 clipped_areas=radius_squared;area=sum(clipped_areas);total=sum(masses);ratio=total/area;assert area==F(1,4) and total==1 and ratio==4
 out={'schema':'r074-exact-top-diagonal-reverse-credit-check-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'reviewed_probe_script_path':str(sp),'reviewed_probe_script_sha256':sha(sp),'reviewed_probe_output_path':str(op),'reviewed_probe_output_sha256':sha(op),'reviewed_probe_receipt_path':str(rp),'reviewed_probe_receipt_sha256':sha(rp),'reviewed_geometry_source_path':str(gp),'reviewed_geometry_source_sha256':sha(gp),'selected_seed':top['seed'],'raw_integer_positions':points,'raw_integer_weights':weights,'exact_scale':str(scale),'exact_positions':top['positions'],'exact_masses':top['masses'],'complete_center_prefixes':prefixes,'complete_nonempty_subsets':subsets,'singleton_radius_squared':list(map(str,radius_squared)),'clipped_singleton_areas':list(map(str,clipped_areas)),'exact_closed_Q_mass':str(total),'exact_clipped_halo_area':str(area),'exact_ratio':str(ratio),'root_area_float':top['area_float'],'root_ratio_float':top['mass_over_area_float'],'floating_ratio_minus_exact':top['mass_over_area_float']-float(ratio),'interpretation':'Exact equality4 for closedQ. The pair mask is ineligible at mass=spanSquared=1; only two singleton halo squares remain, each contributes its quarter area insideQ. Center prefix t=1/2 has equality, so center is outside strict halo.','scope':'Only one selected seed reconstructed and exactly checked. Metadata/source hashes read for1200-case probe; no other halo area replay. No general reverse-credit theorem, finite-class exhaustive result, Borel extension, or open-Q equality assertion.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'raw_weights':weights,'scale':str(scale),'exact_area':str(area),'exact_ratio':str(ratio),'root_float_ratio':top['mass_over_area_float']},indent=2))
if __name__=='__main__':main()
