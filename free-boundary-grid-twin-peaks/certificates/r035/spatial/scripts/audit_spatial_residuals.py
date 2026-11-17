#!/usr/bin/env python3
"""Exact residuals from certified first-five spatial layers; no float Q."""
from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from pathlib import Path
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
 spatial=json.loads((root/'evidence/n7-spatial-crt-audit.json').read_text()); full=json.loads((root/'references/r018-n7-distance-histogram.json').read_text())
 assert spatial['status']=='passed' and full['n']==7 and int(full['A2'])==spatial['global']['A2']
 prefix_E=sum(t*spatial['global']['E'][f'E{t}'] for t in range(1,6)); total_E=int(full['squared_euclidean_total']); residual=total_E-prefix_E;assert residual>=0
 A2=int(full['A2']); counts={int(k):int(v) for k,v in full['distance_counts'].items()};assert sum(counts.values())==A2
 fibre={}
 for dx in range(7):
  for dy in range(7): fibre.setdefault(dx*dx+dy*dy,set()).add(tuple(2*(dx+dy)>q for q in (7,14,21)))
 assert all(len(v)==1 for v in fibre.values())
 far={}
 for k,q in enumerate((7,14,21),1):
  total=sum(c for d,c in counts.items() if 2*d>q);prefix=sum(t*spatial['global']['far'][f'far{k}'][f'h{t}'] for t in range(1,6));rem=total-prefix;assert rem>=0;z=Fraction(rem,A2)
  far[f'far{k}']={'strict_threshold':f'2*d1>{q}','full_numerator':total,'first_five_numerator':prefix,'L_at_least_6_residual_numerator':rem,'L_at_least_6_residual_over_A2':{'numerator':z.numerator,'denominator':z.denominator}}
 result={'status':'passed','scope':'exact n7 Euclidean and Manhattan residuals after L=1..5; no float Q','inputs':{'spatial_audit_sha256':sha(root/'evidence/n7-spatial-crt-audit.json'),'r018_distance_histogram_sha256':sha(root/'references/r018-n7-distance-histogram.json')},'Euclidean_squared':{'full_numerator':total_E,'first_five_numerator':prefix_E,'L_at_least_6_residual_numerator':residual,'full_over_A2':{'numerator':Fraction(total_E,A2).numerator,'denominator':Fraction(total_E,A2).denominator},'residual_over_A2':{'numerator':Fraction(residual,A2).numerator,'denominator':Fraction(residual,A2).denominator}},'Manhattan_far_thresholds':far,'finite_fibre_check':{'r018_histogram_is_Manhattan':True,'squared_Euclidean_fibres_constant_for_events':True,'events':['2*d1>7','2*d1>14','2*d1>21']}}
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result['Euclidean_squared'],indent=2))
if __name__=='__main__':main()
