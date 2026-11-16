#!/usr/bin/env python3
"""Read saved finite screen only; no new halo, caps, pairs or mass points."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import argparse,gzip,hashlib,json,os,sys
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();p=Path(a.input);d=json.loads(gzip.decompress(p.read_bytes()));counts=Counter();top={};position_null=0
 for pair in d['position_pairs']:
  position_null+=pair['all_feasible_mass_points_have_null_formal_triple']
  for m in pair['mass_points']:
   arity=sum(F(b)>0 for b in m['masses']);active=m['formal_triple_Gamma_active'];key=f'arity{arity}_'+('active' if active else'null');counts[key]+=1;record={'pair_index':pair['pair_index'],'category':pair['category'],'positions_units':pair['positions_units'],'V1':pair['V1'],'V2':pair['V2'],'U':pair['U']}|m
   if key not in top or m['Delta_Gamma_float']>top[key]['Delta_Gamma_float']:top[key]=record
 out={'schema':'r076-saved-screen-grouping-v1','python':sys.version,'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'screen_data_path':str(p),'screen_data_sha256':sha(p),'actual_position_pairs':len(d['position_pairs']),'position_pairs_whose_entire_feasible_mass_polytope_has_null_formal_triple':position_null,'counts_by_positive_outside_arity_and_formal_triple':dict(counts),'float_top_by_group':top,'new_position_pairs':0,'new_mass_points':0,'new_halo_evaluations':0,'scope':'Grouping of existing finite screen only. Zero-mass candidate positions may make a formal triple label active, so positive-outside arity is separately reported. No general inequality or rigorous float area bound follows.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'position_null':position_null,'counts':dict(counts),'top_2active':top.get('arity2_active'),'top_2null':top.get('arity2_null')},indent=2))
if __name__=='__main__':main()
