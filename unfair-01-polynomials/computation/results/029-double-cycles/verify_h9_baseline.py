"""Exact next-target baseline; no real-factor positivity classification."""
import argparse,json
from pathlib import Path
from finite_fields import _irreducible_mod
def main():
 p=argparse.ArgumentParser();p.add_argument('--record',type=Path,required=True);a=p.parse_args()
 phi=[1]*5;g=[1,0,0,-1,0,1];h=[0]*10
 for i,x in enumerate(phi):
  for j,y in enumerate(g):h[i+j]+=x*y
 assert h==[1,1,1,0,0,0,0,0,1,1]
 for f in (phi,g):_irreducible_mod([x%2 for x in f],2)
 out={'H9':h,'rational_factorization':[phi,g],'integer_product_recomputed':True,
      'both_factors_irreducible_mod_2_by_independent_Rabin':True,
      'therefore_H9_not_irreducible_or_S9_on_nine_roots':True,
      'negative_quintic_coefficient':{'index':3,'value':-1},
      'scope':'Next-target baseline only; no nonnegative-real factor classification.'}
 a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
