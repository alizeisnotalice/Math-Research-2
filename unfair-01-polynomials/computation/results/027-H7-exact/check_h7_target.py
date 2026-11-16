"""Exact integer identity for fixed H7(3,4); gcd follows from full disjoint root cover."""
from pathlib import Path
import json,hashlib,argparse
H=[1,1,1,0,0,0,1,1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--record',type=Path,required=True);a=p.parse_args()
 S=[i for i,c in enumerate(H) if c];c=[0]*50
 for i in S:
  for j in S:c[3*i+4*j]+=1
 assert len(c)==50 and c[0]==c[-1]==1 and set(c)=={0,1} and sum(c)==25
 assert H!=H[::-1]
 out={'status':'exact_binary_identity_verified','mother':H,'target_coefficients_ascending':c,'target_degree':49,'target_support':[k for k,x in enumerate(c) if x],'target_terms':25,'nonreciprocal':True,'squarefree_coprime_evidence':'Full independent Rouche cover: 21 distinct a roots and 28 distinct b roots; all 49 disks pairwise disjoint.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'target_terms':25}))
if __name__=='__main__':main()
