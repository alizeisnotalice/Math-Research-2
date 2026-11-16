import json
import argparse
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
T={0,1,2,3,6};S={i+8*j for i in T for j in T};assert len(S)==25
N=max(S);A={i:0 for i in range(N+1)if i not in S};B=dict(A)
A.update({0:1,1:0,2:1,3:0,6:0,7:0,8:0});B.update({0:1,1:1,2:0,3:0,6:1,7:0,8:0})
trace=[]
def force_zero(k,i,j,side):
 assert k not in S and i+j==k
 if side=='A':assert B[j]==1;A[i]=0
 else:assert A[i]==1;B[j]=0
 trace.append({'zero_target':k,'pair':[i,j],'forced':side+str(i if side=='A'else j)})
def force_one(k,side):
 pairs=[(i,k-i)for i in range(k+1)if A.get(i)!=0 and B.get(k-i)!=0]
 expected=[(k,0)]if side=='A'else[(0,k)]
 assert k in S and pairs==expected,(k,pairs)
 (A if side=='A'else B)[k]=1;trace.append({'unit_target':k,'only_pair':pairs[0],'forced':side+str(k)})
force_zero(15,9,6,'A');force_zero(12,2,10,'B');force_one(10,'A');force_zero(15,14,1,'A');force_one(14,'B')
assert 16 in S and A[10]*B[6]+A[2]*B[14]==2
out={'status':'exact_integer_pair_trace_verified','support':sorted(S),'trace':trace,'contradiction':{'degree':16,'target':1,'known_nonnegative_contributions':[[10,6,1],[2,14,1]],'lower_bound':2}}
args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
