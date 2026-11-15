#!/usr/bin/env python3
"""Independent small support checks: inclusion-exclusion area and Prim MST.

No author module import, Kruskal code, compressed-cell area or W evaluation.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(e,saved):
    p=[tuple(map(F,v))for v in e['points']];m=list(map(F,e['masses']));n=e['S_count'];S=(1<<n)-1
    eligible=[];parent=list(range(len(p)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,1<<len(p)):
        ids=[i for i in range(len(p))if mask>>i&1];mass=sum(m[i]for i in ids);span=max(max(p[i][j]for i in ids)-min(p[i][j]for i in ids)for j in (0,1))
        if span*span<mass:
            eligible.append((mask,mass))
            for i in ids[1:]:parent[find(i)]=find(ids[0])
    L=max(v for k,v in eligible);assert L==1 and [k for k,v in eligible if v==L]==[S]and len({find(i)for i in range(len(p))})==1
    # All normalized unit squares; inclusion-exclusion has only2^8-1 terms.
    area=F(0)
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1]
        dx=max(F(0),1-(max(p[i][0]for i in ids)-min(p[i][0]for i in ids)));dy=max(F(0),1-(max(p[i][1]for i in ids)-min(p[i][1]for i in ids)))
        area+=(1 if len(ids)%2 else-1)*dx*dy
    def cost(i,j):
        dx,dy=abs(p[i][0]-p[j][0]),abs(p[i][1]-p[j][1]);return dx+dy-dx*dy
    chosen={0};mst=F(0);edges=[]
    while len(chosen)<n:
        c,i,j=min((cost(i,j),i,j)for i in chosen for j in range(n)if j not in chosen);mst+=c;chosen.add(j);edges.append([i,j,str(c)])
    chain=False
    for sx in (-1,1):
        for sy in (-1,1):
            pts=sorted((sx*x,sy*y)for x,y in p[:n])
            if all(pts[i][1]<=pts[i+1][1]for i in range(n-1)):chain=True
    span=max(max(x for x,y in p)-min(x for x,y in p),max(y for x,y in p)-min(y for x,y in p));margin=span*span-F(35,33)*sum(m)
    assert margin>0 and area==F(saved['support_union_area_exact'])and mst==F(saved['MST_cost_exact'])and 1+mst==F(saved['tree_upper_exact'])and chain==saved['chain_after_axis_reflection']==e['expected_chain']
    assert 1+mst==F(e['expected_tree_bound'])
    return dict(id=e['id'],source_subsets=(1<<len(p))-1,unique_maximum_mask=S,nonterminal_margin=str(margin),support_union_area_by_inclusion_exclusion=str(area),MST_cost_by_Prim=str(mst),Prim_edges=edges,chain_after_reflection=chain,tree_pays=mst<=F(5,2),W_states_not_evaluated=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--expected',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not sys.flags.optimize and Path('/Volumes/data').is_mount()
    assert sha(a.input)=='25ccd96bb2e35a37c78ce680a12654a67d25db8adcb87a270d4c22b47cce7aa4'and sha(a.expected)=='489b48eaedb5d642fe56e8fe697f169361d4afebfcf173e24086a68d804bbfc8'
    out=a.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();d=json.loads(a.input.read_text());s=json.loads(a.expected.read_text());by={e['id']:e for e in s['examples']};rows=[one(e,by[e['id']])for e in d['examples']]
    result=dict(status='PASS_R68_INDEPENDENT_SUPPORT_TREE_FIXTURES',python=sys.version,pid=os.getpid(),entry_sha256=sha(Path(__file__)),input_sha256=sha(a.input),expected_sha256=sha(a.expected),examples=rows,elapsed_seconds=time.perf_counter()-start,scope='Three fixed source-domain fixtures; exact inclusion-exclusion support area and Prim MST, independent of author compressed-cell/Kruskal. No W/Wprime computed or universal support claim.')
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
