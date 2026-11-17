#!/usr/bin/env python3
"""Exact module kernels versus tree budgets for every partition of 2,4,8 leaves."""
import argparse,hashlib,json,platform,sys,time
from pathlib import Path

def row_basis(vs):
    piv={}
    for v in vs:
        while v:
            i=v.bit_length()-1
            if i in piv:v^=piv[i]
            else:piv[i]=v;break
    return piv
def reduce_by(v,piv):
    for i in sorted(piv,reverse=True):
        if v>>i&1:v^=piv[i]
    return v
def partitions(n):
    a=[0]
    def rec(maximum):
        if len(a)==n:
            yield tuple(a);return
        for t in range(maximum+2):
            a.append(t);yield from rec(max(maximum,t));a.pop()
    yield from rec(0)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();start=time.monotonic();records=[]
    for ell in [1,2,3]:
        Q=1<<ell;N=Q*Q
        h=[sum(1<<(z+Q*w) for z in range(Q)) for w in range(Q)]
        tri=[];nodes=[]
        for j in range(ell):
            for prefix in range(1<<j):
                leafs=[w for w in range(Q) if w&((1<<j)-1)==prefix]
                tri.append(sum(1<<(z+Q*w) for z in range(Q) for w in leafs if z>>j&1))
                nodes.append(leafs)
        piv=row_basis(h+tri);assert len(piv)==2*Q-1
        summary={};count=0
        for part in partitions(Q):
            m=max(part)+1;cells=[sum(h[w] for w in range(Q) if part[w]==c) for c in range(m)]
            columns=[]
            for g in tri:
                columns.append(sum(reduce_by(g&mask,piv)<<(N*c) for c,mask in enumerate(cells)))
            kernel_tri=Q-1-len(row_basis(columns))
            homogeneous=sum(len({part[w] for w in leaves})==1 for leaves in nodes)
            assert kernel_tri==homogeneous and kernel_tri<=Q-m
            dimB=Q+kernel_tri;assert dimB+m<=2*Q
            key=str(m)
            if key not in summary:summary[key]={'partitions':0,'equality_count':0,'min_dimB':dimB,'max_dimB':dimB}
            s=summary[key];s['partitions']+=1;s['equality_count']+=dimB+m==2*Q
            s['min_dimB']=min(s['min_dimB'],dimB);s['max_dimB']=max(s['max_dimB'],dimB)
            count+=1
        assert count=={2:2,4:15,8:4140}[Q]
        records.append({'Q':Q,'full_code_dimension':2*Q-1,'partitions_checked':count,'by_cell_count':summary})
    out={'status':'passed','arithmetic':'F2 truth-mask kernel computation','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'host':platform.node(),'platform':platform.platform(),'python':sys.version,'records':records,'elapsed_seconds':time.monotonic()-start,
         'scope':'Single base-algebra atom; arbitrary dimensions and multiple atoms are covered by the separate proof.'}
    args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
