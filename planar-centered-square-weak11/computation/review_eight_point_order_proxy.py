#!/usr/bin/env python3
"""Independent exhaustive subset and two-dimensional cell audit of R40 proxy."""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,json,math,platform,sys,time
from itertools import combinations

def cell_union(rects):
    xs=sorted({q for r in rects for q in r[:2]})
    ys=sorted({q for r in rects for q in r[2:]})
    total=0
    for x0,x1 in zip(xs,xs[1:]):
      xm=x0+x1
      for y0,y1 in zip(ys,ys[1:]):
        ym=y0+y1
        if any(2*a<xm<2*b and 2*c<ym<2*d for a,b,c,d in rects):
          total+=(x1-x0)*(y1-y0)
    return total

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    start=time.perf_counter()
    result_path=args.source_dir/'exact-eight-point-proxy.json';ledger_path=args.source_dir/'all-order-signatures.csv'
    data=json.loads(result_path.read_text());rows=data['signatures']
    points=[(int(F(x)*100),int(F(y)*100)) for x,y in data['input']['points']]
    assert points==[(x,y) for x in [-49,0,49] for y in [-49,0,49] if (x,y)!=(0,0)]
    assert data['input']['masses']==['1/8']*8
    # Directly enumerate every submask, independent of the producer's subset DP.
    spans={}
    for mask in range(1,256):
      pp=[points[i] for i in range(8) if mask>>i&1]
      spans[mask]=max(max(x for x,y in pp)-min(x for x,y in pp),max(y for x,y in pp)-min(y for x,y in pp))**2
    choices={};subsets_checked=0
    for j in range(8):
      for old in range(256):
        if old>>j&1:continue
        best=1;sub=old
        while True:
          full=sub|(1<<j); count=full.bit_count();subsets_checked+=1
          if 8*spans[full]<10000*count:best=max(best,count)
          if sub==0:break
          sub=(sub-1)&old
        choices[j,old]=best
    seen={'corner':set(),'edge_midpoint':set()};multiplicity=[0]*len(rows)
    with ledger_path.open(newline='') as f:
      for record in csv.DictReader(f):
        name=record['last_orbit']; order=tuple(map(int,record['order_vertex_indices']));rid=int(record['signature_id'])
        assert sorted(order)==list(range(8)) and order[-1]==(0 if name=='corner' else 1)
        assert order not in seen[name];seen[name].add(order)
        mask=0;typ=[0]*8
        for j in order:typ[j]=choices[j,mask];mask|=1<<j
        assert typ==rows[rid]['types'];multiplicity[rid]+=1
    assert [len(seen[x]) for x in ['corner','edge_midpoint']]==[math.factorial(7)]*2
    assert sum(multiplicity)==10080 and len(rows)==348
    bits=data['bits'];D=100*(1<<bits);assert D*D==data['area_denominator']
    root_bounds={}
    for k in range(1,9):
      q=k*(1<<(2*bits-5));a=math.isqrt(q);b=a+(a*a<q)
      assert a*a<=q<=b*b
      root_bounds[k]=(100*a,100*b)
    all_min=[];all_max=[]
    for rid,row in enumerate(rows):
      assert row['id']==rid and row['multiplicity']==multiplicity[rid]
      pair=[]
      for end in [0,1]:
        rr=[]
        for (px,py),k in zip(points,row['types']):
          cx=px*(1<<bits);cy=py*(1<<bits);r=root_bounds[k][end]
          rr.append((cx-r,cx+r,cy-r,cy+r))
        pair.append(cell_union(rr))
      assert pair==[row['inner_area_numerator'],row['outer_area_numerator']]
      assert pair[0]>2*D*D
      all_min.append(pair[0]);all_max.append(pair[1])
    low=F(min(all_min),D*D);high=F(min(all_max),D*D)
    assert str(low)==data['minimum_proxy_area_lower'] and str(high)==data['minimum_proxy_area_upper']
    candidates=[i for i,v in enumerate(all_min) if v<=min(all_max)]
    assert candidates==[237]
    second_lower=F(min(v for i,v in enumerate(all_min) if i!=237),D*D)
    assert second_lower>high
    winner=[]
    for end in [0,1]:
      rr=[]
      for (px,py),k in zip(points,rows[237]['types']):
        cx=px*(1<<bits);cy=py*(1<<bits);r=root_bounds[k][end]
        rr.append((cx-r,cx+r,cy-r,cy+r))
      winner.append(rr)
    def intersection(rs):
      return max(0,min(r[1] for r in rs)-max(r[0] for r in rs))*max(0,min(r[3] for r in rs)-max(r[2] for r in rs))
    pairgraphs=[]
    for rr in winner:
      pairgraphs.append([list(ij) for ij in combinations(range(8),2) if intersection([rr[i] for i in ij])>0])
    assert pairgraphs[0]==pairgraphs[1]==[[0,1],[1,2],[1,3],[1,4],[3,6],[4,6]]
    assert all(intersection([winner[1][i] for i in ijk])==0 for ijk in combinations(range(8),3))
    assert hashlib.sha256(ledger_path.read_bytes()).hexdigest()==data['ledger_sha256']
    output=dict(status='PASS',python=sys.version,platform=platform.platform(),algorithm='Direct all-submask eligibility; complete distinct-order ledger cardinality; exact 2D endpoint-cell union area (independent of source sweep).',subsets_checked=subsets_checked,orders_checked=10080,signatures_checked=348,inner_outer_union_areas_checked=696,minimum_lower=str(low),minimum_upper=str(high),strict_excess_over_two=str(low-2),unique_optimal_signature=237,next_signature_lower=str(second_lower),winner_positive_overlap_pairs=pairgraphs[0],winner_positive_triple_overlaps=0,winner_exact_area='18121/10000 + 3*sqrt(6)/25 - sqrt(2)/200 (analytic pair-overlap expansion in review)',arithmetic='integers and fractions.Fraction; no random seed or tolerance',elapsed_seconds=time.perf_counter()-start,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [result_path,ledger_path,Path(__file__).resolve()]},scope='All permutations of this exact eight-point proxy via two D4 final-vertex orbits; not a C2 lower bound.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()
