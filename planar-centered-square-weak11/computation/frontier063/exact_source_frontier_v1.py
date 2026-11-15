#!/usr/bin/env python3
"""Finite structural frontier. Integer geometry/cuts; no uniform-C2 claim.

Positions use quarters, masses sixteenths, threshold one. R61's strict
inner/outer construction and exact Hminus flow are preserved. This is a NEW
implementation; running it is required before claiming any case certified.
"""
import argparse
from array import array
from collections import defaultdict, deque
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import resource
import signal
import sys
import time

B = F(1994727599733432151936139, 604462909807314587353088)


def intervals(points, axis):
    """Masks of all closed coordinate intervals, without subset enumeration."""
    groups = defaultdict(int)
    for i, p in enumerate(points):
        groups[p[axis]] |= 1 << i
    values = sorted(groups)
    result = []
    for lo in range(len(values)):
        mask = 0
        for hi in range(lo, len(values)):
            mask |= groups[values[hi]]
            result.append(mask)
    return result


def rectangles(points, mass, bits):
    """Every tight bbox closure occurs as an x/y interval intersection.

    Conversely every nonempty intersection has a tight bbox whose closure
    is that intersection. Thus this replaces R61's 2^N geometric loop exactly.
    """
    closures = {x & y for x in intervals(points, 0)
                for y in intervals(points, 1)} - {0}
    D = 1 << bits
    inner, outer = [], []
    gamma = 0
    for mask in sorted(closures):
        ids = [i for i in range(len(points)) if mask >> i & 1]
        a = min(points[i][0] for i in ids)
        b = max(points[i][0] for i in ids)
        c = min(points[i][1] for i in ids)
        d = max(points[i][1] for i in ids)
        m = sum(mass[i] for i in ids)
        if max(b-a, d-c)**2 >= m:
            continue
        square = m * (D//8)**2
        k = math.isqrt(square)
        assert k*k <= square < (k+1)**2
        areas = []
        for target, r in ((inner, max(0, k-1)), (outer, k+1)):
            rect = (b*(D//4)-r, a*(D//4)+r,
                    d*(D//4)-r, c*(D//4)+r, mask)
            area = max(0, rect[1]-rect[0])*max(0, rect[3]-rect[2])
            areas.append(area)
            if area:
                target.append(rect)
        gamma += areas[1]-areas[0]
    return inner, outer, len(closures), gamma


def cells(inner, outer, n):
    """Exact outer-halo areas labeled by union of INNER eligible sources."""
    ys = sorted({v for r in inner+outer for v in r[2:4]})
    yi = {y:i for i,y in enumerate(ys)}
    events = defaultdict(list)
    ids = {r[4]:[i for i in range(n) if r[4] >> i & 1] for r in inner}
    for kind, rs in ((0, outer), (1, inner)):
        for x0,x1,y0,y1,mask in rs:
            events[x0].append((kind,yi[y0],yi[y1],mask,1))
            events[x1].append((kind,yi[y0],yi[y1],mask,-1))
    xs = sorted(events)
    ny = len(ys)-1
    oc = [0]*ny
    counts = [[0]*n for _ in range(ny)]
    masks = [0]*ny
    weights = [0]*(1 << n)
    dy = [ys[i+1]-ys[i] for i in range(ny)]
    for ix,x in enumerate(xs[:-1]):
        for kind,y0,y1,mask,delta in events[x]:
            if kind == 0:
                for k in range(y0,y1):
                    oc[k] += delta
                    assert oc[k] >= 0
            else:
                for k in range(y0,y1):
                    row = counts[k]
                    for i in ids[mask]:
                        row[i] += delta
                        assert row[i] >= 0
                        if row[i]:
                            masks[k] |= 1 << i
                        else:
                            masks[k] &= ~(1 << i)
        dx = xs[ix+1]-x
        for k in range(ny):
            if oc[k]:
                weights[masks[k]] += dx*dy[k]
    return weights, {'x_intervals':max(0,len(xs)-1),
                     'y_intervals':max(0,ny),
                     'nonzero_masks':sum(v > 0 for v in weights)}


def cuts(weights, points, mass, bits):
    """O(N 2^N), compact extrema/mass arrays, integer ratio comparisons."""
    n = len(points)
    N = 1 << n
    z = weights.copy()
    half = 1
    while half < N:
        block = half*2
        for base in range(0,N,block):
            for s in range(base+half,base+block):
                z[s] += z[s-half]
        half = block
    ms = array('q', [0])*N
    xmin = array('q', [0])*N
    xmax = array('q', [0])*N
    ymin = array('q', [0])*N
    ymax = array('q', [0])*N
    rawnum = geomnum = 0
    rawden = geomden = 1
    rawmask = geommask = 0
    compact = improved = 0
    D2 = 1 << (2*bits)
    Bd = B.denominator
    Bt = B.numerator*D2
    for s in range(1,N):
        bit = s & -s
        i = bit.bit_length()-1
        t = s ^ bit
        x,y = points[i]
        m = ms[t]+mass[i]
        ms[s] = m
        if t:
            xmin[s] = min(xmin[t],x)
            xmax[s] = max(xmax[t],x)
            ymin[s] = min(ymin[t],y)
            ymax[s] = max(ymax[t],y)
        else:
            xmin[s] = xmax[s] = x
            ymin[s] = ymax[s] = y
        span2 = max(xmax[s]-xmin[s],ymax[s]-ymin[s])**2
        compact += span2 <= m
        a = 16*z[s]
        if a*rawden > rawnum*m:
            rawnum,rawden,rawmask = a,m,s
        terminal = Bt*max(m,span2)
        numerical = a*Bd
        improved += terminal < numerical
        u = min(terminal,numerical)
        if u*geomden > geomnum*m:
            geomnum,geomden,geommask = u,m,s
    def critical(s):
        ids = [i for i in range(n) if s >> i & 1]
        span2 = max(max(points[i][a] for i in ids)-min(points[i][a] for i in ids)
                    for a in (0,1))**2
        return {'mask':s,'rank':s.bit_count(),'proper':s != N-1,
                'span_squared_over_mass':str(F(span2,ms[s])),
                'compact':span2 <= ms[s]}
    raw = F(rawnum,D2*rawden)
    geom = F(geomnum,Bd*D2*geomden)
    return {'raw_bound':str(raw),'geometric_bound':str(geom),
            'raw_critical_mask':rawmask,'geometric_critical_mask':geommask,
            'raw_critical':critical(rawmask),'geometric_critical':critical(geommask),
            'compact_cuts':compact,'noncompact_cuts':N-1-compact,
            'terminal_improved_cuts':improved,'all_cuts':N-1,
            'certifies_B':geom <= B}, z, ms


def exact_flow(weights,mass,z,ms,critical):
    """R61 integer Dinic flow on Hminus, unchanged mathematical capacities."""
    n=len(mass);act=[s for s in range(1,1<<n) if weights[s]];start=0;so=1;ro=1+n;sink=ro+len(act);g=[[] for _ in range(sink+1)]
    scale=ms[critical];cap=z[critical];total=sum(weights[s] for s in act)*scale;INF=total+1
    def edge(a,b,c):
        g[a].append([b,c,len(g[b]),c]);g[b].append([a,0,len(g[a])-1,0])
        return len(g[a])-1
    for i,m in enumerate(mass):edge(start,so+i,cap*m)
    arc=[]
    for j,s in enumerate(act):
        edge(ro+j,sink,weights[s]*scale)
        for i in range(n):
            if s>>i&1:arc.append((i,s,so+i,edge(so+i,ro+j,INF)))
    value=0
    while True:
        level=[-1]*len(g);level[start]=0;q=deque([start])
        while q:
            a=q.popleft()
            for b,c,_,_ in g[a]:
                if c and level[b]<0:level[b]=level[a]+1;q.append(b)
        if level[sink]<0:break
        it=[0]*len(g)
        def dfs(a,v):
            if a==sink:return v
            while it[a]<len(g[a]):
                e=g[a][it[a]];b,c,rev,_=e
                if c and level[b]==level[a]+1:
                    sent=dfs(b,min(v,c))
                    if sent:e[1]-=sent;g[b][rev][1]+=sent;return sent
                it[a]+=1
            return 0
        while True:
            v=dfs(start,INF)
            if not v:break
            value+=v
    assert value==total,(value,total)
    records=[];loads=[0]*n;dem=defaultdict(int)
    for i,s,a,k in arc:
        e=g[a][k];v=e[3]-e[1]
        if v:records.append([i,s,str(v)]);loads[i]+=v;dem[s]+=v
    assert all(dem[s]==weights[s]*scale for s in act)
    assert all(loads[i]<=cap*mass[i] for i in range(n))
    return {'integer_scale':scale,'area_denominator':'2^(2*bits)',
            'source_capacity_multiplier':str(cap),'arcs':records,
            'source_load_numerators':[str(v) for v in loads],
            'total_routed_numerator':str(value),
            'scope':'Exact flow on Hminus only; full-H existence uses all cuts.'}


def difference_area(left, right):
    """Integer area of union(left) minus union(right), ignoring mask fields."""
    ys = sorted({v for r in left+right for v in r[2:4]})
    yi = {y:i for i,y in enumerate(ys)}
    ev = defaultdict(list)
    for kind,rs in enumerate((left,right)):
        for x0,x1,y0,y1,_ in rs:
            ev[x0].append((kind,yi[y0],yi[y1],1))
            ev[x1].append((kind,yi[y0],yi[y1],-1))
    xs = sorted(ev)
    counts = [[0]*(len(ys)-1) for _ in range(2)]
    area = 0
    for j,x in enumerate(xs[:-1]):
        for kind,a,b,d in ev[x]:
            for k in range(a,b):
                counts[kind][k] += d
                assert counts[kind][k] >= 0
        height = sum(ys[k+1]-ys[k] for k in range(len(ys)-1)
                     if counts[0][k] and not counts[1][k])
        area += (xs[j+1]-x)*height
    return area


def deletions(points,mass,bits,inner,outer,maxrank,max_mass_compact=False):
    """All proper compact blocks of rank <= maxrank; keep old complement."""
    records=[];candidates={}
    # Every max-mass strictly eligible subset equals its bbox closure: any
    # omitted interior positive source would enlarge its mass while preserving
    # eligibility. The retained rectangle masks enumerate these closures.
    eligible={r[4] for r in outer}
    eligible_mass={s:sum(mass[i] for i in range(len(points)) if s>>i&1) for s in eligible}
    largest=max(eligible_mass.values())
    max_masks=sorted(s for s in eligible if eligible_mass[s]==largest)
    for rank in range(1,min(maxrank,len(points)-1)+1):
        for ids in itertools.combinations(range(len(points)),rank):
            span=max(max(points[i][a] for i in ids)-min(points[i][a] for i in ids)
                     for a in (0,1))
            m=sum(mass[i] for i in ids)
            if span*span>m:continue
            candidates.setdefault(sum(1<<i for i in ids),[]).append('bounded-rank-compact')
    if max_mass_compact:
        for mask in max_masks:
            if mask!=(1<<len(points))-1:
                candidates.setdefault(mask,[]).append('max-mass-strictly-eligible')
    for mask,reasons in sorted(candidates.items()):
            ids=[i for i in range(len(points)) if mask>>i&1];rank=len(ids)
            m=sum(mass[i] for i in ids)
            span=max(max(points[i][a] for i in ids)-min(points[i][a] for i in ids) for a in (0,1))
            keep=[i for i in range(len(points)) if not mask>>i&1]
            inn,out,closed,gamma=rectangles([points[i] for i in keep],[mass[i] for i in keep],bits)
            lo=difference_area(inner,out);hi=difference_area(outer,inn)
            assert 0<=lo<=hi
            denom=(1<<(2*bits))*m
            records.append({'mask':mask,'rank':rank,'selection_reasons':reasons,
                            'mass_sixteenths':m,'span_quarters':span,
                            'loss_numerator_lower':str(lo),'loss_numerator_upper':str(hi),
                            'remaining_inner_area_numerator':str(difference_area(inn,[])),
                            'remaining_outer_area_numerator':str(difference_area(out,[])),
                            'remaining_closed_source_masks':closed,
                            'remaining_inner_rectangles':len(inn),
                            'remaining_outer_rectangles':len(out),
                            'remaining_Gamma_numerator':str(gamma),
                            'ratio_lower':str(F(16*lo,denom)),
                            'ratio_upper':str(F(16*hi,denom))})
    return {'all_proper_compact_blocks_through_rank':maxrank,'blocks':records,
            'max_mass_compact_requested':max_mass_compact,
            'maximum_strictly_eligible_mass_sixteenths':largest,
            'all_maximum_strictly_eligible_masks':max_masks,
            'maximum_strictly_eligible_full_input':((1<<len(points))-1) in max_masks,
            'all_proper_maximum_ties_evaluated':max_mass_compact,
            'max_mass_min_upper':min((r['ratio_upper'] for r in records if 'max-mass-strictly-eligible' in r['selection_reasons']),key=F,default=None),
            'max_mass_min_lower':min((r['ratio_lower'] for r in records if 'max-mass-strictly-eligible' in r['selection_reasons']),key=F,default=None),
            'min_upper':min((r['ratio_upper'] for r in records),key=F,default=None),
            'min_lower':min((r['ratio_lower'] for r in records),key=F,default=None),
            'certified_cheap_block':any(F(r['ratio_upper'])<=F(7,2) for r in records),
            'all_tested_blocks_certified_expensive':bool(records) and all(F(r['ratio_lower'])>F(7,2) for r in records),
            'scope':'Input-specific full-halo deletion bounds; no universal deletion theorem.'}


def eligible_components(points,mass,rects):
    """Components of the true positive-rectangle source hypergraph."""
    components=[1<<i for i in range(len(points))]
    for *_,mask in rects:
        touched=[c for c in components if c & mask]
        joined=mask
        for c in touched:joined |= c
        components=[c for c in components if not c & mask]+[joined]
    result=[]
    for mask in sorted(set(components)):
        ids=[i for i in range(len(points)) if mask>>i&1]
        m=sum(mass[i] for i in ids)
        span=max(max(points[i][a] for i in ids)-min(points[i][a] for i in ids) for a in (0,1))
        result.append({'mask':mask,'rank':len(ids),'span_squared_over_mass':str(F(span*span,m)),
                       'terminal_35_33':33*span*span<=35*m})
    return result


def run_case(inp,bits,deletion_rank=0,mode='allcuts',max_mass_compact=False):
    points=inp['points_quarters'];mass=inp['mass_sixteenths'];n=len(points)
    assert 1<=n<=20 and n==len(mass) and len(set(map(tuple,points)))==n
    assert all(isinstance(m,int) and m>0 for m in mass) and sum(mass)<2**62
    assert all(len(p)==2 and all(isinstance(x,int) and abs(x)<2**30 for x in p) for p in points)
    started=time.perf_counter();times={}
    inner,outer,closed,gamma=rectangles(points,mass,bits)
    times['rectangles']=time.perf_counter()-started;t=time.perf_counter()
    result={}
    if mode=='allcuts':
        w,stats=cells(inner,outer,n);times['cells']=time.perf_counter()-t;t=time.perf_counter()
        result,z,ms=cuts(w,points,mass,bits);times['cuts']=time.perf_counter()-t;t=time.perf_counter()
        flow=exact_flow(w,mass,z,ms,result['raw_critical_mask'])
        times['flow']=time.perf_counter()-t
        result.update({'arrangement':stats,'area_weight_numerators':[[s,str(v)] for s,v in enumerate(w) if v],
                       'outer_area_numerator':str(sum(w)),'uncertain_fringe_numerator':str(w[0]),'flow':flow})
        del z,ms,w
    else:
        assert deletion_rank>0 or max_mass_compact
        result.update({'all_cuts':0,'full_inner_area_numerator':str(difference_area(inner,[])),
                       'full_outer_area_numerator':str(difference_area(outer,[]))})
    components=eligible_components(points,mass,outer)
    # All retained true boxes are nonempty at strict radius. Outer rects have
    # exactly those same source closures, so their masks give this exact graph.
    result.update({'input':inp,'bits':bits,'mode':mode,'closed_source_masks':closed,
                   'inner_rectangles':len(inner),'outer_rectangles':len(outer),
                   'Gamma_numerator':str(gamma),'Gamma_over_min_mass':str(F(16*gamma,(1<<(2*bits))*min(mass))),
                   'eligible_components':components,
                   'all_components_terminal_35_33':all(c['terminal_35_33'] for c in components)})
    if deletion_rank or max_mass_compact:
        t=time.perf_counter();result['deletions']=deletions(points,mass,bits,inner,outer,deletion_rank,max_mass_compact)
        times['deletions']=time.perf_counter()-t
    result['stage_seconds']=times;result['wall_seconds']=time.perf_counter()-started
    return result


def atomic_json(path,payload):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(payload,indent=2)+'\n');tmp.replace(path)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--stratum',type=int,default=0)
    ap.add_argument('--case',type=int)
    ap.add_argument('--limit',type=int)
    ap.add_argument('--bits',type=int,default=24)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--deletion-rank',type=int,choices=(0,1,2,3),default=0)
    ap.add_argument('--mode',choices=('allcuts','deletions'),default='allcuts')
    ap.add_argument('--max-mass-compact',action='store_true')
    ap.add_argument('--seconds-per-case',type=int,default=900)
    ap.add_argument('--memory-mib',type=int,default=768)
    args=ap.parse_args()
    assert os.path.ismount('/Volumes/data') and 4<=args.bits<=128
    assert args.seconds_per_case>0 and args.memory_mib>0
    assert args.mode!='deletions' or args.deletion_rank>0 or args.max_mass_compact
    resource.setrlimit(resource.RLIMIT_AS,(args.memory_mib*2**20,)*2)
    def timeout(signum,frame):raise TimeoutError('case wall budget exceeded')
    signal.signal(signal.SIGALRM,timeout)
    data=json.loads(args.input.read_text());stratum=data['strata'][args.stratum]
    cases=list(enumerate(stratum['inputs']))
    if args.case is not None:cases=[cases[args.case]]
    if args.limit is not None:cases=cases[:args.limit]
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result={'schema':'exact-source-frontier-v1','status':'RUNNING','stratum':args.stratum,
            'stratum_name':stratum['name'],'B':str(B),'script_sha256':source_hash,
            'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
            'python':sys.version,'platform':platform.platform(),
            'budget':{'seconds_per_case':args.seconds_per_case,'address_space_mib':args.memory_mib},
            'arithmetic':'Exact Python integers/Fraction; radii enclosed at 2^-bits.',
            'cases':[]}
    started=time.perf_counter()
    for case,inp in cases:
        try:
            signal.alarm(args.seconds_per_case)
            r=run_case(inp,args.bits,args.deletion_rank,args.mode,args.max_mass_compact);r['case_index']=case
            signal.alarm(0);result['cases'].append(r)
            result['wall_seconds']=time.perf_counter()-started
            result['maxrss_platform_units']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if args.mode=='allcuts':
                bound=F(r['geometric_bound'])
                result['status']='STOP_BUDGET_ABOVE_B' if not r['certifies_B'] else ('STOP_BOUND_ABOVE_7_4' if bound>F(7,4) else 'RUNNING')
            atomic_json(args.out,result)
            print(json.dumps({'case':case,'bound':r.get('geometric_bound'),
                              'minimum_deletion_upper':r.get('deletions',{}).get('min_upper'),
                              'seconds':r['wall_seconds'],'status':result['status']}),flush=True)
            if result['status']!='RUNNING':return 2
        except Exception as exc:
            signal.alarm(0);result['status']='FAILED';result['failed_case_index']=case
            result['error']=repr(exc);atomic_json(args.out,result);raise
    result['status']='PASS';result['wall_seconds']=time.perf_counter()-started
    result['total_cuts']=sum(r['all_cuts'] for r in result['cases'])
    atomic_json(args.out,result)
    return 0


if __name__=='__main__':sys.exit(main())
