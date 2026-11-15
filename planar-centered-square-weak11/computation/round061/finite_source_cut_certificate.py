#!/usr/bin/env python3
"""Exact all-scale finite-source cut certificates; integer root enclosures."""
import argparse, hashlib, json, math, os, platform, random, sys, time
from collections import defaultdict,deque
from fractions import Fraction as F
from pathlib import Path
B=F(1994727599733432151936139,604462909807314587353088)

def rectangles(points,mass,bits):
    n=len(points);D=1<<bits;seen=set();inner=[];outer=[]
    for mask in range(1,1<<n):
        ids=[i for i in range(n) if mask>>i&1]
        a=min(points[i][0] for i in ids);b=max(points[i][0] for i in ids)
        c=min(points[i][1] for i in ids);d=max(points[i][1] for i in ids)
        closure=sum(1<<i for i,(x,y) in enumerate(points) if a<=x<=b and c<=y<=d)
        if closure in seen:continue
        seen.add(closure);m=sum(mass[i] for i in range(n) if closure>>i&1)
        if max(b-a,d-c)**2>=m:continue
        k=math.isqrt(m*(D//8)**2)
        assert k*k<=m*(D//8)**2<(k+1)**2
        for dest,r in ((inner,max(0,k-1)),(outer,k+1)):
            rect=(b*(D//4)-r,a*(D//4)+r,d*(D//4)-r,c*(D//4)+r,closure)
            if rect[0]<rect[1] and rect[2]<rect[3]:dest.append(rect)
    return inner,outer,len(seen)

def cells(inner,outer,n):
    ys=sorted({v for R in inner+outer for v in R[2:4]});yi={y:i for i,y in enumerate(ys)}
    events=defaultdict(list)
    for kind,rs in ((0,outer),(1,inner)):
        for x0,x1,y0,y1,mask in rs:
            for x,delta in ((x0,1),(x1,-1)):events[x].append((kind,yi[y0],yi[y1],mask,delta))
    xs=sorted(events);ny=len(ys)-1
    oc=[0]*ny; counts=[[0]*n for _ in range(ny)]; masks=[0]*ny;weights=[0]*(1<<n)
    dy=[ys[i+1]-ys[i] for i in range(ny)]
    for ix,x in enumerate(xs[:-1]):
        for kind,y0,y1,mask,delta in events[x]:
            if kind==0:
                for k in range(y0,y1):oc[k]+=delta
            else:
                ids=[i for i in range(n) if mask>>i&1]
                for k in range(y0,y1):
                    for i in ids:
                        counts[k][i]+=delta
                        assert counts[k][i]>=0
                        if counts[k][i]:masks[k]|=1<<i
                        else:masks[k]&=~(1<<i)
        dx=xs[ix+1]-x
        for k in range(ny):
            if oc[k]>0:weights[masks[k]]+=dx*dy[k]
    return weights,{'x_intervals':len(xs)-1,'y_intervals':ny,'nonzero_masks':sum(v>0 for v in weights)}

def cuts(weights,points,mass,bits):
    n=len(points);N=1<<n;z=weights[:];ms=[0]*N
    for i in range(n):
        for s in range(N):
            if s>>i&1:z[s]+=z[s^(1<<i)]
    raw=F(0);geom=F(0);rawmask=geommask=0;noncompact=0;compact=0;geom_paid=0
    D2=1<<(2*bits)
    for s in range(1,N):
        bit=s&-s;i=bit.bit_length()-1;ms[s]=ms[s^bit]+mass[i]
        ids=[j for j in range(n) if s>>j&1]
        span=max(max(points[j][0] for j in ids)-min(points[j][0] for j in ids),max(points[j][1] for j in ids)-min(points[j][1] for j in ids))
        ratio=F(16*z[s],D2*ms[s]);terminal=B*max(F(1),F(span*span,ms[s]));u=min(ratio,terminal)
        if span*span<=ms[s]:compact+=1
        else:noncompact+=1
        if terminal<ratio:geom_paid+=1
        if ratio>raw:raw=ratio;rawmask=s
        if u>geom:geom=u;geommask=s
    assert all(F(16*z[s],D2*ms[s])<=raw for s in range(1,N))
    return {'raw_bound':str(raw),'geometric_bound':str(geom),'raw_critical_mask':rawmask,'geometric_critical_mask':geommask,'compact_cuts':compact,'noncompact_cuts':noncompact,'terminal_improved_cuts':geom_paid,'all_cuts':N-1,'certifies_B':geom<=B},z,ms

def exact_flow(weights,mass,z,ms,critical):
    """Integer Dinic certificate for Hminus, at raw full-halo cut bound."""
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
    return {'integer_scale':scale,'area_denominator':'2^(2*bits)','source_capacity_multiplier':str(cap),'arcs':records,'source_load_numerators':[str(v) for v in loads],'total_routed_numerator':str(value),'scope':'Exact nonproportional flow on Hminus; full-H existence follows all actual-source cut certificate, not outward edges.'}

def generate(batch,case):
    rng=random.Random(610000+1000*batch+case);n=8+(batch+case)%5
    if batch%4==0:
        pts=rng.sample([(x,y) for x in range(-4,5,2) for y in range(-4,5,2)],n)
    elif batch%4==1:
        pts=[(3*i+((case+i)%2),0) for i in range(n)]
    elif batch%4==2:
        pts=rng.sample([(x+shift,y) for shift in (-16,16) for x in (-2,0,2) for y in (-2,0,2)],n)
    else:
        pts=rng.sample([(x,y) for x in (-8,-4,0,4,8) for y in (-8,-4,0,4,8)],n)
    mass=[rng.choice([1,4,9,16,25,36,64,100]) for _ in pts]
    if case%16==0:
        pts=[(-2,0),(2,0)]+[(32*i,16*(i%2)) for i in range(1,n-1)];mass=[16]*n
    return {'coordinate_unit':'1/4','mass_unit':'1/16','points_quarters':pts,'mass_sixteenths':mass,'seed':610000+1000*batch+case}

def run_case(inp,bits):
    points=inp['points_quarters'];mass=inp['mass_sixteenths'];t=time.perf_counter()
    inn,out,closed=rectangles(points,mass,bits);w,stats=cells(inn,out,len(points));res,z,ms=cuts(w,points,mass,bits)
    res.update({'input':inp,'bits':bits,'closed_source_masks':closed,'inner_rectangles':len(inn),'outer_rectangles':len(out),'arrangement':stats,'area_weight_numerators':[[s,str(v)] for s,v in enumerate(w) if v],'outer_area_numerator':str(sum(w)),'uncertain_fringe_numerator':str(w[0]),'flow':exact_flow(w,mass,z,ms,res['raw_critical_mask']),'wall_seconds':time.perf_counter()-t})
    return res

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,required=True);ap.add_argument('--cases',type=int,default=12);ap.add_argument('--bits',type=int,default=24);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--input',type=Path);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data');assert args.bits>=4
    t=time.perf_counter();inputs=json.loads(args.input.read_text())['batches'][args.batch]['inputs'] if args.input else [generate(args.batch,c) for c in range(args.cases)];results=[run_case(inp,args.bits) for inp in inputs]
    payload={'schema':'r061-finite-source-cuts-v1','status':'PASS','batch':args.batch,'python':sys.version,'platform':platform.platform(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'B':str(B),'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None,'cases':results,'max_bound':str(max(F(r['geometric_bound']) for r in results)),'all_certify_B':all(r['certifies_B'] for r in results),'total_cuts':sum(r['all_cuts'] for r in results),'wall_seconds':time.perf_counter()-t,'arithmetic':'integer/Fraction','randomness':'deterministic input seeds saved per case; no randomized computation','tolerance':'none; 2^-bits only enclosures, all tests exact'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps({k:payload[k] for k in ('status','batch','max_bound','all_certify_B','total_cuts','wall_seconds')}))
if __name__=='__main__':main()
