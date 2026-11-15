#!/usr/bin/env python3
"""R57 complete rectangle union area dual sweep on a prescribed one-sided order.

All intervals use outward dyadic integer arithmetic. Exact side orders extend
to the wall by strict interval comparisons for every other consecutive pair,
and the unique wall pair is handled by its exact equation and monotonicity.
"""
import argparse,gzip,hashlib,json,os,platform,sys,time
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
sys.dont_write_bytecode=True

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x

def digest(x):return hashlib.sha256(json.dumps(serial(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def ceilq(q):return -((-q.numerator)//q.denominator)
def fixed(q,D):return q.numerator*D//q.denominator,ceilq(q*D)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
    p=[a*b for a in x for b in y];return min(p),max(p)
def scale(x,c):return (x[0]*c,x[1]*c) if c>=0 else (x[1]*c,x[0]*c)
@lru_cache(None)
def sqrtq(q,bits):
    if q<0:raise ValueError('negative radicand')
    D=1<<bits;n=isqrt(q.numerator*D*D//q.denominator)
    assert n*n*q.denominator<=q.numerator*D*D<(n+1)*(n+1)*q.denominator
    return F(n,D),F(n,D) if F(n,D)**2==q else F(n+1,D)
def expr(e):
    a,v,s,q,b=e;a,v,q,b=map(F,(a,v,q,b));s=int(s)
    if not s:return a,v,0,F(0),F(0)
    if not b:
        n,d=isqrt(q.numerator),isqrt(q.denominator)
        if n*n==q.numerator and d*d==q.denominator:return a+s*F(n,d),v,0,F(0),F(0)
    return a,v,s,q,b

def boundq(e,dtlo,dthi,bits):
    a,v,s,q,b=e
    if not s:return (a,a),(F(0),F(0))
    qlo,qhi=sorted((q+b*dtlo,q+b*dthi))
    if qlo<=0:raise ValueError('radicand not uniformly positive')
    rl=sqrtq(qlo,bits)[0];rh=sqrtq(qhi,bits)[1]
    if rl<=0:raise ValueError('positive square root unresolved')
    val=(a+rl,a+rh) if s>0 else (a-rh,a-rl)
    c=s*b/2;der=(c/rh,c/rl) if c>=0 else (c/rl,c/rh)
    return val,der

def ordered(es,dt,bits):
    rows=[]
    for e in es:
        val,_=boundq(e,dt,dt,bits);rows.append((*val,e))
    rows.sort()
    bad=[(a[2],b[2]) for a,b in zip(rows,rows[1:]) if a[1]>=b[0]]
    if bad:raise ValueError('side order unresolved; distinct IDs or actual tie: '+str(serial(bad[0])))
    return [r[2] for r in rows]

class DualUnionTree:
    def __init__(self,y,dy):
        self.y,self.dy=y,dy;self.n=len(y)-1;N=4*self.n+4
        self.cover=[0]*N;self.value=[(0,0)]*N;self.derivative=[(0,0)]*N
    def update(self,a,b,delta):
        if a>=b:return
        def rec(node,lo,hi):
            if a<=lo and hi<=b:self.cover[node]+=delta
            else:
                mid=(lo+hi)//2
                if a<mid:rec(2*node,lo,mid)
                if b>mid:rec(2*node+1,mid,hi)
            if self.cover[node]<0:raise AssertionError('negative cover')
            if self.cover[node]:
                self.value[node]=sub(self.y[hi],self.y[lo]);self.derivative[node]=sub(self.dy[hi],self.dy[lo])
            elif hi-lo==1:self.value[node]=self.derivative[node]=(0,0)
            else:
                self.value[node]=add(self.value[2*node],self.value[2*node+1]);self.derivative[node]=add(self.derivative[2*node],self.derivative[2*node+1])
        rec(1,0,self.n)
    def result(self):return self.value[1],self.derivative[1]

def wall_interval(slab,bits):
    lookup={tuple(r['bbox']):box for r,box in zip(slab['box_records'],slab['boxes'])}
    a=expr(lookup[86,88,29,29][0]);b=expr(lookup[85,88,29,29][1]);c=a[0]-b[0]
    assert a[2]==-1 and b[2]==1 and a[4]==0 and b[4]!=0
    rl,rh=sqrtq(a[3],bits+80);u,v=c-rh,c-rl;assert u>0
    dt=sorted(((u*u-b[3])/b[4],(v*v-b[3])/b[4]));t0=F(slab['t0'])
    return [t0+d for d in dt],{'equation':'sqrt(qA)+sqrt(qBm+qBprime*(t-t0))=C','qA':a[3],'qBm':b[3],'qBprime':b[4],'C':c,'t0':t0,'A_endpoint':a,'B_endpoint':b}

def certify(slab,side,side_t=None,bits=192):
    start=time.perf_counter();D=1<<bits;t0=F(slab['t0']);wall,equation=wall_interval(slab,bits)
    if side_t is None:side_t=sum(wall)/2+(-1 if side=='left' else 1)*F(1,10**12)
    else:side_t=F(side_t)
    assert (side_t<wall[0]) if side=='left' else (side_t>wall[1])
    es=[tuple(expr(e) for e in box) for box in slab['boxes']];xl,xr,yl,yr=map(F,slab['domain']);fixedexpr=lambda q:(q,F(0),0,F(0),F(0));X0,X1,Y0,Y1=map(fixedexpr,(xl,xr,yl,yr))
    X={X0,X1};Y={Y0,Y1}
    for l,r,a,b in es:X.update((l,r));Y.update((a,b));assert a[1]==b[1]==0
    x=ordered(X,side_t-t0,bits);y=ordered(Y,side_t-t0,bits);xrank={e:i for i,e in enumerate(x)};yrank={e:i for i,e in enumerate(y)}
    seglo=min(side_t,wall[0])-t0;seghi=max(side_t,wall[1])-t0;ambiguous={}
    for axis,order in [('x',x),('y',y)]:
        intervals=[boundq(e,seglo,seghi,bits)[0] for e in order];bad=[]
        for i,(a,b) in enumerate(zip(intervals,intervals[1:])):
            if a[1]>=b[0]:bad.append({'rank':i,'left':order[i],'right':order[i+1],'gap_enclosure':sub(b,a)})
        ambiguous[axis]=bad
    # Every consecutive relation other than this exact crossing is strict on
    # the entire closed side-to-wall segment. The fixed A and decreasing B
    # meet only at the explicitly isolated root, proving both limiting orders.
    A,B=equation['A_endpoint'],equation['B_endpoint']
    if ambiguous['y'] or len(ambiguous['x'])!=1:raise ValueError('extra possible endpoint event; independent census required')
    pair=ambiguous['x'][0]
    if {tuple(pair['left']),tuple(pair['right'])}!={A,B}:raise ValueError('ambiguous pair is not the isolated exact wall')
    bder=boundq(B,seglo,seghi,bits)[1]
    assert A[4]==0 and bder[1]<0
    assert (tuple(pair['left']),tuple(pair['right']))==((A,B) if side=='left' else (B,A))
    continuation={'status':'passed','all_nonwall_consecutive_relations_strict_on_closed_side_wall_segment':True,
      'only_ambiguous_pair_is_exact_A_B_wall':True,'A_t_derivative':'0','B_t_derivative_bounds':bder,
      'argument':'A is constant and B strictly decreasing; isolated equation gives the unique equality t*. All other adjacent order relations stay strict throughout the closed segment.',
      'persistent_ties_with_distinct_full_endpoint_id_at_side':False,
      'parameter_side_interval':[min(side_t,wall[0]),max(side_t,wall[1])]}
    vals={};ders={}
    for e in X|Y:
        val,der=boundq(e,wall[0]-t0,wall[1]-t0,bits)
        vals[e]=(fixed(val[0],D)[0],fixed(val[1],D)[1]);ders[e]=(fixed(der[0],D)[0],fixed(der[1],D)[1])
    events=[[] for _ in x];ybounds=[]
    for i,(l,r,a,b) in enumerate(es):
        left,right=max(xrank[l],xrank[X0]),min(xrank[r],xrank[X1]);bottom,top=max(yrank[a],yrank[Y0]),min(yrank[b],yrank[Y1])
        if left>=right or bottom>=top:raise ValueError('empty candidate in side order')
        ybounds.append((bottom,top));events[left].append((1,i));events[right].append((-1,i))
    tree=DualUnionTree([vals[e] for e in y],[ders[e] for e in y]);area=at=(0,0);beta=(0,0);trace=[];active=0;maxactive=0
    for j,e in enumerate(x):
        for sign,i in events[j]:tree.update(*ybounds[i],sign);active+=sign
        maxactive=max(maxactive,active)
        if j+1==len(x):break
        width=sub(vals[x[j+1]],vals[e]);widtht=sub(ders[x[j+1]],ders[e]);height,heightt=tree.result()
        area=add(area,mul(width,height));at=add(at,add(mul(widtht,height),mul(width,heightt)))
        dv=x[j+1][1]-e[1]
        assert dv.denominator==1
        beta=add(beta,scale(height,int(dv)))
        if dv or events[j]:trace.append([j,len(events[j]),active])
    assert active==0 and tree.result()==((0,0),(0,0))
    out={'status':'passed_certified_one_sided_limit','scope':'complete union A,A_t,beta one-sided limit; exact endpoint continuation certified; conditional on supplied whole-mass candidate completeness and positivity','side':side,'side_t':side_t,'side_to_wall_order_certificate':continuation,'wall_t_bounds':wall,'wall_equation':equation,'bits':bits,'box_count':len(es),'x_endpoint_count':len(x),'y_endpoint_count':len(y),'ambiguous_consecutive_segment_pairs':ambiguous,'x_order_sha256':digest(x),'y_order_sha256':digest(y),'event_trace_sha256':digest(trace),'max_active':maxactive,'A_wall_bounds':[F(v,D*D) for v in area],'A_t_wall_bounds':[F(v,D*D) for v in at],'beta_wall_bounds':[F(v,D) for v in beta],'formula':'A_t=sum (xhi_t-xlo_t)*Y + (xhi-xlo)*Y_t; beta=sum(vxhi-vxlo)*Y','arithmetic':'outward exact dyadic integer values and t derivatives, with products at doubled precision','all_radicals_positive_on_side_wall_segment':True,'cpu_wall_seconds':time.perf_counter()-start,'python':sys.version,'platform':platform.platform()}
    return out

def selftest():
    D=1<<80;y=[(0,0),(D,D),(3*D,3*D)];dy=[(0,0),(D,D),(2*D,2*D)]
    t=DualUnionTree(y,dy);t.update(0,1,1);assert t.result()==((D,D),(D,D));t.update(1,2,1);assert t.result()==((3*D,3*D),(2*D,2*D));t.update(0,1,-1);assert t.result()==((2*D,2*D),(D,D));t.update(1,2,-1);assert t.result()==((0,0),(0,0))
    # Moving width and height both contribute: width=2, width_t=3,height=5,height_t=7.
    z=add(mul((2,2),(7,7)),mul((3,3),(5,5)));assert z==(29,29)
    for q in [F(1,9),F(2),F(7,10**30)]:
        a,b=sqrtq(q,80);assert a*a<=q<=b*b
    return {'status':'passed','tests':['dual coverage tree exact insertion/deletion','both product-rule terms contribute','exact dyadic radical inclusions']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--slab');p.add_argument('--side',choices=['left','right']);p.add_argument('--side-t');p.add_argument('--bits',type=int,default=192);p.add_argument('--output',required=True);p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:out=selftest()
    else:
        raw=Path(a.slab).read_bytes();payload=gzip.decompress(raw) if a.slab.endswith('.gz') else raw;out=certify(json.loads(payload),a.side,a.side_t,a.bits);out['slab_file_sha256']=hashlib.sha256(raw).hexdigest();out['slab_payload_sha256']=hashlib.sha256(payload).hexdigest()
    out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out['actual_command']=[sys.executable,*sys.argv]
    assert os.path.ismount('/Volumes/data');Path(a.output).write_text(json.dumps(serial(out),indent=2)+'\n');print(json.dumps(serial({k:v for k,v in out.items() if k not in ('ambiguous_consecutive_segment_pairs',)}),indent=2))
if __name__=='__main__':main()
