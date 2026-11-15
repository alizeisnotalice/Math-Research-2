"""Read-only R39 source/coverage/active-box audit; does not replay 44M boxes."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from datetime import datetime,timezone
from itertools import product
from math import isqrt,lcm
import argparse,hashlib,json,sys,time
PIN='acede2c710655d0c179e7bea3ee9483271caaf37968227e7e9fce6442907056b'
SLABPIN='c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def ceil(q):return (q.numerator+q.denominator-1)//q.denominator

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--slab',type=Path,required=True);ap.add_argument('--velocities',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();d=read(a.slab);base=read(a.input);vel=read(a.velocities);assert sha(a.module)==d['code_sha256']==PIN and sha(a.slab)==SLABPIN
    assert sha(a.input)==d['input_file_sha256']=='ca9f8b0ea3d78036151ed3e5c37ffd558c395fc2de66c273f561f4ef86630f59'and digest(base)==d['input_canonical_sha256']=='7430136b07f641f39e4e26cd12fee6a6fbb5a5ca317a688f4848d90ea03666fa'
    assert sha(a.velocities)==d['direction_file_sha256']and digest(vel)==d['direction_sha256']and vel==d['velocities']and len(vel)==204 and [j for j,v in enumerate(vel)if v]==[130,131]and all(v in(0,1)for v in vel)
    px,py,lam=map(F,(base['Px'],base['Py'],base['lam']));M=sum(F(p[2])for p in base['pts']);s0,H0,H,dt=map(F,(d['s0'],d['requested_H'],d['H'],d['delta_empty']));R=F(d['tail']['radius']);fl,fr,fb,ft=map(F,d['domain']);assert s0==F(5588494330693,500000000000000000)and H0==H==dt==F(328734960629,500000000000000000)and R==10 and lam==1 and(fl,fr,fb,ft)==(0,px,0,py)
    assert min(F(p[2])for p in base['pts'])>0 and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in base['pts'])
    tail=M/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(d['tail']['bound'])<lam and M==F(d['tail']['mass'])and F(d['tail']['lambda'])==lam
    S=[]
    for label,((xx,yy,ww),v)in enumerate(zip(base['pts'],vel)):
        x,y,w=map(F,(xx,yy,ww));xc=x+v*s0
        imin=ceil((fl-R-xc-v*H0)/px);z=(fr+R-xc+v*H0)/px;imax=z.numerator//z.denominator;jmin=ceil((fb-R-y)/py);z=(ft+R-y)/py;jmax=z.numerator//z.denominator
        for i,j in product(range(imin,imax+1),range(jmin,jmax+1)):S.append({'label':label,'period_shift':[i,j],'x':str(xc+i*px),'y':str(y+j*py),'v':v,'w':str(w)})
    assert S==d['source_catalog']and len(S)==3953 and digest(S)==d['source_catalog_sha256']
    xc=sorted({(F(p['x']),p['v'])for p in S});ys=sorted({F(p['y'])for p in S});xs=[z[0]for z in xc];vs=[z[1]for z in xc];assert len(xs)==225 and len(ys)==62 and [[str(x),v]for x,v in xc]==d['x_classes']and list(map(str,ys))==d['y_classes']and digest([d['x_classes'],d['y_classes']])==d['source_order_sha256']
    assert all(x<y for x,y in zip(xs,xs[1:]));ec=Counter();radii=[]
    def event(gap,speed,kind):
        if not speed:ec[kind+('_persistent_equal'if not gap else'_parallel')]+=1;return
        assert gap;radius=abs(gap/speed)/2;assert radius>=H;ec[kind+'_nonpersistent']+=1;radii.append(radius)
    for i,(x,v)in enumerate(xc):
        for y,w in xc[i+1:]:event(y-x,w-v,'source_order');event(y-x-2*R,w-v,'span_2R')
        event(x-(fr+R),v,'reach_right');event(x-(fl-R),v,'reach_left')
    assert dict(ec)==d['source_event_counts']and min([H0,*radii])==H and d['source_limiter']=={'kind':'requested_slab','radius':str(H)}
    xp=[];yp=[];pc=Counter()
    for i,x in enumerate(xs):
        for k in range(i,len(xs)):
            y=xs[k]
            if y-x>2*R:pc['x_pair_span_prefilter']+=1;continue
            if y-R>=fr or x+R<=fl:pc['x_pair_reach_prefilter']+=1;continue
            xp.append((i,k))
    for j,y in enumerate(ys):
        for l in range(j,len(ys)):
            z=ys[l]
            if z-y>2*R:pc['y_pair_span_prefilter']+=1;continue
            if z-R>=ft or y+R<=fb:pc['y_pair_reach_prefilter']+=1;continue
            yp.append((j,l))
    counts=d['counts'];assert len(xp)==d['candidate_x_pairs']==23821 and len(yp)==d['candidate_y_pairs']==1849 and len(xp)*len(yp)==counts['boxes_visited']==44045029 and len(xs)*(len(xs)+1)//2*len(ys)*(len(ys)+1)//2==d['all_bboxes_before_prefilter']==49655025
    assert all(counts.get(k,0)==pc[k]for k in('x_pair_span_prefilter','x_pair_reach_prefilter','y_pair_span_prefilter','y_pair_reach_prefilter'))
    assert counts['active_boxes']+sum(v for k,v in counts.items()if k.startswith(('zero_mass','y_span','y_clip','x_mass','x_cap','x_clip')))==counts['boxes_visited']
    scale=lcm(*(z.denominator for z in[*xs,*ys,fl,fr,fb,ft,R]));T=lcm(*(F(p['w']).denominator for p in S));assert str(scale)==d['coordinate_scale']and str(T)==d['weight_scale']
    # Independent y-major prefix matrix (production uses x-major accumulation).
    xi={p:i for i,p in enumerate(xc)};yi={p:i for i,p in enumerate(ys)};grid=[[0]*len(xs)for _ in ys]
    for p in S:grid[yi[F(p['y'])]][xi[F(p['x']),p['v']]]+=int(F(p['w'])*T)
    P=[[0]*(len(xs)+1)]
    for row in grid:
        c=0;prev=P[-1];out=[0]
        for j,w in enumerate(row):c+=w;out.append(prev[j+1]+c)
        P.append(out)
    def mass(ids):
        i,k,j,l=ids;return F(P[l+1][k+1]-P[j][k+1]-P[l+1][i]+P[j][i],T)
    def expected_endpoint(q0,v,sgn,r2):
        n,t=isqrt(r2.numerator),isqrt(r2.denominator)
        return [str(q0+sgn*F(n,t)),str(v),0,'0']if n*n==r2.numerator and t*t==r2.denominator else[str(q0),str(v),sgn,str(r2)]
    records=d['box_records'];assert len(records)==len(d['boxes'])==counts['active_boxes']==25466;seen=set();xset,yset=set(xp),set(yp)
    for rec,box in zip(records,d['boxes']):
        ids=tuple(rec['bbox']);assert ids not in seen;seen.add(ids);i,k,j,l=ids;assert(i,k)in xset and(j,l)in yset;m=mass(ids);r2=min(R*R,m/(4*lam));assert m==F(rec['mass'])>0 and F(rec['radius_squared'])==r2 and rec['x_velocities']==[vs[i],vs[k]]
        assert (xs[k]-xs[i])**2<4*r2 and(ys[l]-ys[j])**2<4*r2 and max(F(0),fl-xs[i])**2<r2 and max(F(0),xs[k]-fr)**2<r2 and max(F(0),fb-ys[j])**2<r2 and max(F(0),ys[l]-ft)**2<r2
        expected=[expected_endpoint(xs[k],vs[k],-1,r2),expected_endpoint(xs[i],vs[i],1,r2),expected_endpoint(ys[l],0,-1,r2),expected_endpoint(ys[j],0,1,r2)];assert box==expected
    witnesses=[]
    for kind,rec in d['rejection_witnesses'].items():
        w=rec['witness'];ids=w['bbox'];i,k,j,l=ids;m=mass(ids);assert m==F(w['mass'])and list(map(F,w['x_endpoints']))==[xs[i],xs[k]]and w['x_velocities']==[vs[i],vs[k]]and list(map(F,w['y_endpoints']))==[ys[j],ys[l]]
        r2=min(R*R,m/(4*lam));D=xs[k]-xs[i];u=vs[k]-vs[i]
        if 'radius'in rec:
            if kind=='x_mass_span':q=lam*D*D-m;L=lam*abs(u)*(2*abs(D)+abs(u)*H)
            else:
                dist,speed=(fl-xs[i],-vs[i])if kind=='x_clip_left'else(xs[k]-fr,vs[k]);assert dist>=0;q=dist*dist-r2;L=abs(speed)*(2*abs(dist)+abs(speed)*H)
            assert q==F(w['q'])>0 and L==F(w['L'])>0 and F(rec['radius'])==min(H,q/(2*L))>=dt
        elif kind=='zero_mass':assert m==0
        elif kind=='y_span_permanent':assert(ys[l]-ys[j])**2>=4*r2
        elif kind=='y_clip_permanent':assert(ys[j]<=fb and(fb-ys[j])**2>=r2)or(ys[l]>=ft and(ys[l]-ft)**2>=r2)
        elif kind=='x_mass_span_permanent':assert u==0 and lam*D*D>=m
        elif kind.startswith('x_clip_'):
            dist,speed=(fl-xs[i],-vs[i])if kind.startswith('x_clip_left')else(xs[k]-fr,vs[k]);assert speed==0 and dist>=0 and dist*dist>=r2
        else:raise AssertionError('unhandled witness '+kind)
        witnesses.append(kind)
    assert d['min_empty_witness']=={'kind':'source_slab','radius':str(H),'witness':d['source_limiter']}
    out={'status':'passed','scope':'R39 input/catalogue/source-event/coverage counters/all active boxes and rejection summary witnesses; 44M rejected boxes not independently replayed','module_sha256':PIN,'slab_sha256':SLABPIN,'input_file_sha256':sha(a.input),'input_canonical_sha256':digest(base),'velocity_nonzero_labels':[130,131],'direction_file_sha256':sha(a.velocities),'s0':str(s0),'H':str(H),'delta_empty':str(dt),'tail_bound':str(tail),'source_count':len(S),'x_classes':len(xs),'y_classes':len(ys),'source_events_checked':sum(ec.values()),'candidate_pairs':[len(xp),len(yp)],'boxes_visited':counts['boxes_visited'],'all_active_box_masses_radii_endpoints_verified':len(records),'rejection_witness_kinds_verified':witnesses,'all_rejection_classes_individually_replayed':False,'script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
