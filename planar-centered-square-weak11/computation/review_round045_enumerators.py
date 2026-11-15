"""R45 eleven new slab instances; reuse reviewed R39 algorithm, no 44M replay."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import isqrt,lcm
from datetime import datetime,timezone
import argparse,gzip,hashlib,json,sys,time
ENUM='acede2c710655d0c179e7bea3ee9483271caaf37968227e7e9fce6442907056b'
INDEX='855c330b1d36dc212956eefcfba840a831aeabe17d7f8ba6dfa94b5466986a88'
BASE='ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def ceil(q):return (q.numerator+q.denominator-1)//q.denominator
def ser(x):return str(x)if isinstance(x,F)else{k:ser(v)for k,v in x.items()}if isinstance(x,dict)else[ser(v)for v in x]if isinstance(x,(list,tuple))else x

def check_slab(d,base,vel):
    assert d['code_sha256']==ENUM and d['input_file_sha256']==BASE and d['input_canonical_sha256']==digest(base)and d['velocities']==vel and d['direction_sha256']==digest(vel)
    px,py,lam=map(F,(base['Px'],base['Py'],base['lam']));M=sum(F(p[2])for p in base['pts']);s0,H0,H,dt=map(F,(d['s0'],d['requested_H'],d['H'],d['delta_empty']));R=F(d['tail']['radius']);fl,fr,fb,ft=map(F,d['domain']);assert 0<dt<=H<=H0 and R==10 and lam==1 and(fl,fr,fb,ft)==(0,px,0,py)
    tail=M/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(d['tail']['bound'])<lam and M==F(d['tail']['mass'])and F(d['tail']['lambda'])==lam
    S=[]
    for label,((xx,yy,ww),v)in enumerate(zip(base['pts'],vel)):
        x,y,w=map(F,(xx,yy,ww));xc=x+v*s0
        if not w:continue
        imin=ceil((fl-R-xc-v*H0)/px);z=(fr+R-xc+v*H0)/px;imax=z.numerator//z.denominator;jmin=ceil((fb-R-y)/py);z=(ft+R-y)/py;jmax=z.numerator//z.denominator
        for i,j in product(range(imin,imax+1),range(jmin,jmax+1)):S.append({'label':label,'period_shift':[i,j],'x':str(xc+i*px),'y':str(y+j*py),'v':v,'w':str(w)})
    assert S==d['source_catalog']and digest(S)==d['source_catalog_sha256'];xc=sorted({(F(p['x']),p['v'])for p in S});ys=sorted({F(p['y'])for p in S});xs=[z[0]for z in xc];vs=[z[1]for z in xc]
    assert [[str(x),v]for x,v in xc]==d['x_classes']and list(map(str,ys))==d['y_classes']and digest([d['x_classes'],d['y_classes']])==d['source_order_sha256']and all(x<y for x,y in zip(xs,xs[1:]))
    ec=Counter();h=H0;limiter={'kind':'requested_slab','radius':h}
    def event(gap,speed,kind,w):
        nonlocal h,limiter
        if not speed:ec[kind+('_persistent_equal'if not gap else'_parallel')]+=1;return
        assert gap;radius=abs(gap/speed)/2;assert radius>=H;ec[kind+'_nonpersistent']+=1
        if radius<h:h=radius;limiter={'kind':kind,'radius':h,'gap':gap,'relative_speed':speed,'witness':w}
    for i,(x,v)in enumerate(xc):
        for k in range(i+1,len(xc)):
            y,w=xc[k];event(y-x,w-v,'source_order',{'xclasses':[i,k]});event(y-x-2*R,w-v,'span_2R',{'xclasses':[i,k]})
        event(x-(fr+R),v,'reach_right',{'xclass':i,'source_x':x});event(x-(fl-R),v,'reach_left',{'xclass':i,'source_x':x})
    assert dict(ec)==d['source_event_counts']and h==H and ser(limiter)==d['source_limiter'];xp=[];yp=[];pc=Counter()
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
    counts=d['counts'];assert len(xp)==d['candidate_x_pairs']and len(yp)==d['candidate_y_pairs']and len(xp)*len(yp)==counts['boxes_visited']and len(xs)*(len(xs)+1)//2*len(ys)*(len(ys)+1)//2==d['all_bboxes_before_prefilter']
    assert all(counts.get(k,0)==pc[k]for k in('x_pair_span_prefilter','x_pair_reach_prefilter','y_pair_span_prefilter','y_pair_reach_prefilter'))and counts['active_boxes']+sum(v for k,v in counts.items()if k.startswith(('zero_mass','y_span','y_clip','x_mass','x_cap','x_clip')))==counts['boxes_visited']
    scale=lcm(*(z.denominator for z in[*xs,*ys,fl,fr,fb,ft,R]));T=lcm(*(F(p['w']).denominator for p in S));assert str(scale)==d['coordinate_scale']and str(T)==d['weight_scale'];xi={p:i for i,p in enumerate(xc)};yi={p:i for i,p in enumerate(ys)};grid=[[0]*len(xs)for _ in ys]
    for p in S:grid[yi[F(p['y'])]][xi[F(p['x']),p['v']]]+=int(F(p['w'])*T)
    P=[[0]*(len(xs)+1)]
    for row in grid:
        c=0;prev=P[-1];out=[0]
        for j,w in enumerate(row):c+=w;out.append(prev[j+1]+c)
        P.append(out)
    def mass(ids):
        i,k,j,l=ids;return F(P[l+1][k+1]-P[j][k+1]-P[l+1][i]+P[j][i],T)
    def endpoint(q0,v,sgn,r2):
        n,t=isqrt(r2.numerator),isqrt(r2.denominator)
        return[str(q0+sgn*F(n,t)),str(v),0,'0']if n*n==r2.numerator and t*t==r2.denominator else[str(q0),str(v),sgn,str(r2)]
    records=d['box_records'];assert len(records)==len(d['boxes'])==counts['active_boxes'];seen=set();xset,yset=set(xp),set(yp)
    for rec,box in zip(records,d['boxes']):
        ids=tuple(rec['bbox']);assert ids not in seen;seen.add(ids);i,k,j,l=ids;assert(i,k)in xset and(j,l)in yset;m=mass(ids);r2=min(R*R,m/(4*lam));assert m==F(rec['mass'])>0 and F(rec['radius_squared'])==r2 and rec['x_velocities']==[vs[i],vs[k]]
        assert(xs[k]-xs[i])**2<4*r2 and(ys[l]-ys[j])**2<4*r2 and max(F(0),fl-xs[i])**2<r2 and max(F(0),xs[k]-fr)**2<r2 and max(F(0),fb-ys[j])**2<r2 and max(F(0),ys[l]-ft)**2<r2
        assert box==[endpoint(xs[k],vs[k],-1,r2),endpoint(xs[i],vs[i],1,r2),endpoint(ys[l],0,-1,r2),endpoint(ys[j],0,1,r2)]
    radii=[H];witnesses=[]
    for kind,rec in d['rejection_witnesses'].items():
        w=rec['witness'];ids=w['bbox'];i,k,j,l=ids;m=mass(ids);assert m==F(w['mass'])and list(map(F,w['x_endpoints']))==[xs[i],xs[k]]and w['x_velocities']==[vs[i],vs[k]]and list(map(F,w['y_endpoints']))==[ys[j],ys[l]];r2=min(R*R,m/(4*lam));D=xs[k]-xs[i];u=vs[k]-vs[i]
        if 'radius'in rec:
            if kind=='x_mass_span':q=lam*D*D-m;L=lam*abs(u)*(2*abs(D)+abs(u)*H)
            else:
                dist,speed=(fl-xs[i],-vs[i])if kind=='x_clip_left'else(xs[k]-fr,vs[k]);assert dist>=0;q=dist*dist-r2;L=abs(speed)*(2*abs(dist)+abs(speed)*H)
            assert q==F(w['q'])>0 and L==F(w['L'])>0 and F(rec['radius'])==min(H,q/(2*L));radii.append(F(rec['radius']))
        elif kind=='zero_mass':assert m==0
        elif kind=='y_span_permanent':assert(ys[l]-ys[j])**2>=4*r2
        elif kind=='y_clip_permanent':assert(ys[j]<=fb and(fb-ys[j])**2>=r2)or(ys[l]>=ft and(ys[l]-ft)**2>=r2)
        elif kind.startswith('x_mass_span'):assert u==0 and lam*D*D>=m
        elif kind.startswith('x_cap_span'):assert u==0 and D*D>=4*R*R
        elif kind.startswith('x_clip_'):
            dist,speed=(fl-xs[i],-vs[i])if kind.startswith('x_clip_left')else(xs[k]-fr,vs[k]);assert speed==0 and dist>=0 and dist*dist>=r2
        else:raise AssertionError('unhandled witness '+kind)
        witnesses.append(kind)
    assert dt==min(radii);z=d['min_empty_witness']
    if z['kind']=='source_slab':assert z=={'kind':'source_slab','radius':str(H),'witness':d['source_limiter']}
    else:assert z=={'kind':z['kind'],**d['rejection_witnesses'][z['kind']]}and F(z['radius'])==dt
    return {'source_count':len(S),'x_classes':len(xs),'y_classes':len(ys),'source_events_checked':sum(ec.values()),'candidate_pairs':[len(xp),len(yp)],'boxes_visited':counts['boxes_visited'],'active_boxes_verified':len(records),'rejection_summary_kinds_verified':witnesses,'H':str(H),'delta_empty':str(dt),'tail_bound':str(tail)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--geometry-topic',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;rr=t/'computation/results/round045';indexp=raw/'results_index.json';assert sha(indexp)==sha(rr/'results_index.json')==INDEX;idx=read(indexp);assert len(idx['runs'])==11
    bp=t/'computation/inputs/r041_final_0036.json';b8p=t/'computation/inputs/r008_final_0611.json';base,b8=read(bp),read(b8p);assert sha(bp)==BASE and sha(b8p)=='6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910';assert bp.read_bytes()==(raw.parent/'round041/cross-group-probe/candidate_0036.json').read_bytes()
    for name in('primary_old_R10.json','independent_old_R10.json'):assert read(t/'computation/results/round041'/name)['input_sha256']==BASE
    assert len(b8['pts'])==51 and len(base['pts'])==204 and F(base['Px'])==2*F(b8['Px'])and F(base['Py'])==2*F(b8['Py'])and all(F(p[2])>0 for p in base['pts'])and sum(F(p[2])for p in base['pts'])==4*sum(F(p[2])for p in b8['pts'])
    assert F(b8['pts'][26][0])==F(b8['pts'][27][0]);P8=F(b8['Px']);assert P8==F(328734960629,100000000000);modulep=Path(read(raw/'local-slope/run_manifest.json')['commands'][0][1]);assert sha(modulep)==ENUM
    archive=read(rr/'source-manifest.json');assert sha(rr/'source-manifest.json')=='7b444c5046d07b8e83c4101c02d6e5eb46e8cbcaca8d104df584a12130673058';afiles=0
    for z in archive['files']:
        src=Path(z['source']);dst=t.parent/z['destination'];assert src.stat().st_size==z['raw_size']and sha(src)==z['raw_sha256']and dst.stat().st_size==z['archived_size']and sha(dst)==z['archived_sha256']
        if dst.suffix=='.gz':assert hashlib.sha256(gzip.decompress(dst.read_bytes())).hexdigest()==z['raw_sha256']
        else:assert src.read_bytes()==dst.read_bytes()
        afiles+=1
    centers={'local-slope':F(-1,1280000),'local-slope-q0':F(-47,57600000),'local-slope-q2':F(-47,57600000),'local-slope-q3':F(-47,57600000),'local-slope-all':F(-47,57600000),'negative-slope-q0':F(-1,1280000),'negative-slope-q1':F(1,19200000),'negative-slope-q2':F(-1,1280000),'negative-slope-q3':F(-1,1280000),'negative-slope-all':F(-43,57600000),'negative-slope-all-next':F(-41,57600000)};assert {z['name']for z in idx['runs']}==set(centers);runs=[]
    for z in idx['runs']:
        name=z['name'];rp=raw/name;m=read(rp/'run_manifest.json');vel=read(rp/'velocities.json');d=read(rp/'slab.json');slope=read(rp/'slope.json')
        for fn,h in z['hashes'].items():assert sha(rp/fn)==h
        q=1 if name=='local-slope'else None if 'all'in name else int(name[-1]);labels=[51*k+j for k in(range(4)if q is None else[q])for j in(26,27)];expected=[int(j in labels)for j in range(204)];assert vel==expected and d['direction_file_sha256']==sha(rp/'velocities.json')==m['velocity_sha256'];tau=centers[name];HH=F(1,100000000);physical=tau*P8;physicalH=HH*P8
        assert F(z['normalized_s0'])==F(m['normalized_s0'])==tau and F(z['normalized_requested_H'])==F(m['normalized_H'])==HH and F(m['physical_s0'])==F(d['s0'])==physical and F(m['physical_H'])==F(d['requested_H'])==physicalH
        assert m['input_sha256']==BASE and m['base_sha256']==sha(b8p)and m['enum_sha256']==ENUM and m['sweep_sha256']=='a4c44c806c533bcefb6fbdcda05ee98918caa4fbf585a9ffc64eb5727dd09f56'
        runner='run_r45_slope.py'if name=='local-slope'else 'run_r45_'+('slope_'+name.removeprefix('local-slope-')if name.startswith('local-slope-')else'negative_'+name.removeprefix('negative-slope-').replace('-','_'))+'.py';runnerp=a.geometry_topic/'computation'/runner;assert sha(runnerp)==m['runner_sha256']
        cmd=m['commands'][0];assert cmd[cmd.index('--s0')+1]==str(physical)and cmd[cmd.index('--H')+1]==str(physicalH)and cmd[cmd.index('--R')+1]=='10'and sha(Path(cmd[2]))==BASE and Path(cmd[cmd.index('--output')+1])==rp/'slab.json'
        result=check_slab(d,base,vel);assert z['counts']==d['counts']and z['tail']==d['tail']and slope['slab_sha256']==z['hashes']['slab.json']and slope['slab_code_sha256']==ENUM and F(slope['s0'])==physical and F(z['physical_delta'])==F(slope['delta'])>0 and F(slope['delta'])<=F(d['delta_empty'])and slope['beta_lower']==z['beta_lower']and slope['beta_upper']==z['beta_upper']
        result.update({'name':name,'labels':labels,'normalized_s0':str(tau),'physical_s0':str(physical),'normalized_requested_H':str(HH),'physical_requested_H':str(physicalH),'hashes':z['hashes'],'runner':str(runnerp),'runner_sha256':sha(runnerp),'slope_beta_sign_recorded':1 if F(z['beta_lower'])>0 else-1 if F(z['beta_upper'])<0 else 0,'slope_recomputed_by_this_script':False});runs.append(result);print(name,'catalogue/active/margins passed',flush=True)
    report=t/'rounds/045.md';assert sha(report)=='0017a70bee0bf505f37735150b363bf6fdb27d84f8ba0828f37c5b1677cc99d7'
    out={'status':'passed','scope':'eleven frozen R45 new source/slab instances; R39 frozen enumerator review reused; no large rejected-bbox loop or slope sort/merge replay','source_module_path':str(modulep),'source_module_sha256':sha(modulep),'index_sha256':INDEX,'R41_input_sha256':BASE,'R41_input_canonical_sha256':digest(base),'R8_base_period':str(P8),'archive_files_identity_verified':afiles,'runs':runs,'total_active_boxes_verified':sum(z['active_boxes_verified']for z in runs),'total_saved_bbox_visits':sum(z['boxes_visited']for z in runs),'all_rejected_boxes_individually_replayed':False,'bracket_refine_runs_included':False,'report_sha256':sha(report),'source_manifest_sha256':sha(rr/'source-manifest.json'),'script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='runs'},indent=2))
if __name__=='__main__':main()
