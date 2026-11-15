#!/usr/bin/env python3
"""Independent exact finite-instance audit of R52 position/mass bbox margins."""
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,os,platform,random,sys
from pathlib import Path

def load(path,name):
    sp=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(implementation):
    mp=Path(implementation); impl=load(mp,'audit_impl'); en=load(mp.with_name('position_slab_enumerator.py'),'audit_en')
    rng=random.Random(520052); reports=[]; tested=active_samples=nonempty_samples=0; rejected_centers=[]
    for case in range(60):
        n=2+case%4
        xx=rng.sample(range(-24,25),n); yy=rng.sample(range(-24,25),n)
        pts=[[str(F(x,7)),str(F(y,11)),str(F(rng.randrange(5,31),5))] for x,y in zip(xx,yy)]
        vx=[rng.randrange(2) for _ in pts]; mu=[F(rng.randrange(-6,7),7) for _ in pts]
        periodic=case>=40
        data={'lam':'1','pts':pts}; domain=list(map(str,[0,20,0,21] if periodic else [-10,10,-10,10]))
        if periodic:data.update(Px='20',Py='21')
        slab=en.enumerate_slab(data,vx,F(1,13),F(1,100),5,domain=domain,periodic=periodic)
        out=impl.build_position_mass_slab(data,mu,vx,slab,T=F(1,100),t0=F(1,2))
        rho,H,T=map(F,(out['rho'],out['H'],out['T'])); records={tuple(r['bbox']):r for r in out['box_records']}
        boxes=dict(zip((tuple(r['bbox']) for r in out['box_records']),out['boxes']))
        xs=[(F(x),int(v)) for x,v in slab['x_classes']];ys=list(map(F,slab['y_classes']))
        sources=slab['source_catalog']; fl,fr,fb,ft=map(F,out['domain']); R=F(out['R']);lam=F(out['lambda'])
        # Direct source summation, independent of the two implementation prefix arrays.
        def bound_ok(D,j,q): return D<0 or D*D<j*j*q
        for i in range(len(xs)):
          for k in range(i,len(xs)):
           for j in range(len(ys)):
            for l in range(j,len(ys)):
             key=(i,k,j,l)
             subset=[s for s in sources if xs[i][0]<=F(s['x'])<=xs[k][0] and ys[j]<=F(s['y'])<=ys[l]]
             m=sum((F(s['w']) for s in subset),F(0));dm=sum((mu[s['label']] for s in subset),F(0))
             tested+=1
             if key in records:
                rec=records[key]; assert F(rec['m0'])==m and F(rec['mu'])==dm
                assert F(rec['q0'])==m/(4*lam) and F(rec['q1'])==dm/(4*lam)
                for e in boxes[key]:
                    if int(e[2]):assert F(e[3])==m/(4*lam) and F(e[4])==dm/(4*lam)
             for fz in [-1,F(-1,3),0,F(1,3),1]:
              for ftau in [-1,0,1]:
                dz=fz*rho*H;dt=ftau*rho*T;a=xs[i][0]+xs[i][1]*dz;b=xs[k][0]+xs[k][1]*dz;c=ys[j];d=ys[l]
                q=min(R*R,(m+dm*dt)/(4*lam))
                positive=(q>0 and bound_ok(b-a,2,q) and bound_ok(d-c,2,q) and bound_ok(fl-a,1,q) and bound_ok(b-fr,1,q) and bound_ok(fb-c,1,q) and bound_ok(d-ft,1,q))
                if positive:
                    nonempty_samples+=1
                    assert key in records,('omitted box activated',case,key,dz,dt)
                if key in records:active_samples+=1
        reports.append({'case':case,'periodic':periodic,'sources':len(sources),'rho':str(rho),'candidates':len(records),'visited':out['counts']['boxes_visited']})
    # Check integer lifting independently for general positive rational parameters.
    lifts=0
    for _ in range(1000):
        S,W,ln,ld=[rng.randrange(1,10) for _ in range(4)];B=4*W*ln
        D,u,j=rng.randrange(1,100),rng.randrange(-3,4),rng.choice([1,2]);M,MU=rng.randrange(1,30),rng.randrange(-20,21)
        H,T=F(rng.randrange(1,9),rng.randrange(2,20)),F(rng.randrange(1,9),rng.randrange(2,20))
        RI=rng.randrange(1,30);K=min(M*ld*S*S,RI*RI*B);Q=D*D*B-j*j*K
        if Q<=0:continue
        hn,hd=H.numerator,H.denominator;tn,td=T.numerator,T.denominator
        Lint=abs(u)*(2*abs(D)*hd+abs(u)*hn*S)*(B*S*hn*td)+j*j*abs(MU)*(S*S*ld*tn*hd*hd)
        q=F(Q,j*j*S*S*B); v=F(u,j);d0=F(D,j*S)
        L=abs(v)*(2*abs(d0)+abs(v)*H)*H+F(abs(MU)*ld,B)*T
        assert L==F(Lint,j*j*S*S*B*hd*hd*td)
        if L: assert q/(2*L)==F(Q*hd*hd*td,2*Lint)
        lifts+=1
    # Independent min-Lipschitz tests include crossing and persistent cap regimes.
    lipschitz=crossings=0
    for _ in range(1000):
        cap=F(rng.randrange(1,30),rng.randrange(1,10)); x=F(rng.randrange(0,100),rng.randrange(1,10)); y=F(rng.randrange(0,100),rng.randrange(1,10))
        assert abs(min(cap,x)-min(cap,y))<=abs(x-y)
        lipschitz+=1;crossings+=int(min(x,y)<cap<max(x,y))
    # At Q0=0, u=0, and lower mass above cap, the failed strict condition
    # remains equality for every dt even with nonzero MU. These are algebra
    # tests, independent of whether this branch occurs in a legal tail slab.
    cap_equalities=0
    for cap in [F(1,4),F(9,4),F(25,9)]:
      for dm in [F(-7,3),F(2,5)]:
        T=F(1,10);m=cap+abs(dm)*T
        for tau in [-1,F(-1,3),0,F(1,3),1]:
            assert min(cap,m+dm*tau*T)==cap
        cap_equalities+=1
    # Fixed-geometry affine mass path cannot activate an endpoint-empty box.
    interpolation=0
    for _ in range(1000):
        d=F(rng.randrange(0,30),rng.randrange(1,10));R=F(rng.randrange(1,30),rng.randrange(1,10));lam=F(rng.randrange(1,10),rng.randrange(1,10))
        m0=F(rng.randrange(0,100),rng.randrange(1,10));m1=F(rng.randrange(0,100),rng.randrange(1,10))
        active=lambda m:R>d and m>4*lam*d*d
        ends=[active(m0),active(m1)]
        for tau in [F(0),F(1,7),F(1,3),F(1,2),F(6,7),F(1)]:
            inside=active((1-tau)*m0+tau*m1)
            assert not inside or any(ends)
            if ends[0]==ends[1]:assert inside==ends[0]
            interpolation+=1
    # Zero-margin mass activation is retained, never silently omitted.
    ydata={'lam':'1','pts':[['0','0','9'],['0','5','16']]}
    yslab=en.enumerate_slab(ydata,[0,0],0,F(1,10),3,domain=['-4','4','-4','9'],periodic=False)
    yo=impl.build_position_mass_slab(ydata,[1,0],[0,0],yslab,T=F(1,1000))
    equal=[r for r in yo['box_records'] if r['center_status']=='empty_equality_wall']
    assert any(F(r['mu'])==1 for r in equal)
    return {'status':'passed','scope':'exact finite test instances and scaling identities, not the R52 full catalogue or endpoint-order certificate','seed':520052,'instance_cases':reports,'direct_bbox_checks':tested,'candidate_parameter_samples':active_samples,'nonempty_parameter_samples':nonempty_samples,'integer_lifting_checks':lifts,'min_lipschitz_checks':lipschitz,'cross_cap_lipschitz_checks':crossings,'abstract_persistent_cap_equality_checks':cap_equalities,'fixed_geometry_affine_mass_checks':interpolation,'zero_margin_mass_activation_retained':True,'implementation_sha256':digest(mp),'enumerator_sha256':digest(mp.with_name('position_slab_enumerator.py')),'audit_code_sha256':digest(__file__),'python':sys.version,'executable':sys.executable,'platform':platform.platform()}

def real_audit(path):
    p=Path(path);out=json.loads(p.read_text());cmd=out['command']
    paths={'center_data':Path(cmd[2]),'mass_velocity':Path(cmd[cmd.index('--mass-velocity')+1]),'x_velocity':Path(cmd[cmd.index('--x-velocity')+1]),'source_slab':Path(cmd[cmd.index('--source-slab')+1])}
    for name,ip in paths.items():assert digest(ip)==out['input_file_sha256'][name]
    data,mus,vxs,slab=(json.loads(paths[k].read_text()) for k in ['center_data','mass_velocity','x_velocity','source_slab'])
    xs=[(F(x),int(v)) for x,v in slab['x_classes']];ys=list(map(F,slab['y_classes']));ix={x:i for i,x in enumerate(xs)};iy={y:j for j,y in enumerate(ys)}
    W=int(out['weight_derivative_common_scale']);lam=F(out['lambda']);R=F(out['R']);H=F(out['H']);T=F(out['T']);rho=F(out['rho'])
    mu=list(map(F,mus));cols0=[[0]*(len(ys)+1) for _ in xs];cols1=[[0]*(len(ys)+1) for _ in xs]
    for s in slab['source_catalog']:
        i=ix[F(s['x']),s['v']];j=iy[F(s['y'])]+1;cols0[i][j]+=int(F(s['w'])*W);cols1[i][j]+=int(mu[s['label']]*W)
    for c0,c1 in zip(cols0,cols1):
        for j in range(1,len(ys)+1):c0[j]+=c0[j-1];c1[j]+=c1[j-1]
    def masses(key):
        i,k,j,l=key
        return F(sum(c[l+1]-c[j] for c in cols0[i:k+1]),W),F(sum(c[l+1]-c[j] for c in cols1[i:k+1]),W)
    records=out['box_records'];keys={tuple(r['bbox']) for r in records};oldkeys={tuple(r['bbox']) for r in slab['box_records']}
    assert keys==oldkeys and len(keys)==len(records)==26327
    nonzero=0
    for r,box in zip(records,out['boxes']):
        m,dm=masses(r['bbox']);nonzero+=int(dm!=0)
        assert F(r['m0'])==m and F(r['mu'])==dm
        assert F(r['q0'])==m/(4*lam) and F(r['q1'])==dm/(4*lam)
        assert r['radius_branch']=='root' and F(r['q_full_domain_lower'])==m/(4*lam)-abs(dm)/(4*lam)*rho*T
        assert F(r['q_full_domain_lower'])>0 and F(r['q_full_domain_upper'])<R*R
        i,k,j,l=r['bbox'];expected=[(xs[k][0],xs[k][1],-1),(xs[i][0],xs[i][1],1),(ys[l],0,-1),(ys[j],0,1)]
        for e,(a,v,sign) in zip(box,expected):
            assert list(map(F,[e[0],e[1],e[3],e[4]]))==[a,F(v),m/(4*lam),dm/(4*lam)] and e[2]==sign
    fl,fr,fb,ft=map(F,out['domain']);witnesses=0
    for kind,w in out['rejection_witnesses'].items():
        i,k,j,l=w['bbox'];m,dm=masses(w['bbox']);a,va=xs[i];b,vb=xs[k];c,d=ys[j],ys[l]
        D,V={'x_span':((b-a)/2,F(vb-va,2)),'y_span':((d-c)/2,F(0)),'x_clip_left':(fl-a,F(-va)),'x_clip_right':(b-fr,F(vb)),'y_clip_below':(fb-c,F(0)),'y_clip_above':(d-ft,F(0))}[kind]
        q=D*D-min(R*R,m/(4*lam));L=abs(V)*(2*abs(D)+abs(V)*H)*H+abs(dm)/(4*lam)*T
        assert F(w['m0'])==m and F(w['mu'])==dm and F(w['d0'])==D and F(w['d_velocity'])==V and F(w['q'])==q and F(w['L'])==L
        r=F(w['rho_bound']);assert q>0 and q-r*L>=q/2
        if w['clip_distance_positivity_required']:assert D-abs(V)*r*H>0
        witnesses+=1
    counts=out['counts'];assert counts['candidate_boxes']+sum(v for k,v in counts.items() if k.startswith('reject_'))==counts['boxes_visited']==45633320
    assert F(out['uniform_tail']['bound'])<lam and rho==1
    return {'status':'passed','scope':'all actual candidate mass slopes/endpoints via independent column-sum method; six rejection witnesses and ledger, not full omitted-box re-enumeration or endpoint ordering','full_output_sha256':digest(p),'full_output_bytes':p.stat().st_size,'candidates':len(records),'nonzero_mass_derivative_candidates':nonzero,'rejection_witness_checks':witnesses,'visited_ledger':counts['boxes_visited'],'rho':str(rho),'H':str(H),'T':str(T),'active_keys_equal_source_slab':True,'input_file_sha256':{k:digest(v) for k,v in paths.items()}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--implementation',required=True);ap.add_argument('--output',required=True);ap.add_argument('--real-output');args=ap.parse_args()
    assert os.path.ismount('/Volumes/data'); out=run(args.implementation)
    if args.real_output:out['real_catalogue_audit']=real_audit(args.real_output)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('instance_cases','python','platform')},indent=2))
