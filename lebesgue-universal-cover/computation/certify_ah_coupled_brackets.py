"""Shared contact-order polytopes at three directions and common angle brackets on dangerous AH bins."""
import argparse,contextlib,itertools,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb
from certify_ah_projected_jets import geometry


def support_bound(ns,j,a,b):
    A,add,mul,dot,cross,sub,reflection=(ns[k] for k in ('A','add','mul','dot','cross','sub','reflection')); zero=A(0); E=ns['deficit_by_bin'][j]
    v=ns['centered_points'](a,b); edges=[sub(v[k],v[l]) for k,l in ns['orders']]; z=(zero,zero)
    ca,cb=[ns['chart'](add(ns['z'],mul(A(w),ns['D']))) for w in (a,b)]; l,t=ca[0].union(cb[0]),ca[1].union(cb[1]); lp=ns['coord_hull'](ns['psi'],l); jets={}
    for k,ell,T in ((4,l,t),(2,lp,-t)):
        x,lc,nc=ns['chart_bounds'](ell+(-A(F(2,3))*E).union(A(F(2,3))*E),T+(-E).union(E));u,w=ns['derivs'](x,lc,nc);jets[k]=(u,w) if k==4 else (reflection(u),reflection(w))
    def choices(c,u):
        cr=[cross(e,u) for e in edges];result=[]
        for k,(vertex,arc) in enumerate(ns['sector']):
            m=(k+1)%10
            if cr[k]<0 or cr[m]>0:continue
            sign=-1 if arc else 1; weight=sign*c; det=cross(edges[k],edges[m]); assert det>0
            l,r=-cr[m]/det,cr[k]/det; left,right=ns['orders'][k][0 if arc else 1],ns['orders'][m][0 if arc else 1]
            if not weight>0:result.append([(vertex,mul(sign*weight,u))])
            if not weight<0:result.append([(left,mul(sign*weight*l,edges[k])),(right,mul(sign*weight*r,edges[m]))])
        assert result
        return result
    def extremal(vals):
        return max(x.lower() for x in vals).union(max(x.upper() for x in vals))
    def support_chain(coefs,lo,hi):
        return extremal([sum(c*(lo if i<k else hi) for i,c in enumerate(coefs)) for k in range(len(coefs)+1)])
    def event(action):
        vec=[z]*5
        for i,u in action:vec[i]=add(vec[i],u)
        ll,lr=dot(vec[0],ns['el']),dot(vec[1],reflection(ns['el'])); dt=dot(vec[3],ns['en']); a2,t2=(dot(vec[2],q) for q in jets[2]);a4,t4=(dot(vec[4],q) for q in jets[4]);
        def value(x,eta,xi,first,de,partner):
            return support_chain((x/A('.95'),-eta,-xi),zero,A(F(2,3)))+support_chain((first,de,-partner),A(-1),A(1))
        return extremal((value(ll,lr,a2,t4,dt,t2),value(lr,ll,a4,t2,-dt,t4)))
    def bound(terms,detail=False):
        groups=[choices(c,u) for c,u in terms]; common=extremal([event(sum(g,[])) for g in itertools.product(*groups)])
        if detail:
            # For singleton cones (the case used for the strict incompatibility example), both sides are finite maxima of explicit linear functions.
            sep=sum(extremal([event(action) for action in group]) for group in groups)
            return common,sep,[len(g) for g in groups]
        return common.upper()
    def incompatibility(terms,coattain=False):
        groups=[choices(c,u) for c,u in terms]
        assert all(len(g)==1 for g in groups)
        actions=[g[0] for g in groups]
        # State coordinates are (L0,L1,ell2,t2,DE,ell4,t4) in their intrinsic axes.
        states=[]
        for swap,lk,tk in itertools.product((False,True),range(4),range(4)):
            lx,le,lp=[F(0) if k<lk else F(2,3) for k in range(3)]
            x,eta,xi=lx/F('.95'),-le,-lp
            first,de,npartner=[F(-1) if k<tk else F(1) for k in range(3)];partner=-npartner
            coeff=(x,eta,xi,partner,de,0,first) if not swap else (eta,x,0,first,-de,xi,partner)
            states.append((swap,lk,tk,coeff))
        axes=[ns['el'],reflection(ns['el']),*jets[2],ns['en'],*jets[4]]
        vertex_axes=((0,),(1,),(2,3),(4,),(5,6))
        def value(actions,coefficients):
            directions=[z for _ in axes]
            for action,state in zip(actions,coefficients,strict=True):
                for index,u in action:
                    for m in vertex_axes[index]:
                        if state[m]: directions[m]=add(directions[m],mul(A(state[m]),u))
            return sum(dot(u,a) for u,a in zip(directions,axes,strict=True))
        if coattain:
            for shared in states:
                good=True
                for action in actions:
                    for other in states:
                        diff=tuple(x-y for x,y in zip(shared[3],other[3]))
                        if value([action],[diff]).lower()<0: good=False;break
                    if not good:break
                if good:return shared[:3]
            return None
        independent=[max(states,key=lambda state:value([action],[state[3]]).mid()) for action in actions]
        diffs=[]
        for shared in states:
            delta=[tuple(x-y for x,y in zip(st[3],shared[3])) for st in independent]
            diffs.append(value(actions,delta).lower())
        return min(diffs),[st[:3] for st in independent]
    return E,v,edges,bound,incompatibility


def scan(ns,j,begin,end,left,right,angular=16):
    A,rotate,trig=ns['A'],ns['rotate'],ns['trig']; n0,n1=ns['n0'],ns['n1'];turn,sign,orient=(0,1,-1) if j==34 else (300,1,1); bestroot=bestedge=A(1); rootwhere=edgewhere=None; gain=A(0); gainrow=None
    for i in range(begin,end):
        a=F(j,100)+F(i,7500); b=a+F(1,7500);E,v,edges,bound,_=support_bound(ns,j,a,b)
        for s,phi in ((orient,left),(-orient,right)):
            terms=[(A(-s*(-1)**k),rotate(n0,trig(A(phi+60*k)))) for k in range(3)]
            common,sep,opts=bound(terms,True); slack=(s*ns['F_at'](v,edges,phi)-E*common.upper()).lower()
            if slack < bestroot: bestroot,rootwhere=slack,[i,str(phi)]
            saving=(E*(sep.lower()-common.upper())).lower()
            if opts==[1,1,1] and saving>gain: gain,gainrow=saving,{'subcell':i,'phi':str(phi),'sign':s,'separate_per_E':str(sep),'joint_per_E':str(common),'E':str(E)}
        for p in range(angular):
            al=left+(right-left)*F(p,angular);bl=left+(right-left)*F(p+1,angular);ms=ns['placement'](v,edges,al,bl,turn,sign);sn=trig(ns['box'](al-turn,bl-turn));u0=rotate(n0,sn);u1=rotate((n1[0],sign*n1[1]),sn)
            for m,(name,u,_) in zip(ms,ns['rows'],strict=True):
                rough=(m-A('3.21')*E).lower()
                if rough>A('.00025'): val=rough
                else:
                    d0=u[0]-u[1]/ns['r3'];d1=2*u[1]/ns['r3'];ur=rotate((u[0],sign*u[1]),sn)
                    val=(m-E*bound([(A(1),ur),(-d0,u0),(-d1,u1)])).lower()
                if val<bestedge:bestedge,edgewhere=val,[i,p,name]
    return {'reference_bin':j,'subcell_halfopen_indices':[begin,end],'closed_w_endpoints':[str(F(j,100)+F(k,7500)) for k in (begin,end)],'common_bracket_degrees':[str(left),str(right)],'angle_cells_per_parameter_cell':angular,'root_slack_lower':str(bestroot),'worst_root':rootwhere,'placement_slack_lower':str(bestedge),'worst_pose':edgewhere,'certified_strict_conflict_gain':str(gain),'separate_vs_shared_example':gainrow}


def run(archive,inp,output):
    settings=json.loads(inp.read_bytes()); extra=[tuple(map(F,r)) for r in settings['lines']];ns,source=geometry(archive,extra); rows=[]
    for j,begin,end,left,right in settings['brackets']:
        r=scan(ns,j,begin,end,F(left),F(right));rows.append(r);print(j,begin,end,r['root_slack_lower'],r['placement_slack_lower'],r['certified_strict_conflict_gain'],flush=True)
    assert all(arb(r['root_slack_lower'])>0 and arb(r['placement_slack_lower'])>0 for r in rows)
    for j in (31,34):
        selected=sorted(r['subcell_halfopen_indices'] for r in rows if r['reference_bin']==j);assert selected[0][0]==0 and selected[-1][1]==75 and all(s[1]==t[0] for s,t in zip(selected,selected[1:]))
    conflicts=[]
    for j,i,phi,name in settings.get('conflict_examples',[]):
        a=F(j,100)+F(i,7500); b=a+F(1,7500); E,v,edges,bo,inc=support_bound(ns,j,a,b)
        A=ns['A'];turn=0 if j==34 else 300;sn=ns['trig'](A(F(phi)-turn)); rotate=ns['rotate']; u=next(u for label,u,_ in ns['rows'] if label==name); d0=u[0]-u[1]/ns['r3'];d1=2*u[1]/ns['r3']
        terms=[(A(1),rotate(u,sn)),(-d0,rotate(ns['n0'],sn)),(-d1,rotate(ns['n1'],sn))]
        gain,states=inc(terms); assert gain>0,(j,gain); common,separate,_=bo(terms,True)
        conflicts.append({'bin':j,'subcell':i,'w_endpoints':[str(a),str(b)],'phi_degrees':phi,'edge':name,'model_independent_states':states,'verified_excess_over_every_common_state_per_E':str(gain),'excess_after_multiplying_E':str((E*gain).lower()),'shared_error_upper':str(common.upper())})
    common_examples=[]
    for j,i,phi,name in settings.get('common_state_examples',[]):
        a=F(j,100)+F(i,7500);b=a+F(1,7500);E,v,edges,bo,inc=support_bound(ns,j,a,b);A=ns['A'];turn=0 if j==34 else 300;sn=ns['trig'](A(F(phi)-turn));rot=ns['rotate'];u=next(u for label,u,_ in ns['rows'] if label==name);d0=u[0]-u[1]/ns['r3'];d1=2*u[1]/ns['r3'];terms=[(A(1),rot(u,sn)),(-d0,rot(ns['n0'],sn)),(-d1,rot(ns['n1'],sn))];state=inc(terms,True);assert state is not None,(j,i,phi,name)
        common_examples.append({'bin':j,'subcell':i,'phi_degrees':phi,'edge':name,'one_joint_maximizing_relaxed_state':state})
    output.write_text(json.dumps({'input':settings,'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'geometry_sha256':sha256(source.encode()).hexdigest(),'geometry_and_placement_bits':128,'python':sys.version,'flint':ns['flint_version'],'brackets':rows,'strict_nonattainment_examples':conflicts,'common_relaxed_extremizers':common_examples},indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
