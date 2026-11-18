"""Balance rational facet normals for mixed-area and whole-neighborhood exclusions."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys

from flint import arb, ctx

K = 1024
EDGE_COUNT = 24
AREA_LOWER = F('.83735')
WALL_TARGET = F('.83740')
WALL_RADIUS = F('.00035')
JOINT_GATE = F('.00022')
PI_UP = F(22,7)
MIXED_GOALS = {
    ('new_P30','Q11'):F('.00025'),
    ('old_P30','Q11'):F('.00024'),
    ('escape108','Q11'):F('.00030'),
    ('contact102','P30'):F('.00014'),
    ('contact102','Q11'):F('.00010'),
}
WALL_CASES = ('new_P30','old_P30','escape108')
CENTERED_P30_POSES = {
    'new_P30': ('1/3','.012538374171301372','-.01767909572630291'),
    'old_P30': ('1/3','.009649153645133033','-.019734111618847052'),
    'escape108': ('1/3','.007806446413339705','-.019114591132998335'),
}


def ar(x):
    q = F(x)
    return arb(q.numerator) / q.denominator


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def rotate(p,c,s):
    x,y = p
    return x*c-y*s,x*s+y*c


def upper(x, scale=10**12):
    v = F(round(float(x.mid())*scale)+2,scale)
    assert x < ar(v)
    return v


class Host:
    def __init__(self, pose, pi):
        from certify_nine import S, append_point, arc_body, hull

        nodes=[]
        for j in range(8*K):
            angle=pi*j/(4*K)
            append_point(nodes,angle.cos()/2,angle.sin()/2)
        for n,index in ((3,-1),(5,2),(7,5),(9,8)):
            rho=ar(pose[index]) if index>=0 else arb(0)
            arc_body(nodes,K,pi,n,rho,*map(ar,pose[index+1:index+3]))
        poly=hull(nodes)
        self.material=poly
        self.den=S+4
        assert max(dot(x,x) for x in poly)*25 < 9*self.den**2
        area_twice=0
        edges=[]
        for p,q in zip(poly,poly[1:]+poly[:1],strict=True):
            normal=(q[1]-p[1],p[0]-q[0])
            sq=dot(normal,normal)
            margin=dot(normal,p)
            assert margin>0 and sq>0
            area_twice+=det(p,q)
            edges.append((sq,normal,margin))
        self.area=F(area_twice,2*self.den**2)
        self.hull_hash=sha256(json.dumps(poly,separators=(',',':')).encode()).hexdigest()
        self.edges=sorted(edges,reverse=True)[:EDGE_COUNT]
        self.normals=[n for _,n,_ in self.edges]
        self.nums=[(ar(F(n[0],self.den)),ar(F(n[1],self.den))) for n in self.normals]
        self.eta=F(6,S)+PI_UP**2/(72*K*K)
        self.wall={}

    def offsets(self,mode):
        source=[ar(F(m,self.den*self.den)) for _,_,m in self.edges]
        if mode=='mixed':
            return source
        assert mode=='wall'
        s=1-2*WALL_RADIUS
        depth=ar(WALL_RADIUS+self.eta+F(3,5)*(1-s))
        delta=ar(WALL_TARGET-AREA_LOWER)
        out=[]
        capped=[]
        for b,(square,_,_) in zip(source,self.edges,strict=True):
            length=arb(square).sqrt()/self.den
            n=ar(s)*length
            e=upper((delta+(delta**2+2*n*delta*depth).sqrt())/n)
            assert n*ar(e)**2/(2*(ar(e)+depth))>delta
            out.append(b+length*ar(WALL_RADIUS+self.eta+e))
            capped.append(str(e))
        self.wall['protrusion_upper']=capped
        self.wall['fit_budget']=str(WALL_TARGET-AREA_LOWER)
        return out

    def fits(self,cloud,where,pi):
        turn,tx,ty=map(F,where)
        angle=ar(2*turn)*pi
        c,s=angle.cos(),angle.sin()
        points=[]
        for x,y in cloud:
            px,py=rotate((ar(x),ar(y)),c,s)
            points.append((px+ar(tx),py+ar(ty)))
        least=arb(1)
        for p,q in zip(self.material,self.material[1:]+self.material[:1],strict=True):
            normal=(q[1]-p[1],p[0]-q[0])
            margin=arb(dot(normal,p))/self.den
            scale=arb(dot(normal,normal)).sqrt()
            for pt in points:
                value=(margin-dot(normal,pt))/scale
                assert value>ar(F('.0005'))
                least=min(least,value.lower())
        return str(least)


def choose_weights(normals,raw):
    """Exactify a float dual by correcting two entries in a rational null relation."""
    import numpy as np

    d=10**8
    w=[F(int(x),d) for x in np.rint(np.asarray(raw)*d*.9999)]
    sx=-sum(t*a for t,(a,_) in zip(w,normals,strict=True))
    sy=-sum(t*b for t,(_,b) in zip(w,normals,strict=True))
    if sx or sy:
        choices=[(abs(det(normals[i],normals[j])),i,j)
                 for i in range(len(w)) for j in range(i)
                 if w[i]>F(1,100000) and w[j]>F(1,100000)]
        if not choices:
            return None
        for mag,i,j in sorted(choices,reverse=True):
            if not mag:
                break
            a,b=normals[i],normals[j]
            detab=det(a,b)
            r=w[i]+det((sx,sy),b)/detab
            t=w[j]+det(a,(sx,sy))/detab
            if 0<=r<=F(1,2) and 0<=t<=F(1,2):
                w[i],w[j]=r,t
                break
        else:
            return None
    assert all(0<=q<=F(1,2) for q in w)
    assert sum(v*n[0] for v,n in zip(w,normals,strict=True))==0
    assert sum(v*n[1] for v,n in zip(w,normals,strict=True))==0
    return w


def cloud(points,denominator):
    maximum=max(sum((a-b)**2 for a,b in zip(p,q,strict=True)) for p in points for q in points)
    assert maximum<denominator**2
    length=len(points)
    mx,my=(sum(p[j] for p in points) for j in (0,1))
    return ([(F(length*x-mx,length*denominator),F(length*y-my,length*denominator))
             for x,y in points],maximum)


class Problem:
    def __init__(self,host,points,mode,pi,gate):
        self.host,self.points,self.pi,self.gate=host,points,pi,gate
        self.offsets=host.offsets(mode)

    def evaluate(self,left,right,which,rows):
        pi=self.pi
        angle=ar(left+right)*pi
        assert 0<right-left<=F(1,2)
        half=ar(right-left)*pi
        c,s=angle.cos(),angle.sin()
        A,B,C=arb(0),arb(0),arb(0)
        nx,ny=F(0),F(0)
        seen=set()
        for j,i,a,b in rows:
            assert 0<=j<len(self.host.normals) and j not in seen
            seen.add(j)
            assert 0<=i<len(self.points)
            w=F(a,b)
            assert 0<w<=F(1,2)
            n=self.host.normals[j]
            nx+=w*n[0]
            ny+=w*n[1]
            x,y=self.points[i]
            turn=rotate((ar(x),which*ar(y)),c,s)
            uj=self.host.nums[j]
            A+=ar(w)*dot(turn,uj)
            B+=ar(w)*det(turn,uj)
            C+=ar(w)*self.offsets[j]
        assert nx==ny==0
        assert A>0
        return A*half.cos()-abs(B)*half.sin()-C

    def proposal(self,left,right,which):
        import numpy as np

        coords=np.array(self.points,dtype=float)
        coords[:,1]*=which
        n=np.array(self.host.normals,dtype=float)/self.host.den
        c,s=np.cos(float(left+right)*np.pi),np.sin(float(left+right)*np.pi)
        support=n@(coords@np.array([[c,s],[-s,c]])).T
        z=support.max(axis=1)-np.array([float(x.mid()) for x in self.offsets])
        return z,np.argmax(support,axis=1),n

    def suggest(self,left,right,which):
        from scipy.optimize import linprog

        z,which_pt,n=self.proposal(left,right,which)
        opt=linprog(-z,A_eq=n.T,b_eq=[0,0],bounds=(0,.5),method='highs')
        assert opt.success
        if -opt.fun<=float(self.gate):
            return None
        w=choose_weights(self.host.normals,opt.x)
        if not w:
            return None
        return [[j,int(which_pt[j]),x.numerator,x.denominator] for j,x in enumerate(w) if x]

    def make(self,which):
        def interval(i,k):
            left,right=F(i,2**k),F(i+1,2**k)
            candidate=self.suggest(left,right,which)
            if candidate and self.evaluate(left,right,which,candidate)>ar(self.gate):
                yield [i,k,candidate]
                return
            assert k<15,(which,i,k,self.gate)
            yield from interval(2*i,k+1)
            yield from interval(2*i+1,k+1)
        return [row for j in range(64) for row in interval(j,6)]

    def check(self,which,rows):
        previous=F(0)
        floor=arb(1)
        peak=0
        for i,k,part in rows:
            left,right=F(i,2**k),F(i+1,2**k)
            assert previous==left<right<=1
            value=self.evaluate(left,right,which,part)
            assert value>ar(self.gate)
            floor=min(floor,value.lower())
            previous=right
            peak=max(peak,len(part))
        assert previous==1
        return dict(cells=len(rows), max_active_rows=peak, strict_floor=str(floor))


def closed_cell(i,k):
    return F(i,2**k),F(i+1,2**k)


def split_cell(q):
    i,k=q
    return (2*i,k+1),(2*i+1,k+1)


class Joint:
    """Exact shared edge budgets for independently rotated and translated inputs."""
    def __init__(self,host,P,Q,pi):
        assert {(x,-y) for x,y in P} == set(P)
        self.facets=len(host.normals)
        self.problems=(Problem(host,P,'mixed',pi,JOINT_GATE),
                       Problem(host,Q,'mixed',pi,JOINT_GATE))

    def one(self,prob,which):
        def visit(q):
            left,right=closed_cell(*q)
            rows=prob.suggest(left,right,which)
            if rows and prob.evaluate(left,right,which,rows)>ar(JOINT_GATE):
                yield [*q,rows]
            elif q[1] == 9:
                yield [*q,None]
            else:
                for child in split_cell(q):
                    yield from visit(child)
        return [r for j in range(64) for r in visit((j,6))]

    def propose(self,p,q,which):
        import numpy as np
        from scipy.linalg import block_diag
        from scipy.optimize import linprog

        a,ai,n=self.problems[0].proposal(*closed_cell(*p),1)
        b,bi,_=self.problems[1].proposal(*closed_cell(*q),which)
        equality=block_diag(n.T,n.T)
        capacity=np.c_[np.eye(len(n)),np.eye(len(n))]
        opt=linprog(-np.r_[a,b],A_eq=equality,b_eq=[0]*4,A_ub=capacity,
                    b_ub=[.5]*len(n),bounds=(0,.5),method='highs')
        assert opt.success
        if -opt.fun<=float(JOINT_GATE):
            return None
        normals=self.problems[0].host.normals
        x,y=(choose_weights(normals,w) for w in (opt.x[:len(n)],opt.x[len(n):]))
        if x is None or y is None or any(a+b>F(1,2) for a,b in zip(x,y,strict=True)):
            return None
        return [[[j,int(i),a.numerator,a.denominator] for j,(a,i) in enumerate(zip(w,idx,strict=True)) if a]
                for w,idx in ((x,ai),(y,bi))]

    def verify(self,p,q,which,rows):
        weights=[F(0)]*self.facets
        assert len(rows)==2 and any(rows)
        value=arb(0)
        for problem,cell,sign,data in zip(self.problems,(p,q),(1,which),rows,strict=True):
            if data:
                for j,_,a,b in data:
                    assert 0<=j<len(weights)
                    weights[j]+=F(a,b)
                value+=problem.evaluate(*closed_cell(*cell),sign,data)
        assert max(weights)<=F(1,2)
        return value

    def make(self,which):
        sides=[self.one(p,sign) for p,sign in zip(self.problems,(1,which),strict=True)]
        bad=[[row[:2] for row in side if row[2] is None] for side in sides]
        def visit(p,q):
            rows=self.propose(p,q,which)
            if rows is not None and self.verify(p,q,which,rows)>ar(JOINT_GATE):
                return [rows]
            k=0 if p[1]<=q[1] else 1
            assert max(p[1],q[1])<17, (p,q,which)
            return [k,*([visit(child,q) for child in split_cell(p)] if k==0 else
                        [visit(p,child) for child in split_cell(q)])]
        return {'sides':sides, 'weak_products':[visit(p,q) for p in bad[0] for q in bad[1]]}

    def rep(self,which,data):
        unresolved=[]
        floor=arb(1)
        for part,prob,sign in zip(data['sides'],self.problems,(1,which),strict=True):
            before=F(0)
            lacking=[]
            for i,k,rows in part:
                l,r=closed_cell(i,k)
                assert l==before<r<=1
                if rows is None:
                    lacking.append((i,k))
                else:
                    b=prob.evaluate(l,r,sign,rows)
                    assert b>ar(JOINT_GATE)
                    floor=min(floor,b.lower())
                before=r
            assert before==1
            unresolved.append(lacking)
        reached=0
        def visit(p,q,part):
            nonlocal reached,floor
            if len(part)==1:
                b=self.verify(p,q,which,part[0])
                assert b>ar(JOINT_GATE)
                floor=min(floor,b.lower())
                reached+=1
                return
            assert len(part)==3 and part[0] in (0,1)
            k,*children=part
            original=p if k==0 else q
            for child,item in zip(split_cell(original),children,strict=True):
                if k==0:
                    visit(child,q,item)
                else:
                    visit(p,child,item)
        roots=[(p,q) for p in unresolved[0] for q in unresolved[1]]
        assert len(data['weak_products'])==len(roots)
        for (p,q),tree in zip(roots,data['weak_products'],strict=True):
            visit(p,q,tree)
        return dict(one_angle_parts=[len(item) for item in data['sides']],
                    weak_roots=len(roots), joint_terminals=reached, strict_floor=str(floor))


def run(a):
    ctx.prec=160
    sys.path.insert(0,str(a.computation))
    files={k:a.computation/(name+'.json') for k,name in {
        'new_P30':'free-p30-reoptimized-pose','old_P30':'free-p30-placement',
        'escape108':'competing-fits-input','contact102':'five-body-contact-basin',
        'Q11':'multi-host-points','P30':'finite-point-leaf-input',
    }.items()}
    d={k:read(p) for k,p in files.items()}
    poses={k:d[k][:11] for k in ('new_P30','old_P30')}
    poses['escape108']=d['escape108']['zero_gain_five']['params']
    poses['contact102']=d['contact102']['raw_center']
    pi=arb.pi()
    assert pi<ar(PI_UP)
    R=ar(F(703,1000))
    point_area=pi/4+(R*R-ar(F(1,4))).sqrt()/2-(1/(2*R)).acos()/4
    assert point_area>ar(WALL_TARGET)
    hosts={k:Host(list(map(F,p)),pi) for k,p in poses.items()}
    cloud_data={k:cloud(d[k]['points'],d[k][key]) for k,key in (('Q11','point_denominator'),('P30','denominator'))}
    inputs={k:sha(p) for k,p in files.items()}
    config=dict(inputs=inputs, host_inner={k:v.hull_hash for k,v in hosts.items()},
                K=K,edge_count=EDGE_COUNT,wall_radius=str(WALL_RADIUS),
                wall_target=str(WALL_TARGET),area_lower=str(AREA_LOWER),joint_gate=str(JOINT_GATE),
                mixed_gates={f'{h}/{q}':str(v) for (h,q),v in MIXED_GOALS.items()})
    spec=list((k,p,'mixed',g) for (k,p),g in MIXED_GOALS.items())
    spec.extend((k,'Q11','wall',F(0)) for k in WALL_CASES)
    if a.generate and not a.certificate.exists():
        recorded={}
    else:
        data=read(a.certificate)
        assert data['config']==config
        recorded=data['proofs']
    output={}
    for host,key,mode,gate in spec:
        tag=f'{host}-{key}-{mode}'
        problem=Problem(hosts[host],cloud_data[key][0],mode,pi,gate)
        if a.generate and tag not in recorded:
            cases={str(hand):problem.make(hand) for hand in (1,-1)}
            recorded[tag]=cases
            a.certificate.write_text(json.dumps(dict(config=config,proofs=recorded),separators=(',',':'))+'\n')
            print(tag, [len(cases[str(h)]) for h in (1,-1)], file=sys.stderr, flush=True)
        else:
            cases=recorded[tag]
        results={str(hand):problem.check(hand,cases[str(hand)]) for hand in (1,-1)}
        output[tag]=dict(results=results, gate=str(gate), inner_area=str(hosts[host].area),
                         augmented_area_floor=str(hosts[host].area+gate) if mode=='mixed' else None,
                         walls=hosts[host].wall if mode=='wall' else None)
    joint=Joint(hosts['contact102'],cloud_data['P30'][0],cloud_data['Q11'][0],pi)
    for hand in (1,-1):
        tag=f'joint-{hand}'
        if a.generate and tag not in recorded:
            recorded[tag]=joint.make(hand)
            a.certificate.write_text(json.dumps(dict(config=config,proofs=recorded),separators=(',',':'))+'\n')
        result=joint.rep(hand,recorded[tag])
        output[tag]=dict(**result,gate=str(JOINT_GATE),
                         augmented_area_floor=str(hosts['contact102'].area+JOINT_GATE))
        if a.generate:
            print(tag,result,file=sys.stderr,flush=True)
    fits={key:{'P30_pose':p,'support_clearance_lower':hosts[key].fits(cloud_data['P30'][0],p,pi)}
          for key,p in CENTERED_P30_POSES.items()}
    return dict(config=config, output=output, diameters={k:x[1] for k,x in cloud_data.items()}, fits=fits,
                certificate_sha256=sha(a.certificate),source_sha256=sha(Path(__file__)),radius_bound='703/1000')


if __name__=='__main__':
    arguments=ArgumentParser(description=__doc__)
    arguments.add_argument('computation',type=Path)
    arguments.add_argument('certificate',type=Path)
    arguments.add_argument('--generate',action='store_true')
    print(json.dumps(run(arguments.parse_args()),separators=(',',':')))
