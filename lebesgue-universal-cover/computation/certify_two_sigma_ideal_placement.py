"""Certify a common complete placement for both root brackets of a closed ideal family."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
from flint import arb_series
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
UPSTREAM = HERE if (HERE/'certify_sigma_paired_rigidity.py').is_file() else HERE/'upstream'
sys.path.insert(0, str(UPSTREAM))
from certify_sigma_paired_rigidity import rig, g, degree
from certify_eh_upper_chain import front
from certify_eh_sprague_relocation import arc_data, voronoi_signs

h,sol=g.HALF,g.arb
sq3=sol(3).sqrt()
add,sub,mul,rot,dot,cross,neg=g.plus,g.minus,g.scale,g.rot,g.dot,g.cross,g.neg
l0=(sq3/2-h,(sq3/2-h-1)/sq3)
mid=(-h/2,-sq3/4)
z=g.circle_pair(l0,mid,-1)
directions=tuple(g.unit(i*g.PI/3) for i in range(6))
def qa(v):
    return ((-v[0]+sq3*v[1])/2,(sq3*v[0]+v[1])/2)

def lerp(a,b,s):
    return add(a,mul(s,sub(b,a)))

def xy(sigma):
    g.SIGMA=degree(sigma)
    n,m,c3,b3,_,_,_,t_eh,_=g.setup()
    e3=g.line_pair(n[4],m[3])
    c=(-e3[0],e3[1])
    v=sub(z,c)
    vv=dot(v,v)
    ta=add(c,mul(1/vv,sub(v,mul((vv-1).sqrt(),rot(v)))))
    D=sub(ta,z)
    return m,c3,b3,c,ta,D,t_eh

def anchors(D,w):
    x=add(z,mul(w,D))
    def inverse(u,e):
        off=dot(x,u)-h
        return sub(sub(x,mul(off,u)),mul((1-off**2).sqrt(),e))
    ll=inverse(directions[5],(sq3/2,h))
    nn=inverse(directions[4],(-sq3/2,h))
    t=ll[0]+h
    b0=(t+sq3*(1-t*t).sqrt())/2-h
    b=(b0,(b0-1)/sq3)
    f,N=qa(ll),qa(nn)
    q=g.circle_pair(b,N,-1)
    return [b,sub(b,directions[5]),f,add(f,directions[0]),N,sub(N,directions[4]),q,qa(x)]

def ea(a,b):
    return g.exact(a).union(g.exact(b))

def edges(data):
    n,m,c3,b3,_,_,_,t_eh,_=g.setup()
    p=tuple(g.arb(x) for x in ('.429','.009'))
    _,_,_,c,ta,D,_=data
    v=neg(mul(1/g.norm(D),D))
    ca,sa=g.unit(p[1])
    pivot_dir=add(mul(ca,v),mul(sa,rot(v)))
    pivot=add(z,mul(p[0],D))
    def qp(w):
        return anchors(D,w)[6]
    def cap_excess(w):
        return -cross(pivot_dir,sub(qp(g.exact(w)),pivot))
    a,b=F('.20'),F('.40')
    assert cap_excess(a)>0 and cap_excess(b)<0
    for _ in range(60):
        m0=(a+b)/2
        value=cap_excess(m0)
        if value>0:
            a=m0
        else:
            assert value<0
            b=m0
    root=ea(a,b)
    # Endpoint of the new AH chord and genuine retained paired arc.
    start=qp(root)
    qps=[qp(root+(1-root)*g.exact(F(j,16))) for j in range(17)]
    e3=g.line_pair(n[4],m[3])
    y=(sol(0),c[1]+(1-c[0]**2).sqrt())
    eh=front(json.loads((UPSTREAM/'eh-two-line-cover-certificate.json').read_bytes()),
             json.loads((UPSTREAM/'eh-upper-chain-parameters.json').read_bytes()))
    hn,ho=eh['normals'],eh['offsets']
    pairs=list(zip(eh['active'],eh['active'][1:]))
    chain=[eh['v']]+[g.line_pair(hn[i],hn[j],ho[i],ho[j]) for i,j in pairs]+[eh['t']]
    outline=[(y,ta,'YTa inner chord'),(ta,pivot,'retained ZTa line'),
             (pivot,start,'new AH chord')]
    outline.extend((u,v,'paired AH inner chord') for u,v in zip(qps,qps[1:]))
    outline.extend((x,y,f'EH line {i+1}') for i,(x,y) in enumerate(zip(chain,chain[1:])))
    line_rows=[]
    for a,b,name in outline:
        t=sub(b,a)
        u=neg(mul(1/g.norm(t),rot(t)))
        line_rows.append((name,u,dot(u,a)))
    line_rows.extend((f'm{j}',m[j],h) for j in (3,5))
    # Inherited five arcs are included by the P-placement lemma. Recheck the
    # witness endpoint inequalities at this sigma rather than reuse fixed numbers.
    vs=voronoi_signs(arc_data())
    inherited=arc_data()
    gcenter=inherited[3][1]; acenter=inherited[4][1]
    qp1=qps[-1]
    q0=add(acenter,n[1]); p0=add(gcenter,n[2])
    major=[(acenter,q0,y,inherited[4][4],False),(gcenter,qp1,p0,inherited[3][4],True),
           (c3,add(c3,n[3]),chain[0],g.setup()[8],False),
           (b3,t_eh,add(b3,m[3]),g.setup()[8],True)]
    # Check each retained span begins no earlier and ends no later than the
    # R90 Sprague arc, with a common open semicircle for these radial tests.
    for center,first,last,strict,short_left in major:
        f,l,s=sub(first,center),sub(last,center),sub(strict,center)
        assert cross(f,l)>0 and dot(f,l)>0
        if short_left:
            assert cross(s,f)>0 and dot(s,f)>0
        else:
            assert cross(l,s)>0 and dot(l,s)>0
    assert len(eh['active'])==11
    outline_proof=certify_outline(n,m,c3,b3,c,e3,y,ta,pivot,qps,chain,D,root)
    return line_rows, vs, (a,b), outline_proof

def certify_outline(n,m,c3,b3,c,e3,y,ta,pivot,qps,chain,D,root):
    b1=g.line_pair(n[0],n[1])
    f1=g.line_pair(n[2],n[3])
    d1=g.line_pair(n[4],n[5])
    f3=g.line_pair(n[3],m[2])
    center_g=g.circle_line(f3,n[5],1)
    Q=add(e3,n[1]); P=add(center_g,n[2])
    H=add(c3,n[3]); I=add(b3,m[3]); K=add(f3,m[5])
    def derivative(w):
        q=anchors(D,arb_series([w,1,0],3))[6]
        return (q[0][1],q[1][1]),(2*q[0][2],2*q[1][2])
    curv=[]
    cells=131072
    for j in range(cells):
        v,dv=derivative(ea(F(j,cells),F(j+1,cells)))
        assert v[0]<0 and v[1]<0,(j,v)
        curv.append(cross(v,dv).lower())
    assert min(curv)>0
    v0,_=derivative(root)
    v1,_=derivative(sol(1))
    def circ(a,b,o):
        for p in (a,b):
            assert (dot(sub(p,o),sub(p,o))-1).contains(0)
        return rot(sub(a,o)),rot(sub(b,o)),True
    def segment(a,b):
        d=sub(b,a)
        return d,d,False
    outline=[
        ('n1',segment(b1,Q)),('QY',circ(Q,y,e3)),('YTa',circ(y,ta,c)),
        ('TaW',segment(ta,pivot)),('new AH',segment(pivot,qps[0])),
        ('paired AH',(v0,v1,True)),('paired G',circ(qps[-1],P,center_g)),
        ('n2',segment(P,f1)),('n3',segment(f1,H)),('C3',circ(H,chain[0],c3)),
    ]
    outline.extend((f'EH {i+1}',segment(a,b)) for i,(a,b) in enumerate(zip(chain,chain[1:])))
    outline.extend([
        ('B3',circ(chain[-1],I,b3)),('m3',segment(I,e3)),
        ('n4',segment(e3,d1)),('n5',segment(d1,center_g)),
        ('F3',circ(center_g,K,f3)),('m5',segment(K,c3)),('n0',segment(c3,b1)),
    ])
    tangent={'QY','TaW','paired G','n2','C3','m3','m5'}
    total=sol(0); joint=[]
    def turn(a,b,label,zero=False):
        y0,x0=cross(a,b),dot(a,b)
        assert x0>0,(label,x0)
        if zero:
            assert y0.contains(0),(label,y0)
            return sol(0)
        assert y0>0,(label,y0)
        return sol.atan2(y0,x0)
    for j,(label,(a,b,curved)) in enumerate(outline):
        if curved:
            total+=turn(a,b,label+' internal')
        at_join=turn(outline[j-1][1][1],a,label,label in tangent)
        total+=at_join
        if label not in tangent:
            joint.append(at_join.lower())
    assert total.contains(2*g.PI) and total>g.PI and total<3*g.PI
    e_l=(sq3/2,h)
    # At w=1 the inverse DC root is L=c. The condition (Ta-c)·D=0,
    # with (Ta-c)·e_l ≠0, then proves L'=0 and b'=0 analytically.
    # The unit q--b condition gives exact tangent matching to the G circle.
    assert (dot(sub(ta,c),D)).contains(0) and dot(sub(ta,c),e_l)>0
    return {'curvature_cells':len(curv),'curvature_cross_lower':str(min(curv)),
            'positive_join_turn_lower_radians':str(min(joint)),
            'total_turn_radians':str(total),
            'exact_tangent_joins': sorted(tangent)}


def lens(points):
    vp=[points[i] for i in (0,2,7,4,6)]
    possible=[]
    def admiss(v,skip=()):
        return all(not dot(sub(v,c),sub(v,c))>1 for j,c in enumerate(vp) if j not in skip)
    for i,j in combinations(range(5),2):
        for s in (-1,1):
            w=g.circle_pair(vp[i],vp[j],s)
            if admiss(w,(i,j)):
                possible.append(w)
    wedges=[]
    for i,p in enumerate(vp):
        a=sub(p,vp[i-1]); b=sub(p,vp[(i+1)%5])
        if cross(a,b)<0:
            a,b=b,a
        assert cross(a,b)>0
        wedges.append((p,a,b))
    def support(u):
        low=max(dot(p,u).lower() for p in points)
        for p,a,b in wedges:
            if cross(a,neg(u))>=0 and cross(neg(u),b)>=0:
                low=max(low,(1+dot(p,u)).lower())
        high=max(dot(v,u).upper() for v in possible)
        for i,p in enumerate(vp):
            if admiss(add(p,u),(i,)):
                high=max(high,(1+dot(p,u)).upper())
        assert low<=high
        return low.union(high)
    return support

def validate_cell(a,b,theta0,theta1,D,rows):
    p=anchors(D,ea(a,b))
    H=lens(p)
    c,s=g.unit(degree(theta0-300).union(degree(theta1-300)))
    def pulled(u):
        u=(u[0],-u[1])
        return add(mul(c,u),mul(s,rot(u)))
    right0=h-H(pulled(directions[0]))
    right2=h-H(pulled(directions[2]))
    shift=g.line_pair(directions[0],directions[2],right0,right2)
    vals=[]
    for name,u,bd in rows:
        margin=bd-H(pulled(u))-dot(u,shift)
        vals.append((margin.lower(),name))
    return min(vals)

def bracket(sigma,a,b):
    vals=[a,b]
    return rig(sigma,sigma,*PARAM,probe_angles=vals)

PARAM=(F('.42849'),F('.42851'))
SIGMAS={
 '1.54':[('53.47','53.50'),('54.18','54.21')],
 '1.5497':[('53.69','53.73'),('53.99','54.03')],
}

def run(sigma_label,brackets):
    sigma=F(sigma_label)
    w0,w1=PARAM
    data=xy(sigma)
    D=data[5]
    bounds,arcs,crossing,outline=edges(data)
    out=[]
    for i,(left,right) in enumerate(brackets):
        lut,rr=F(left),F(right)
        angles=[lut,rr+60] if i==0 else [lut+60,rr]
        lower=[]
        for j in range(32):
            a,b=[w0+(w1-w0)*F(z,32) for z in (j,j+1)]
            lower.extend(rig(sigma,sigma,a,b,probe_angles=angles))
        min_sign=min(lower)
        assert min_sign>0,(sigma_label,i,str(min_sign))
        scores=[]
        for j in range(16):
            a,b=[w0+(w1-w0)*F(z,16) for z in (j,j+1)]
            scores.append(validate_cell(a,b,lut,rr,D,bounds))
        min_slack,wall=min(scores)
        assert min_slack>0,(sigma_label,i,str(min_slack),wall,[(str(x),y) for x,y in scores])
        out.append({'bracket_degrees':[left,right],
                    'minimum_signed_endpoint_lower':str(min_sign),
                    'minimum_inward_normal_slack_lower':str(min_slack),
                    'weakest_row':wall,'root_endpoint_parameter_cells':32,
                    'same_pose_parameter_cells':16,'support_rows_per_box':len(bounds)})
    return {'sigma_degrees':sigma_label,'branches':out,'Sprague_voronoi':arcs,'convex_outline':outline,
            'AH_crossing_rational_bracket':[str(q) for q in crossing]}

if __name__=='__main__':
    data={s:run(s,brackets) for s,brackets in SIGMAS.items()}
    inputs=('certify_sigma_paired_rigidity.py','certify_eh_upper_chain.py',
            'certify_eh_sprague_relocation.py','eh-two-line-cover-certificate.json',
            'eh-upper-chain-parameters.json')
    hashes={f:hashlib.sha256((UPSTREAM/f).read_bytes()).hexdigest() for f in inputs}
    payload={'complete':True,'input_sha256':hashes,'w_closed_interval_rational':list(map(str,PARAM)),
             'new_AH_pivot':'.429','new_AH_angle_radians':'.009','EH_facets':11,
             'orthogonal_part':'F_y R_(300deg-phi)', 'interval_precision_bits':g.ctx.prec,
             'certificates':data,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(sys.argv[1]).write_text(json.dumps(payload,indent=2)+'\n')
