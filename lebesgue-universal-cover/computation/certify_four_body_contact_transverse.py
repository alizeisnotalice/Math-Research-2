"""Supporting fan velocities in six directions transverse to the four-body manifold."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import flint
from flint import arb,arb_mat,ctx
from certify_four_body_contact_stationary import BOX,add,chain,det,dot,interval,left,normal,q,r,scale,sub

EVENT_ORDER=('A','Aopp','B','Bopp','C','D')
SIMPLEX=(83,261,267,326,339,343,791)


def plus(a,b):
    return [x+y for x,y in zip(a,b)]


def minus(a,b):
    return [x-y for x,y in zip(a,b)]


def rigid_velocities(owner,point,t):
    v=[(0,0) for _ in range(6)]
    if owner==5:
        v[:3]=[left(sub(point,t[3:5])),(1,0),(0,1)]
    elif owner==7:
        v[3:]=[left(sub(point,t[6:8])),(1,0),(0,1)]
    return v


def connector(i,P,j,Q,t):
    W=map(add,rigid_velocities(i,P,t),rigid_velocities(j,Q,t))
    return [det(v,sub(Q,P))/2 for v in W]


def face_orders(geometry,ns):
    p=geometry
    a,b,c,d=(ns[k] for k in ('A','B','C','D'))
    M={'A':[(0,scale(r,a)),(3,p[3][0]),(5,add(p[5][3],a)),(7,p[7][0])],
       'B':[(0,scale(r,b)),(3,p[3][1]),(5,p[5][1]),(7,p[7][1])],
       'C':[(0,scale(r,c)),(5,add(p[5][2],c)),(7,add(p[7][2],c))],
       'D':[(0,scale(r,d)),(3,p[3][0]),(7,add(p[7][2],d))]}
    for label,points in M.items():
        n=ns[label];tangent=left(n)
        points.sort(key=lambda row:float(dot(tangent,row[1]).mid()))
        assert all(dot(tangent,p2)>dot(tangent,p1) for (_,p1),(_,p2) in zip(points,points[1:]))
    for label in ('A','B'):
        M[label+'opp']=[(i,sub(v,ns[label])) for i,v in reversed(M[label])]
    return M


def gradients(a,b):
    pieces,t,ns,p=chain(a,b)
    orders=face_orders(p,ns)
    base=[arb(0)]*6; increments={}
    start=arb(0)
    for (tag,stop,i,rad,center),(_,_,j,rad2,center2) in zip(pieces,pieces[1:]+pieces[:1]):
        n0,n1=normal(start),normal(stop)
        delta=sub(n1,n0)
        base=plus(base,[rad*det(v,delta) for v in rigid_velocities(i,center,t)])
        P=add(center,scale(rad,n1));Q=add(center2,scale(rad2,n1))
        g=connector(i,P,j,Q,t)
        base=plus(base,g)
        if tag in orders:
            seq=orders[tag]
            assert seq[0][0]==i and seq[-1][0]==j
            increments[tag]=[]
            for select in product((0,1),repeat=len(seq)-2):
                items=[seq[0],*(item for yes,item in zip(select,seq[1:-1]) if yes),seq[-1]]
                inc=[arb(0)]*6
                for old,new in zip(items,items[1:]):
                    inc=plus(inc,connector(*old,*new,t))
                increments[tag].append(minus(inc,g))
        start=stop
    all_gradients=[]
    for keys in product(*(range(len(increments[key])) for key in EVENT_ORDER)):
        total=list(base)
        for i,k in zip(keys,EVENT_ORDER):total=plus(total,increments[k][i])
        all_gradients.append(total)
    return all_gradients,{k:[i for i,_ in orders[k]] for k in EVENT_ORDER}


def locate():
    import numpy as np
    from scipy.optimize import linprog
    g,_=gradients(q('.512504445326185003543919'),q('1.577473015380594982391575'))
    G=np.array([[float(v.mid()) for v in row] for row in g],dtype=float)
    A=np.vstack((np.ones(len(G)),G.T)); b=np.zeros(7); b[0]=1
    rng=np.random.default_rng(20261118)
    best=None
    for _ in range(64):
        res=linprog(rng.normal(size=len(G)),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
        if not res.success: continue
        cols=np.where(res.x>1e-9)[0]
        if len(cols)!=7: continue
        w=np.linalg.solve(A[:,cols],b)
        if best is None or min(w)>best[0]:best=(min(w),cols.tolist(),w.tolist(),np.linalg.det(A[:,cols]))
    print(best)


def certify():
    nodes=tuple(interval(*row) for row in BOX)
    _,_,ns,p=chain(*nodes)
    assert det(ns['A'],ns['B'])>0
    assert all(abs(dot(left(ns[key]),p[3][index]))>q('.1') for key,index in [('A',0),('B',1)])
    g,order=gradients(*nodes)
    A=arb_mat([[1]*7,*[[g[i][j] for i in SIMPLEX] for j in range(6)]])
    det_A=A.det();assert det_A>0 or det_A<0
    rhs=arb_mat([[1]]+[[0]]*6)
    w=A.solve(rhs)
    for i in range(7): assert w[i,0]>q('.001')
    B=A.inv()
    # Rows give barycentric coordinate variations on the entire centred cube.
    cube=q('.0001')
    for i in range(7):
        assert w[i,0]>cube*sum(abs(B[i,j]) for j in range(1,7))
    path=Path(__file__)
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(path.read_bytes()).hexdigest(),'geometry_script_sha256':sha256(path.with_name('certify_four_body_contact_stationary.py').read_bytes()).hexdigest(),
        'stationary_angle_box':BOX,'sorted_contact_order':order,'number_of_available_fans':len(g),
        'simplex_fan_indices':SIMPLEX,'simplex_determinant':str(det_A),'origin_barycentric_weights':[str(w[i,0]) for i in range(7)],
        'support_gradient_cube_radius':'1/10000','scope':'Uniform transverse first-order support in the verified two-angle box.'},indent=2))

if __name__=='__main__':
    parser=ArgumentParser();parser.add_argument('--locate',action='store_true');args=parser.parse_args()
    locate() if args.locate else certify()
