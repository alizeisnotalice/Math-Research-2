"""Full-variable necessary-point LP for a frozen spline kernel family.

prepare: run in existing python-flint environment; solve: run in SciPy env.
No LP result from this file is a new weak-type bound.
"""
import argparse, hashlib, importlib.util, json, os, platform, time, struct
from fractions import Fraction as F
from pathlib import Path


def orbit(i,j):
    return sorted({(s*a,t*b) for a,b in ((i,j),(j,i)) for s in (-1,1) for t in (-1,1)})


def load_core(root):
    spec=importlib.util.spec_from_file_location('floor_core',root/'hl_core.py')
    core=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    core.ctx.prec=192
    return core


def family_spec(core,data):
    """Only frozen family parameters; no candidate coefficient values enter a proof."""
    sp=object.__new__(core.Spec)
    sp.alpha,sp.R,sp.h=(F(data[k]) for k in ('alpha','R','h'))
    assert 1<sp.alpha<2 and sp.R>0 and sp.h>0
    sp.acls={tuple(v):F(1) for v in data['classes']}
    sp.maxidx=max(i for i,j in sp.acls)
    sp.a_=core.A(sp.alpha)
    a=sp.a_
    sp.e=3-a
    sp.kappa=-2/(a*(1-a)*(2-a)*(3-a))
    sp.hA=core.A(sp.h)
    sp.sc=core.pw(sp.hA,-a)*sp.kappa
    sp.RA=core.A(sp.R)
    sp.Ra=core.pw(sp.RA,-a)
    sp.Ca=-8*(1-a).gamma()*(2*a).gamma()*(core.arb.pi()*a/2).cos()**2/(a+1).gamma()
    return sp


def prepare(args):
    import flint
    start=time.time()
    root=Path(args.root)
    core=load_core(root)
    if args.certificate:
        saved=json.loads(Path(args.certificate).read_text())
        sp=family_spec(core,saved['variables'])
        candidate_sha=saved['candidate_sha256']
    else:
        sp=core.Spec(root/'candidate.json')
        candidate_sha=hashlib.sha256((root/'candidate.json').read_bytes()).hexdigest()
    classes=sorted(sp.acls)
    end=sp.maxidx+2
    # C1/C2 at full knot grid, excluding the singular origin.
    kvals=[(F(i)*sp.h,F(j)*sp.h) for i in range(end+1) for j in range(i+1) if i+j]
    kvals += [(sp.h/2,sp.h/2)]
    if args.step==1:
        kvals += [((F(i)+F(1,2))*sp.h,(F(j)+F(1,2))*sp.h)
                  for i in range(end) for j in range(i+1)]
        kvals=sorted(set(kvals))
    # C4 exact nonaxial points. No variable is dropped at this coarse first stage.
    avals=[((F(i)+F(1,2))*sp.h,(F(j)+F(1,2))*sp.h)
           for i in range(0,end,args.step) for j in range(0,i+1,args.step)]
    avals += [(F(35,64),F(5,256))]
    if args.step==1:
        avals += [(F(i)*sp.h,F(j)*sp.h) for i in range(1,end+1) for j in range(1,i+1)]
        avals += [((F(i)+F(1,2))*sp.h,sp.h/64) for i in range(end)]
    avals=sorted(set(avals))
    rows=[]
    for u,v in kvals:
        value=sp.Psi(core.A(u+v))
        rows.append({'kind':'K','point':[str(u),str(v)],'rhs':int(u+v<=1),
                     'singular_mid':float(value.mid()),'singular_ball':str(value)})
    for u,v in avals:
        assert u>0 and v>0 and u+v!=sp.R
        value=sp.APsi_point(u,v)
        assert value.is_finite()
        rows.append({'kind':'A','point':[str(u),str(v)],'rhs':0,
                     'singular_mid':float(value.mid()),'singular_ball':str(value)})
    m0=sp.a_*core.pw(sp.RA,2-sp.a_)/(2-sp.a_)
    data={'scope':'necessary finite constraints; all 1181 frozen spline classes plus c retained',
          'candidate_sha256':candidate_sha,
          'core_sha256':hashlib.sha256((root/'hl_core.py').read_bytes()).hexdigest(),
          'alpha':str(sp.alpha),'R':str(sp.R),'h':str(sp.h),
          'classes':[list(x) for x in classes], 'variable_count':len(classes)+1,
          'rows':rows,'mass_singular_mid':float(m0.mid()),'mass_singular_ball':str(m0),
          'K_rows':len(kvals),'A_rows':len(avals),'generator_step_in_knots':args.step,
          'python':platform.python_version(),'python_flint':flint.__version__,
          'precision_bits':192,'prepare_seconds':time.time()-start}
    Path(args.out).write_text(json.dumps(data,indent=2))
    print(json.dumps({k:v for k,v in data.items() if k not in ('rows','classes')},indent=2))


def make_matrix(data):
    import numpy as np
    alpha=float(F(data['alpha'])); h=float(F(data['h']))
    rows=data['rows']; classes=data['classes']
    coords=np.array([[float(F(v))/h for v in row['point']] for row in rows])
    maxidx=max(i for i,j in classes)
    indices=np.arange(-maxidx,maxidx+1)
    def B(x):
        ax=np.abs(x)
        return np.where(ax<1,(4-6*ax**2+3*ax**3)/6,np.where(ax<2,(2-ax)**3/6,0.))
    def LB(x):
        kappa=-2/(alpha*(1-alpha)*(2-alpha)*(3-alpha))
        out=np.zeros_like(x)
        for k,w in zip(range(-2,3),(1,-4,6,-4,1)):
            out+=w*np.abs(x-k)**(3-alpha)
        return h**(-alpha)*kappa*out
    x=coords[:,0,None]-indices[None,:]
    y=coords[:,1,None]-indices[None,:]
    bx,by=B(x),B(y)
    lx,ly=LB(x),LB(y)
    is_a=np.array([row['kind']=='A' for row in rows])
    G=np.empty((len(rows),len(classes)+1))
    G[:,0]=[row['singular_mid'] for row in rows]
    mass=np.empty(len(classes)+1); mass[0]=data['mass_singular_mid']
    for col,(i,j) in enumerate(classes,1):
        orb=orbit(i,j)
        kval=np.zeros(len(rows)); aval=np.zeros(len(rows))
        for a,b in orb:
            ia,ib=a+maxidx,b+maxidx
            kval+=bx[:,ia]*by[:,ib]
            aval+=lx[:,ia]*by[:,ib]+bx[:,ia]*ly[:,ib]
        G[:,col]=np.where(is_a,aval,kval)
        mass[col]=h*h*len(orb)/2
    return G,np.array([row['rhs'] for row in rows],dtype=float),mass


def solve_lp(args):
    import numpy as np, scipy
    from scipy.optimize import linprog
    start=time.time()
    data=json.loads(Path(args.input).read_text())
    G,b,mass=make_matrix(data)
    # c>=0; all spline coefficients unrestricted. No zero-feature pricing omission.
    # Normalize rows to improve conditioning; map dual back to original rows.
    scale=np.maximum(1.,np.max(np.abs(G),axis=1))
    res=linprog(mass,A_ub=-G/scale[:,None],b_ub=-b/scale,
                bounds=[(0,None)]+[(None,None)]*(len(mass)-1),method='highs',
                options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,
                         'threads':1})
    result={'status':int(res.status),'message':res.message,'success':bool(res.success),
            'variable_count':len(mass),'row_count':len(b),'all_frozen_variables_retained':True,
            'input_sha256':hashlib.sha256(Path(args.input).read_bytes()).hexdigest(),
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
            'solve_seconds':time.time()-start,'rigorous_certificate':False,
            'mathematical_role':'candidate lower floor for this frozen kernel family; never a C2 bound'}
    if res.success:
        dual=-res.ineqlin.marginals/scale
        positive=np.flatnonzero(dual>0)
        residual=G.T@dual-mass
        result.update({'float_objective':float(res.fun),'float_dual_objective':float(b@dual),
                       'max_abs_column_residual':float(np.max(np.abs(residual))),
                       'positive_dual_count':len(positive),
                       'positive_dual':[{'row':int(k),'weight_float':float(dual[k])} for k in positive],
                       'primal_coefficients_float':res.x.tolist(),
                       'c_lower_bound_dual':float(res.lower.marginals[0]),
                       'minimum_primal_point_slack':float(np.min(G@res.x-b))})
        # Regenerable floating matrix kept only in caller's data-temp directory.
        if args.matrix:
            np.savez_compressed(args.matrix,G=G,b=b,mass=mass,dual=dual,primal=res.x)
    Path(args.out).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ('positive_dual','primal_coefficients_float')},indent=2))


def certify(args):
    """Rebuild the exact computable-real active matrix, then certified Arb solve."""
    import flint
    from flint import arb,arb_mat
    start=time.time()
    data=json.loads(Path(args.input).read_text())
    if 'active_constraints' in data:
        data={**data,**data['variables']}
        active=data['active_constraints']
        lp={'input_sha256':data['source_input_sha256']}
    else:
        lp=json.loads(Path(args.lp).read_text())
        assert lp['success'] and lp['all_frozen_variables_retained']
        assert lp['input_sha256']==hashlib.sha256(Path(args.input).read_bytes()).hexdigest()
        active=[data['rows'][v['row']] for v in lp['positive_dual']]
        if lp['c_lower_bound_dual']>1e-10:
            active.append({'kind':'c_nonnegative','rhs':0})
    root=Path(args.root)
    core=load_core(root)
    assert hashlib.sha256((root/'hl_core.py').read_bytes()).hexdigest()==data['core_sha256']
    sp=family_spec(core,data)
    classes=[tuple(v) for v in data['classes']]
    assert len(classes)==1181 and len(set(classes))==len(classes)
    assert all(0<=j<=i for i,j in classes)
    assert classes==sorted(sp.acls)
    assert (sp.alpha,sp.R,sp.h)==tuple(F(data[k]) for k in ('alpha','R','h'))
    n=len(classes)+1
    assert len(active)==n, (len(active),n)
    matrix=arb_mat(n,n)
    mass=arb_mat(n,1)
    mass[0,0]=sp.a_*core.pw(sp.RA,2-sp.a_)/(2-sp.a_)
    for k,(i,j) in enumerate(classes,1):
        mass[k,0]=core.A(sp.h**2*len(orbit(i,j))/2)
    def B(x):
        x=abs(x)
        if x>=2:return F(0)
        if x>=1:return (2-x)**3/6
        return (4-6*x*x+3*x**3)/6
    power_cache={}
    lb_cache={}
    orbits=[orbit(i,j) for i,j in classes]
    all_indices=range(-sp.maxidx,sp.maxidx+1)
    def LB(x):
        x=abs(x)
        if x not in lb_cache:
            value=arb(0)
            for k,w in zip(range(-2,3),(1,-4,6,-4,1)):
                q=abs(x-k)
                if q not in power_cache:power_cache[q]=core.apow(core.A(q),sp.e)
                value+=w*power_cache[q]
            lb_cache[x]=sp.sc*value
        return lb_cache[x]
    for col,row in enumerate(active):
        if row['kind']=='c_nonnegative':
            assert row['rhs']==0
            matrix[0,col]=1
            continue
        u,v=map(F,row['point'])
        assert 0<=v<=u and u+v>0
        if row['kind']=='K':
            assert row['rhs']==int(u+v<=1)
            matrix[0,col]=sp.Psi(core.A(u+v))
        else:
            assert row['kind']=='A' and row['rhs']==0
            assert v>0 and u+v!=sp.R
            matrix[0,col]=sp.APsi_point(u,v)
        X,Y=u/sp.h,v/sp.h
        bxs={a:B(X-a) for a in all_indices}
        bys={b:B(Y-b) for b in all_indices}
        if row['kind']=='A':
            lxs={a:LB(X-a) for a in all_indices}
            lys={b:LB(Y-b) for b in all_indices}
        for k,orb in enumerate(orbits,1):
            val=arb(0)
            for a,b in orb:
                bx,by=bxs[a],bys[b]
                if row['kind']=='K':
                    if bx and by:val+=core.A(bx*by)
                else:
                    if by:val+=lxs[a]*core.A(by)
                    if bx:val+=core.A(bx)*lys[b]
            if not val.is_zero():matrix[k,col]=val
        if col%200==0:
            print(json.dumps({'stage':'matrix','columns_done':col+1,'total':n}),flush=True)
    matrix_seconds=time.time()-start
    if args.neumann:
        saved=json.loads(Path(args.certificate or args.input).read_text())
        assert saved['variables']['classes']==data['classes']
        assert len(saved['active_constraints'])==n
        y0=[]
        for row,previous in zip(active,saved['active_constraints']):
            assert row['kind']==previous['kind'] and row.get('point')==previous['point'] and row['rhs']==previous['rhs']
            y0.append(F(previous['y0']) if 'y0' in previous else (F(previous['dual']['lo'])+F(previous['dual']['hi']))/2)
        Y0=arb_mat(n,1,[core.A(v) for v in y0])
        raw=Path(args.preconditioner).read_bytes()
        assert len(raw)==n*n*8
        R=arb_mat(n,n,[core.A(F(v[0])) for v in struct.iter_unpack('<d',raw)])
        print(json.dumps({'stage':'neumann_matrix_product','n':n}),flush=True)
        E=R*matrix
        for i in range(n):E[i,i]-=1
        row_norms=[sum(abs(E[i,j]) for j in range(n)).upper() for i in range(n)]
        q=max(row_norms)
        assert q<1
        residual=mass-matrix*Y0
        correction=R*residual
        eta=max(abs(correction[i,0]).upper() for i in range(n))
        delta=eta/(1-q)
        assert all(core.A(v)>delta for v in y0)
        objective0=sum(row['rhs']*v for row,v in zip(active,y0))
        lower=core.A(objective0)-sum(row['rhs'] for row in active)*delta
        def endpoint(x):
            return {'lo':str(x.lower().fmpq()),'hi':str(x.upper().fmpq()),'display':str(x)}
        result={'certified':True,'method':'Neumann residual proof using matrix multiplication only; no certified inverse or solve',
                'scope':'fixed kernel family cost floor only; not a C2 lower or upper bound',
                'variables':saved['variables'],
                'active_constraints':[{'kind':row['kind'],'point':row.get('point'),'rhs':row['rhs'],'y0':str(v)} for row,v in zip(active,y0)],
                'all_1182_columns_checked':True,'free_spline_columns_require_equalities':True,
                'preconditioner_format':'little-endian float64 row-major; entries interpreted as exact dyadic rationals',
                'preconditioner_sha256':hashlib.sha256(raw).hexdigest(),
                'q_upper':endpoint(q),'eta_upper':endpoint(eta),'uniform_solution_error':endpoint(delta),
                'minimum_y0':str(min(y0)),'objective_at_y0':str(objective0),
                'certified_objective_lower':endpoint(lower),
                'precision_bits':192,'python':platform.python_version(),'python_flint':flint.__version__,
                'matrix_build_seconds':matrix_seconds,'total_seconds':time.time()-start,
                'source_input_sha256':saved['source_input_sha256'],'candidate_sha256':saved['candidate_sha256'],'core_sha256':saved['core_sha256']}
        Path(args.out).write_text(json.dumps(result,indent=2))
        print(json.dumps({k:v for k,v in result.items() if k not in ('variables','active_constraints','objective_at_y0')},indent=2))
        return
    print(json.dumps({'stage':'certified_solve','n':n,'matrix_seconds':matrix_seconds}),flush=True)
    y=matrix.solve(mass,algorithm='precond')
    assert all(y[i,0].is_finite() for i in range(n))
    negative=[i for i in range(n) if not (y[i,0]>=0)]
    objective=sum(row['rhs']*y[i,0] for i,row in enumerate(active))
    def ball(v):
        return {'lo':str(v.lower().fmpq()),'hi':str(v.upper().fmpq()),'display':str(v)}
    cert={'scope':'cost floor for frozen 1182-variable kernel family only; not C2 lower/upper bound',
          'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':192,
          'algorithm':'arb_mat.solve(precond), rigorous interval solution of full active transpose matrix',
          'variables':{'alpha':data['alpha'],'R':data['R'],'h':data['h'],'classes':data['classes']},
          'source_input_sha256':lp['input_sha256'],'candidate_sha256':data['candidate_sha256'],
          'core_sha256':data['core_sha256'],'active_constraints':[
              {'kind':row['kind'],'point':row.get('point'),'rhs':row['rhs'],'dual':ball(y[i,0])}
              for i,row in enumerate(active)],
          'objective':ball(objective),'all_dual_nonnegative':not negative,
          'unproved_nonnegative_indices':negative,'certified':not negative,
          'matrix_build_seconds':matrix_seconds,'total_seconds':time.time()-start}
    if not negative:
        lo=F(str(objective.lower().fmpq()))
        cert['certified_family_floor']=str(F((lo*1000000).__floor__(),1000000))
    Path(args.out).write_text(json.dumps(cert,indent=2))
    print(json.dumps({k:v for k,v in cert.items() if k not in ('active_constraints','variables')},indent=2))


def preconditioner(args):
    import numpy as np
    from scipy.linalg import inv
    start=time.time()
    data=json.loads(Path(args.input).read_text())
    lp=json.loads(Path(args.lp).read_text())
    G,b,m=make_matrix(data)
    ids=[row['row'] for row in lp['positive_dual']]
    A=G[ids,:].T
    if lp['c_lower_bound_dual']>1e-10:
        extra=np.zeros((len(m),1)); extra[0,0]=1
        A=np.hstack([A,extra])
    assert A.shape==(len(m),len(m))
    R=inv(A,check_finite=True)
    raw=R.astype('<f8').tobytes(order='C')
    Path(args.out).write_bytes(raw)
    print(json.dumps({'stage':'untrusted_float_preconditioner','dimension':len(m),
                      'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),
                      'seconds':time.time()-start,'format':'little-endian float64 row-major'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['prepare','solve','certify','preconditioner'])
    parser.add_argument('--root')
    parser.add_argument('--input')
    parser.add_argument('--out',required=True)
    parser.add_argument('--step',type=int,default=2)
    parser.add_argument('--matrix')
    parser.add_argument('--lp')
    parser.add_argument('--neumann',action='store_true')
    parser.add_argument('--certificate')
    parser.add_argument('--preconditioner')
    args=parser.parse_args()
    if args.mode=='prepare':prepare(args)
    elif args.mode=='solve':solve_lp(args)
    elif args.mode=='certify':certify(args)
    else:preconditioner(args)
