"""Certify R221's row-specific point displacement on every allowed angle."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
import sys
import tarfile
import tempfile
import time

from flint import arb,ctx,__version__ as flint_version

PREC=160
DEN=10**16

def sha(b):
    return sha256(b).hexdigest()

def as_arb(value):
    v=F(value)
    return arb(v.numerator)/v.denominator

def rounded(x, upper):
    z=(x*DEN).upper().ceil() if upper else (x*DEN).lower().floor()
    return F(int(z.unique_fmpz()),DEN)

def check_paths(comp,work):
    public=comp/'six-input-083745'
    archive_path=public/'volumes/native45-receipts.tar.gz'
    with tarfile.open(archive_path,'r:gz') as archive:
        index={path:(size,digest) for path,size,digest in json.load(archive.extractfile('MANIFEST.json'))}
        inputs=('inputs/ordinary/produce_p30_083745.py','inputs/ordinary/finite-point-lift-083744.json',
                'inputs/ordinary/scope.json')
        raws=[archive.extractfile(name).read() for name in inputs]
    for name,raw in zip(inputs,raws):
        assert (len(raw),sha(raw))==index[name]
    ordinary,last,scope=raws
    last=json.loads(last)
    scope=json.loads(scope)
    paths=('finite-point-leaf-input.json','finite-point-global-input.json','tomar_six_083730_guard.cpp',
           'tomar_p30_mesh_guard.py','tomar_p30_scaled_material.py','replay_scaled_p30_high_leaf.py',
           'replay_finite_point_083732.py','certify_finite_point_leaf.py','finite-point-lift-083742.json',
           'certify_tomar_six_high.py',
           'six-input-083745/records/integration-audit.json','six-input-083745/records/dense45.json')
    source={k:(comp/k).read_bytes() for k in paths}
    hashes={k:sha(v) for k,v in source.items()}
    for name,key74,key72 in (('finite-point-leaf-input.json','point_sha256','point_sha256'),
                             ('finite-point-global-input.json','bank_sha256','bank_sha256'),
                             ('tomar_six_083730_guard.cpp','weighted_sha256','weighted_sha256')):
        assert last[key74]==json.loads(source['finite-point-lift-083742.json'])[key72]==hashes[name]
    assert len(scope['leaves'])==3071 and len(scope['deferred'])==202
    tally=json.loads(source['six-input-083745/records/integration-audit.json'])
    assert tally['ordinary_closed']==3053
    assert tally['ordinary_rescued']==18 and tally['five_only']==196
    d=json.loads(source['six-input-083745/records/dense45.json'])
    assert d['count']==24 and d['all_full_field_equal_and_closed']
    assert b'bank=PointProof(read("finite-point-leaf-input.json"),read("finite-point-global-input.json")["cells"],100000000008865717121)' in ordinary
    assert b'support=bank.presolve(box)' in ordinary and b'len(support[3])-1,2,out' in ordinary
    dense=source['replay_scaled_p30_high_leaf.py']
    for s in (b'proof=PointProof(cloud,bank[\'cells\'],100000000008865717121)',
              b'support=piecewise_support(proof,box,Q,U,ceil)',b'len(support[3])-1,2,st',
              b'variant==\'mesh64-dense\'',b'mesh.dense_material = scaled_material'):
        assert s in dense,s
    assert b'return proof.directions,proof.disk,np.ascontiguousarray(groups,dtype=np.int64).ravel(),proof.pointers,proof.triples' in source['tomar_p30_mesh_guard.py']
    assert b'return (self.directions,self.disk,np.array(ends,dtype=np.int64).ravel(),self.pointers,self.triples)' in source['certify_tomar_six_high.py']
    sys.path.insert(0,str(comp))
    from tomar_p30_mesh_guard import write_kernel
    with tempfile.TemporaryDirectory(prefix='row-scope-',dir=work) as d:
        p=Path(d)
        (p/'fan-exposed.hpp').write_text('hi[index]>lo[index] && (hi[index]+lo[index])%2==0')
        mesh=write_kernel(p,comp/'tomar_six_083730_guard.cpp').read_text()
    ordinary_src=source['tomar_six_083730_guard.cpp'].decode()
    row_logic=lambda s:s[s.index('        bool needs_collective=false;'):s.index('    bool stop(')]
    assert row_logic(ordinary_src)==row_logic(mesh)
    assert ordinary_src.count('const L* q=triples+8*k;')==2
    assert ordinary_src.count('for(int j=0;j<3;++j)')>=1
    return (json.loads(source['finite-point-leaf-input.json']),
            json.loads(source['finite-point-global-input.json']),
            dict(hashes=hashes,archive_sha256=sha(archive_path.read_bytes()),
                 archive_members_sha256={name:sha(data) for name,data in zip(inputs,raws)},
                 ordinary_whole_leaf_success=3053,dense_whole_leaf_success=24,
                 unused_p30_five_only=196,row_logic_matches_generated_dense=True))

def peak(A,B,left,right):
    r"""Enclose max(A cos(t)+B sin(t)) in a closed cell of width less than pi."""
    cl,sl=left;cr,sr=right
    f0,f1=A*cl+B*sl,A*cr+B*sr
    d0,d1=-A*sl+B*cl,-A*sr+B*cr
    edge=(max(f0.lower(),f1.lower()),max(f0.upper(),f1.upper()))
    if d0<0 or d1>0:
        return (*edge,'endpoint')
    v=(arb(A*A+B*B)).sqrt()
    if d0>0 and d1<0:
        return (v.lower(),v.upper(),'interior')
    return (edge[0],max(edge[1],v.upper()),'boundary_interval')

def enum(cloud,bank,source):
    from certify_finite_point_leaf import point_tables
    ctx.prec=PREC
    P,V=cloud['points'],cloud['directions']
    cells=bank['cells']
    dia,pointers,flat=point_tables({**cloud,'cells':cells})
    assert (len(P),len(V),len(cells),pointers[-1])==(30,105,384,2493)
    assert dia==100000000008865717121
    assert {tuple(p) for p in P}=={(x,-y) for x,y in P}
    arc=[((2*arb.pi()*j/len(cells)).cos(),(2*arb.pi()*j/len(cells)).sin()) for j in range(len(cells)+1)]
    dot={(i,j):(a*x+b*y,b*x-a*y) for i,(a,b) in enumerate(V) for j,(x,y) in enumerate(P)}
    lens=[arb(x*x+y*y).sqrt() for x,y in V]
    rows=[]
    choices=0
    peak_cases={}
    distinct=set()
    for i,entries in enumerate(cells):
        theta=(2*arb.pi()*i/len(cells)).union(2*arb.pi()*(i+1)/len(cells))
        rc,rs=theta.cos(),theta.sin()
        available={}
        for direction in {j for row in entries for j in row[:3]}:
            whole=[A*rc+B*rs for j in range(len(P)) for A,B in (dot[direction,j],)]
            floor=max(x.lower() for x in whole)
            available[direction]=[j for j,x in enumerate(whole) if x.upper()>=floor]
            assert available[direction]
        assert len(entries)==pointers[i+1]-pointers[i]
        for j,row in enumerate(entries):
            k=pointers[i]+j
            ids,material=row[:3],row[3:]
            sampled=flat[8*k:8*k+8]
            assert sampled[:3]==ids
            weights=sampled[3:6]
            assert all(x>0 for x in weights)
            assert all(sum(w*V[l][d] for w,l in zip(weights,ids))==0 for d in (0,1))
            denom=cloud['denominator']*sum(w*lens[v] for w,v in zip(weights,ids))
            distinct.add(tuple(ids))
            def value(js):
                A,B=(sum(w*dot[l,p][ab] for w,l,p in zip(weights,ids,js)) for ab in (0,1))
                mi,ma,how=peak(A,B,arc[i],arc[i+1])
                peak_cases[how]=peak_cases.get(how,0)+1
                return (mi/denom).lower(),(ma/denom).upper()
            projection=[value(js) for js in product(*(available[l] for l in ids))]
            choices+=len(projection)
            fixed=value(material)
            retained=(max(v[0] for v in projection),max(v[1] for v in projection))
            assert retained[1]>=fixed[0]
            rows.append(dict(cell=i,row=j,ids=ids,witness=material,
                             lo=rounded(retained[0],False),hi=rounded(retained[1],True),
                             point_lo=rounded(fixed[0],False),point_hi=rounded(fixed[1],True)))
    assert len(rows)==pointers[-1]
    top=lambda k: sorted(rows,key=lambda r:r[k],reverse=True)[:8]
    global_bounds={k:(str(max(r['lo']for r in rows)),str(max(r['hi']for r in rows))) if k=='support'
                   else (str(max(r['point_lo']for r in rows)),str(max(r['point_hi']for r in rows)))
                   for k in ('support','fixed_witness')}
    q=F(global_bounds['support'][1])
    assert q<F('.528619')
    return dict(complete=True,precision=ctx.prec,python=platform.python_version(),flint=flint_version,
                cells=len(cells),row_angle_instances=len(rows),distinct_triples=len(distinct),binary_rows=0,
                represented_handednesses=2,incidence_angle_width='2*pi/384',
                extreme_support_combinations=choices,trig_cases=peak_cases,
                worst_intervals=global_bounds,top_support=top('hi'),top_fixed_witness=top('point_hi'),
                shrink_radius_was='.528619',source=source)

def main():
    p=ArgumentParser(description=__doc__)
    p.add_argument('comp',type=Path)
    p.add_argument('work',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    assert __debug__
    start=time.monotonic()
    cloud,bank,scope=check_paths(a.comp.resolve(),a.work.resolve())
    ans=enum(cloud,bank,scope)
    ans['elapsed_seconds']=time.monotonic()-start
    ans['source']['verifier_sha256']=sha(Path(__file__).read_bytes())
    a.out.write_text(json.dumps(ans,indent=2,default=str)+'\n')
    print(ans['worst_intervals'],ans['elapsed_seconds'],flush=True)

if __name__=='__main__':
    main()
