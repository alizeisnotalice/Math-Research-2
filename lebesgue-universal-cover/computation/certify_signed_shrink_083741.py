"""Check the R140 contracted material, the signed endpoint kernel and the .83741 shallow caps."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
from flint import arb, ctx


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("computation","tree","table32","table64","work","cache"):
        p.add_argument(name,type=Path)
    p.add_argument("--margin-kit",type=Path)
    p.add_argument("--base-kit",type=Path)
    args=p.parse_args()
    assert __debug__ and (args.base_kit is None)==(args.margin_kit is None)
    comp=args.computation.resolve()
    sys.path.insert(0,str(comp))
    from certify_polygon_083731 import LAMBDA, ORDERS, ar, lower_radius, verified_table, verify_core_table
    from certify_tomar_six_high import SourceTree
    from select_tomar_083731_scores import radius_rows
    from certify_shallow_cap_transfer import certify as verify_previous_shallow
    import replay_tomar_083731_scored as baseline
    from replay_tomar_signed_translation_screen import compile_screen
    ctx.prec=160
    work=args.work.resolve()
    work.mkdir(exist_ok=False)
    seed=load(comp/'tomar-083741-seed.json')
    stage=load(comp/'tomar-083741-stage-ledger.json')
    screen=load(comp/'tomar-083741-signed-screen-audit.json')
    ledger=load(comp/'tomar-083741-score-layers-audit.json')
    terminal=load(comp/'finite-point-lift-083741.json')
    roots=load(comp/'tomar-083741-retained.json')
    assert set(screen['files']).issuperset({'check_tomar_signed_translation_identity.py','tomar_signed_translation_screen.inc'})
    for name,h in screen['files'].items():
        assert sha(comp/name)==h
    algebra=json.loads(subprocess.check_output([sys.executable,"-I","-B",str(comp/'check_tomar_signed_translation_identity.py')]))
    assert algebra==screen['algebra']
    assert all(F(row['target'])==F('.83741') for row in (seed,stage,screen,ledger,terminal))
    assert stage['complete'] and screen['complete'] and stage['parts']==1024
    assert ledger['retained']==roots and sha(comp/'tomar-083741-stage-ledger.json')==ledger['stage_sha256']==terminal['outer_audit_sha256']
    assert sha(comp/'tomar-083741-retained.json')==terminal['retained_sha256']
    assert sha(comp/'tomar-083741-seed.json')==ledger['first_seed_sha256']==stage['seed_receipt_sha']
    assert len(stage['audits'])==4 and stage['audits'][-1]['remaining']==2083
    code_work=work/'signed';code_work.mkdir()
    _,signed_sha,_=compile_screen(baseline,baseline.kernel,comp,code_work,args.cache)
    assert signed_sha==stage['source_sha']==screen['geometry_source_sha256']
    tree=SourceTree(args.tree)
    assert sha(args.tree/'SHA256SUMS')==seed['tree_manifest_sha256']==screen['manifest_2048']
    base=load(comp/'tomar-six-high-input.json')
    assert terminal['base_sha256']==sha(comp/'tomar-six-high-input.json')
    assert sha(comp/'tomar_six_083730_guard.cpp')==terminal['weighted_sha256']
    assert sha(comp/'finite-point-leaf-input.json')==terminal['point_sha256'] and sha(comp/'finite-point-global-input.json')==terminal['bank_sha256']
    radial=load(comp/'five-body-radial-affine.json')
    assert sha(comp/'five-body-radial-affine.json')==seed['provenance']['radial_certificate']
    assert radius_rows(radial,F(seed['radial_target']))==seed['radii']
    table_shas={}
    for n,path in [(32,args.table32),(64,args.table64)]:
        cnt,pts,stamp=verified_table(path,tree,n,base['tables'][str(n)])
        verify_core_table(tree,cnt,pts,arb.pi())
        table_shas[str(n)]=stamp
        if n==32:
            assert stamp==radial['source_table_sha256']
            f=comp/'five-body-radial-material.npz'
            assert sha(f)==radial['table_sha']
            with np.load(f,allow_pickle=False) as row:
                assert np.array_equal(row['counts'],cnt) and np.array_equal(row['points'],pts[0])
    G=F(seed['target_integer'],2**57)
    assert G==F(stage['integer_goal'],2**57)>F('.83741')
    assert all(F(a['integer_radius'],2**44)<F('.069') and F(a['integer_area_lower'],2**57)>G for a in seed['radii'])
    s=LAMBDA
    assert s==F('.999999035')
    r,R=F(1,2),F('.703')
    for i,n in enumerate(ORDERS[1:]):
        k=3*i
        dist=min(F(tree.root['hi'][k]),F(tree.root['hi'][k+1]),-F(tree.root['lo'][k]),-F(tree.root['lo'][k+1]))
        ins=ar(s)*lower_radius(n,arb.pi())
        val=arb.pi()*(ar(s*r)**2+ins**2)/2+ar(dist)*(ar(s*r)+ins)
        assert val>ar(G)
    prior=verify_previous_shallow(work/'rational-points.json')
    factor=F(prior['common_scale'])
    assert sha(work/'rational-points.json')==prior['sha256']==sha(comp/'shallow-cap-points.json')
    oldfive=load(comp/'five-global-minimum-single-worker-replay.json')
    assert oldfive['complete'] and arb(oldfive['minimum_interval'])>ar(F('.837295'))
    radius_area=(arb.pi()-(ar(s*r)/ar(R)).acos())*ar(s*r)**2+ar(s*r)*(ar(R)**2-ar(s*r)**2).sqrt()
    assert radius_area>ar(G)
    height=F(prior['cap_height']);excess=G-s*s*F('.837295');area=s*r*height**2/R
    assert height==F('.013') and area>excess and F(prior['comparison_margin'])>0
    loss=F(prior['four_loss'])
    assert loss==F(22,7)*(2*R*F('.6')*(1-s)+(F('.6')*(1-s))**2)
    assert F('.837413')-loss>G
    pure=set(map(tuple,terminal['five_only']))
    ordinary={tuple(row[:2]) for row in terminal['rows']}
    keys={(int(k),val) for k,v in roots.items() for val in v}
    assert len(ordinary)==2049 and len(pure)==34 and not ordinary&pure and ordinary|pure==keys
    result=dict(complete=False,signed_kernel_sha256=signed_sha,shrunk_table_sha256=table_shas,
                terminal_sha256=sha(comp/'finite-point-lift-083741.json'),curve_scale=str(s),
                integer_goal=str(G),cap_height=str(height),max_allowed_excess=str(excess),
                cap_area_at_height=str(area),loss=str(loss),transfer_margin=str(F('.837413')-loss-G),
                shallow_cap_margin=prior['comparison_margin'],coordinate_sha256=prior['sha256'],
                rational_scale=str(factor),conditional_floor=str(factor**2*G))
    if args.margin_kit:
        assert sha(args.base_kit)==terminal['five_kit_sha256']
        strong=[]
        with ZipFile(args.base_kit) as original, ZipFile(args.margin_kit) as archive:
            for i,j in sorted(pure):
                stem=f'expected/part{i}-pos{j}'
                a=json.loads(original.read(stem+'.json'))
                b=json.loads(archive.read(stem+'.json'))
                for field in ('root','upstream','original_tree_sha256','original_header_sha256'):
                    assert a[field]==b[field]
                assert (b['part'],b['pos'])==(i,j) and F(b['threshold'])>=F('.837413') and not b['result']['remaining']
                assert sha256(archive.read(stem+'.bin')).hexdigest()==b['trace_sha256']
                strong.append([i,j,b['threshold'],b['result']['vis'],b['trace_sha256']])
        floor=factor**2*G
        assert floor>F('.83741123225')
        result.update(complete=True,high_margin_sha256=sha(args.margin_kit),strong=strong,
                      area_floor=str(floor),decimal_floor='.83741123225')
    (work/'replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='strong'},indent=2))


if __name__=='__main__':
    main()
