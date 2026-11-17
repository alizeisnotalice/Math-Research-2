#!/usr/bin/env python3
"""Portable finite replay; no asymptotic inference and no frozen-file edits."""
import argparse,hashlib,json,os,pathlib,shutil,subprocess,sys

p=argparse.ArgumentParser()
p.add_argument('--repo-root',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--package-root',type=pathlib.Path)
a=p.parse_args(); repo=a.repo_root.resolve();out=a.output.resolve()
cert=repo/'free-boundary-grid-twin-peaks/certificates'
pkg=(a.package_root or cert/'r030').resolve();r29=cert/'r029'
assert pkg!=out and pkg not in out.parents and out not in pkg.parents
out.mkdir(parents=True,exist_ok=True)
commands=[]
def run(args,stdout=None):
    args=[str(x) for x in args]; commands.append(args)
    env=dict(os.environ,CLANG_MODULE_CACHE_PATH='/Volumes/data/CodexTemp/cache/clang-modules')
    if stdout:
        with open(stdout,'w') as f:subprocess.run(args,cwd=out,env=env,stdout=f,check=True)
    else:subprocess.run(args,cwd=out,env=env,check=True,stdout=subprocess.DEVNULL)
def py(rel,*args,stdout=None):run([sys.executable,'-B',pkg/rel,*args],stdout)
py('r030-short-tail-crt/audit_short_tail_crt.py','--raw',pkg/'r030-short-tail-crt/raw',
   '--pratt',cert/'r028/r028-punctured-layer/n7-center-crt/pratt-certificates.json',
   '--reference',r29/'r029-punctured-orbits-n7/verified/v24-verified.json',
   '--r29-source',r29/'r029-short-tail-mod/short_tail_mod.cpp','--output',out/'crt.json')
py('r030-integrated-tail/verify_prefix_tail_identity_n3.py',
   r29/'r029-tail-full-n5/tail-full-n3.json',out/'prefix.json',
   '--short-ie',r29/'r029-short-tail-ie/short_tail_ie_check.py')
py('r030-integrated-tail/audit_quadratic_envelope_finite.py','--r29',
   r29/'r029-tail-full-n5',stdout=out/'quadratic.json')
py('r030-integrated-tail/audit_local_qe8_n7_counterexample.py',out/'crt.json',out/'counterexample.json')
run(['clang++','-std=c++17','-O3',pkg/'r030-gap-geometry/gap_geometry_dp.cpp','-o',out/'gap_geometry_dp'])
run([out/'gap_geometry_dp','3','all'],out/'gap-n3-all.json')
py('r030-gap-geometry/audit_n3_joint.py',out/'gap-n3-all.json',cert/'r021-root-n7/reference/n3-roots.json')
for v in [0,1,2,6,7,12]:
    name=f'gap-n5-v{v}.json';run([out/'gap_geometry_dp','5',str(v)],out/name)
    assert json.loads((out/name).read_text())==json.loads((pkg/'r030-gap-geometry'/name).read_text())
py('r030-gap-geometry/summarize_n5_geometry.py','--reference',
   cert/'r018-root-moments/n5-d4-root-moments.json','--input-dir',out,'--output',out/'n5-gap-geometry-summary.json')
py('r030-gap-geometry/crosscheck_r29_tail.py','--geometry',out/'n5-gap-geometry-summary.json','--r29',
   r29/'r029-tail-full-n5','--output',out/'tail-crosscheck.json')
for rel in ['r030-corner-far-parent/audit_corner_insertion.py',
            'r030-corner-far/check_corner_geometry.py',
            'r030-spatial-geometry/audit_geodesic_path_kernel.py',
            'r030-integrated-tail/check_qe8_transfer_fractions.py']:
    dest=out/pathlib.Path(rel).name;shutil.copy2(pkg/rel,dest);run([sys.executable,'-B',dest])
receipt={'status':'passed','scope':'finite replay only; all six n5 holes independently recomputed',
         'commands':commands,'python':sys.version,'driver_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
(out/'replay-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS: CRT, prefix identity, finite tail checks, n3 exhaustive joint table, six n5 holes, corner insertion, path kernel, transfer Fraction checks')
