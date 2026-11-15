#!/usr/bin/env python3
"""Portable exact R66 archive replay. Read only pinned archive-relative dependencies.

--output must name a NEW JSON file below /Volumes/data/CodexTemp.
Requires Python 3.10+ standard library and all32 pinned archive dependencies.
Frozen audit code is imported; no producer or original workspace is accessed.
"""
import argparse,hashlib,importlib.util,json,os,sys,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
sys.dont_write_bytecode=True

PINS={'planar-centered-square-weak11/computation/results/round066/theta/r64-final-theta.json': '9a68cc30dac1cd030da2ae8b0b2673fb963c5f32bd86473c4a81d4fa6d05eb9e', 'planar-centered-square-weak11/computation/results/round066-erosion/r64-final-erosion.json': '0dc19c78d01b341f94d7c94648e5531850e0c2df9a088a9683adbb458fb3af11', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-00.json': 'f292dabea972784f9332c52d3a2f4cfebcc0e295699417f73b6477c0d2c3982c', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-00-erosion.json': '9a11b6bb5a9789b7f622f25ba57d1dc12073e47b2b2c791eea28684f06a1ea03', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-01.json': 'ecdcf3a0a263f5c8ecb06f1c9e550333b3d3e4356f66eb40efc9c904a9c6868b', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-01-erosion.json': 'f9c1fae9a1bb8a91965b30f42f6cc8e1633227240a4879b5a20a9b11bb4f8188', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-02.json': '921751c7e3de5efddb1c55e57cacfadd22cf982b7a631db71ba1bd9ae66cb70e', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-02-erosion.json': '7be956ce7d9b0f480d2f6c2f4209e673b42bf7802277c6f56a8ac80a0d014673', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-03.json': 'b72b159e1ff6c4d96c8626f97a198816f2a065f37713b63a5ededff4343398d4', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-03-erosion.json': '216fdd0b2115288794a549a70ff2fb0536fd60b5f6f505203e25201b1bbbf260', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-04.json': 'c6f4df7ef1b5eca533888dadcb6c9878c899e25ef1e59fe01131e1450ae61ef6', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-04-erosion.json': '60897162b601475ead95eb99f54030b5e634de02be914838eba78c9b73b12bd6', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-05.json': '4f55b93ccdeb623578a395aeb7b32efa4197822e753a4920f5de3ce3d7d879d6', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-05-erosion.json': 'fae33889d7123fe15acc95c7f75a24fbb1082d6a55ebad12f8171e2bae0ce5c1', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-06.json': '4b7c38f6636f799dbe29811f2f6ec67ee5fecd58d182651d702cda38fbeab50d', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-06-erosion.json': '8579381e04539bbc1d35d2512a55bb636551d12309cfe5325ac6b93ac8589d94', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-07.json': '4671d4a0f669f01b28c22adaac3545ab21adbd2f89c3c68330cbea5d12a845b6', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-07-erosion.json': 'a91598d75ca7adf920c55a1a421b8dcc5bd5253c7d95c4545551c5992535e89b', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-08.json': '782e622f7e078255cdcbd0bfe0ccef12abb8fca08377f581e1836790f3fc84ed', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-08-erosion.json': 'e1aa14d3eea73d22fad649f8c3db4af27ce52838b23e1ea5ad234f91b057706a', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-09.json': '822af31eed800d1c24d9e0291bfd8817730e91400a885180610c50159865cab6', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-09-erosion.json': 'a66078a5ba9548a1bac8e170f4a04e15293a60a8509f7c47e6cb9b747d72ff28', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-10.json': '96d8e4ebbb282be4a36d2d29c785fffc0b3e05e6f78a405ed422d3ad8cd628f0', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-10-erosion.json': '8ff082787b48fd9a05bb7ee8795e485c6baaf2a2208b777eabf4831437788019', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-11.json': '885e80ec4172fd21c14c70ecbdd138fc4063c246c72657adabfd8edd41d68542', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-11-erosion.json': 'ad4b9312d06f8e12ff47349ff7131047db7dcc6c6da903c4a4e45fb56dd27526', 'planar-centered-square-weak11/computation/results/round066-erosion/five-source-control/budget.json': '42851809b3bb86c42a0981cc429287561dd67844357abb932a3051e92c242408', 'planar-centered-square-weak11/computation/results/round066-erosion/five-source-control/erosion.json': '113656907e2298a5f6d07695ebccd72b0e93c3894bb0b85a74b6f3edc7333c3d', 'planar-centered-square-weak11/computation/results/round066-erosion/three-source-fat-control/budget.json': '1762e58524542dc332073522b2bfc3a82c28c87b981c9b60bfdaf5c0d7f900d0', 'planar-centered-square-weak11/computation/results/round066-erosion/three-source-fat-control/erosion.json': '68389678b1cadc1a3938c002f4a516cef0cad0b8dd183fe574f429665f24ab50', 'planar-centered-square-weak11/computation/review_r066_grid_budgets.py': 'c6524c436d8efd15a1b941c3e37448608622f6723fc51e3325a6eaa44aa159b1', 'planar-centered-square-weak11/computation/review_r064_mass_variation.py': '267b430e167801314bba81200250945bedcb1473918992f57ac03c1feb94cad4'}
PAIRS=[('R64', 'planar-centered-square-weak11/computation/results/round066/theta/r64-final-theta.json', 'planar-centered-square-weak11/computation/results/round066-erosion/r64-final-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-00.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-00-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-01.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-01-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-02.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-02-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-03.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-03-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-04.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-04-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-05.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-05-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-06.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-06-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-07.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-07-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-08.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-08-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-09.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-09-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-10.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-10-erosion.json'), ('R63', 'planar-centered-square-weak11/computation/results/round066/theta/stratum-11.json', 'planar-centered-square-weak11/computation/results/round066-erosion/stratum-11-erosion.json'), ('THIN_CONTROL', 'planar-centered-square-weak11/computation/results/round066-erosion/five-source-control/budget.json', 'planar-centered-square-weak11/computation/results/round066-erosion/five-source-control/erosion.json'), ('FAT_CONTROL', 'planar-centered-square-weak11/computation/results/round066-erosion/three-source-fat-control/budget.json', 'planar-centered-square-weak11/computation/results/round066-erosion/three-source-fat-control/erosion.json')]

CORE_REL='planar-centered-square-weak11/computation/review_r066_grid_budgets.py'
BUDGET_SHA='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
EROSION_SHA='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
GEOMETRY_SHA='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def require(ok,message):
    if not ok:raise RuntimeError(message)
def boot(core_path):
    global AUDIT
    s=importlib.util.spec_from_file_location('frozen_r66_independent_grid',core_path)
    AUDIT=importlib.util.module_from_spec(s);s.loader.exec_module(AUDIT)
    AUDIT.init_worker()
def evaluate(job):return AUDIT.evaluate_case(job)
def run(args):
    require(sys.flags.optimize==0,'Assertions required: do not use Python -O')
    require(1<=args.workers<=12,'workers must be between1 and12')
    root=args.archive_root.resolve()
    if root.name=='planar-centered-square-weak11':root=root.parent
    out=args.output.expanduser().resolve()
    require(Path('/Volumes/data').is_mount(),'/Volumes/data must be mounted')
    require(out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve()),'--output must be below /Volumes/data/CodexTemp')
    require(out.suffix=='.json','--output must name a JSON file')
    rec=out.with_suffix('.receipt.json')
    require(not out.exists()and not rec.exists(),'Output/receipt already exists: choose new paths')
    out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();receipt=dict(status='RUNNING',entry_sha256=digest(Path(__file__)),python=sys.version,
        command=sys.argv,archive_root=str(root),output=str(out),pid=os.getpid(),workers=args.workers)
    try:
        pins=[]
        for rel,expected in sorted(PINS.items()):
            p=(root/rel).resolve();require(p.is_relative_to(root),'Dependency escapes archive root')
            require(p.is_file(),f'Missing archived dependency: {rel}')
            actual=digest(p);require(actual==expected,f'Dependency SHA mismatch: {rel}')
            pins.append(dict(path=rel,sha256=actual,bytes=p.stat().st_size))
        jobs=[]
        for role,brel,erel in PAIRS:
            bp,ep=root/brel,root/erel;b=json.loads(bp.read_text());e=json.loads(ep.read_text())
            require(b['status']==e['status']=='PASS','Non-PASS archived input')
            require(b['source_sha256']==BUDGET_SHA and e['source_sha256']==EROSION_SHA,'Producer source chain mismatch')
            require(b['geometry_module_sha256']==e['geometry_sha256']==GEOMETRY_SHA,'Geometry source chain mismatch')
            require(e['budget_result_sha256']==digest(bp),'Budget/erosion pair mismatch')
            emap={v['id']:v for v in e['cases']}
            require(len(emap)==len(e['cases'])and set(emap)=={v['id']for v in b['cases']},'Case coverage mismatch')
            jobs.extend((role,c,emap[c['id']])for c in b['cases'])
        core=str(root/CORE_REL)
        if args.workers==1:
            boot(core);rows=[evaluate(j)for j in jobs]
        else:
            with ProcessPoolExecutor(max_workers=args.workers,initializer=boot,initargs=(core,))as pool:rows=list(pool.map(evaluate,jobs))
        blocks=[b for r in rows for b in r['blocks']];old=[b for r in rows if r['role']=='R63'for b in r['blocks']]
        stats=dict(cases=len(rows),blocks=len(blocks),R63_cases=sum(r['role']=='R63'for r in rows),R63_blocks=len(old),
            R63_budget_pass=sum(b['budget_pass']for b in old),R63_strong_pass=sum(b['strong_pass']for b in old),
            R63_strong_strict_fail=sum(b['strong_strict_fail']for b in old),high_labels=sum(b['high_labels']for b in blocks),
            total_compressed_grid_cells=sum(b['compressed_grid_cells']for b in blocks))
        require(tuple(stats[k]for k in ('cases','blocks','R63_cases','R63_blocks','R63_budget_pass','R63_strong_pass','R63_strong_strict_fail'))==(75,202,72,199,199,181,18),'Final replay coverage mismatch')
        fat=[r for r in rows if r['role']=='FAT_CONTROL'];require(len(fat)==1,'Fat control coverage')
        f=fat[0];b=f['blocks'][0]
        from fractions import Fraction
        require(f['all_maximum_masks']==[3]and f['connected']and f['whole_nonterminal']and b['fat_both_axes'],'Fat counterexample domain')
        require(Fraction(b['budget']['lower'])>Fraction(7,2)and Fraction(b['deletion']['upper'])<1,'Fat counterexample signs')
        elapsed=time.perf_counter()-start
        result=dict(status='PASS_PORTABLE_R66_ARCHIVE_REPLAY',entry_sha256=receipt['entry_sha256'],python=sys.version,
            workers=args.workers,elapsed_seconds=elapsed,statistics=stats,dependency_pins=pins,cases=rows,
            scope='Independent exact2D difference-grid replay of all202 archived blocks. Frozen mathematical audit pins unchanged. Fixed theta envelope refuted; no general Delta/C2 claim.')
        out.write_text(json.dumps(result,indent=2)+'\n')
        receipt.update(status=result['status'],elapsed_seconds=elapsed,statistics=stats,output_sha256=digest(out),dependency_count=len(pins))
        rec.write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({k:v for k,v in receipt.items()if k not in ('command',)},indent=2))
    except BaseException as exc:
        receipt.update(status='FAIL_PORTABLE_R66_ARCHIVE_REPLAY',elapsed_seconds=time.perf_counter()-start,error=str(exc))
        rec.write_text(json.dumps(receipt,indent=2)+'\n');raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[2],help='Repository root or its planar-centered-square-weak11 topic directory')
    p.add_argument('--output',type=Path,required=True,help='New .json result path below /Volumes/data/CodexTemp; paired .receipt.json written beside it')
    p.add_argument('--workers',type=int,default=12)
    run(p.parse_args())
if __name__=='__main__':main()
