#!/usr/bin/env python3
"""Replay finite checks of the short-gap and spatial identities, outside Git."""
import argparse,hashlib,json,pathlib,shutil,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--certificates',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
a=p.parse_args();pkg=pathlib.Path(__file__).resolve().parent;out=a.output.resolve();cert=a.certificates.resolve()
assert pkg not in out.parents and out!=pkg and out not in pkg.parents
out.mkdir(parents=True,exist_ok=True)
groups={
 'r031-short-gap-identity':[
 ('enumerate_five_vertex_prefixes.py','five-vertex-prefix-enumeration.json'),
 ('enumerate_short_gaps_small_graphs.py','small-graph-short-gap.json'),
 ('check_rank3_identity_small_graphs.py','rank3-identity-small-graphs.json')],
 'r031-short-gap-spatial':[
 ('audit_joint_h2_h3.py','joint_h2_h3_audit.json'),
 ('audit_n3_rootpair_flux.py','n3_rootpair_flux_audit.json'),
 ('audit_n3_top_prefix_flux.py','n3_top_prefix_flux_audit.json'),
 ('audit_top5_polyomino_flux.py','top5_polyomino_flux_audit.json')]}
checks=[]
for group,items in groups.items():
    run_dir=out/group;run_dir.mkdir(exist_ok=True)
    for script,result in items:
        shutil.copy2(pkg/group/script,run_dir/script)
        args=[sys.executable,'-B',str(run_dir/script)]
        if script=='audit_joint_h2_h3.py':
            args+=['--geometry-dir',str(cert/'r030/r030-gap-geometry')]
        subprocess.run(args,cwd=run_dir,check=True,stdout=subprocess.DEVNULL)
        actual=json.loads((run_dir/result).read_text())
        expected=json.loads((pkg/group/result).read_text())
        assert actual==expected,(group,script)
        checks.append({'script':group+'/'+script,'complete_json_matches':True,
                       'result_sha256':hashlib.sha256((run_dir/result).read_bytes()).hexdigest()})
receipt={'status':'passed','scope':'finite integer crosschecks; not an asymptotic proof',
         'checks':checks,'all_subprocesses_waited':True,'python':sys.version,
         'driver_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
(out/'proof-replay-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS: all seven finite proof checks and complete output JSON comparisons')
