#!/usr/bin/env python3
"""Final byte/provenance binding for separately described R71 math audits."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,action='append',required=True);ap.add_argument('--review-root',type=Path,required=True);ap.add_argument('--archive-root',type=Path);ap.add_argument('--report',type=Path);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    packets=[];author=None
    for p in a.packet:
        d=json.loads(p.read_text());checked=[]
        topic=Path(d.get('root',next((str(parent) for parent in p.parents if parent.name=='planar-centered-square-weak11'),str(p.parents[2]))))
        for row in d['files']+d.get('existing_dependencies',[]):
            rel=row.get('destination','planar-centered-square-weak11/'+row.get('path',''))
            target=a.archive_root/rel if a.archive_root else Path(row['source']) if 'source' in row else topic/row['path']
            assert target.is_file() and sha(target)==row['sha256'],target
            if 'bytes'in row:assert target.stat().st_size==row['bytes']
            checked.append(dict(destination=rel,sha256=sha(target),bytes=target.stat().st_size))
            if rel.endswith('/proofs/r071-shared-outside-mass-and-residual-duality.md'):author=sha(target)
        for pin in d.get('source_pins',[]):
            rel='planar-centered-square-weak11/'+pin['path'].split('/planar-centered-square-weak11/',1)[1]
            target=a.archive_root/rel if a.archive_root else Path(pin['path'])
            assert sha(target)==pin['sha256']
            checked.append(dict(destination=rel,sha256=sha(target),bytes=target.stat().st_size,role='inherited source'))
        ext=[r['external_input'] for r in d.get('external_oracle_input_bindings',[])]+[d[k] for k in ('external_robust_input','external_robust_compared_claim') if k in d]
        for pin in ext:
            rel='planar-centered-square-weak11/'+pin['path']
            target=a.archive_root/rel if a.archive_root else Path(d['producer_workspace'])/pin['path']
            assert sha(target)==pin['sha256'] and target.stat().st_size==pin['bytes']
            checked.append(dict(destination=rel,sha256=sha(target),bytes=target.stat().st_size,role='external exact input/compared claim'))
        packets.append(dict(manifest=p.name,sha256=sha(p),checked_rows=checked))
    assert author=='d61554809a35ef393b735985fea4c8732a463df47ba51ff69ef837ca86fc9770'
    r=a.review_root;own=[]
    for stem,script in [('boundary-controls','review_r071_joint_measure_boundaries.py'),('production-certificates','review_r071_production_certificates.py')]:
        result=r/'computation/results/round071-independent'/f'{stem}.json';receipt=result.with_name(f'{stem}-receipt.json');d=json.loads(result.read_text());rec=json.loads(receipt.read_text())
        assert rec['exit_status']==0 and rec['owned_pid_absent'] and not rec['stderr']
        assert rec['output_sha256']==sha(result)
        assert d['source_sha256']==rec['source_sha256']==sha(r/'computation'/script)
        if stem=='production-certificates':
            cmd=rec['command']
            for name,h in d['input_sha256'].items():
                original=Path(cmd[cmd.index('--'+name.replace('_','-'))+1])
                if a.archive_root:
                    suffix=str(original).split('/planar-centered-square-weak11/',1)[1];target=a.archive_root/'planar-centered-square-weak11'/suffix
                else:target=original
                assert sha(target)==h,(name,target)
        own.append(dict(result_destination='planar-centered-square-weak11/'+str(result.relative_to(r)),result_sha256=sha(result),receipt_sha256=sha(receipt),checker_sha256=d['source_sha256'],recorded_pid=rec['pid'],recorded_exit0_absent=True))
    result=dict(status='PASS_R71_FINAL_BYTE_BINDINGS_ONLY',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(Path(__file__)),packets=packets,author_proof_sha256=author,independent_result_chains=own,report=None if a.report is None else dict(destination='planar-centered-square-weak11/rounds/071.md',sha256=sha(a.report),bytes=a.report.stat().st_size),scope='Byte/provenance binding only; proof reading, actual independent exact controls, producer code review, source-bound author controls, and root MATLAB/Wolfram work remain distinct. No LP or old math calculation rerun by this binding check.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],packets=len(packets),checked_packet_rows=sum(len(p['checked_rows']) for p in packets),own_result_chains=len(own))))


if __name__=='__main__':main()
