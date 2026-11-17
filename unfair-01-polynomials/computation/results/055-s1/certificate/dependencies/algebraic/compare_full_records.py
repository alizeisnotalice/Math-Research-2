"""Exact all-core comparison of independent direct-AB output to supplied schemas."""
import hashlib
import json
from pathlib import Path

OWN=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round055-s1-independent/full-v1')
PARENT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round055-degree8-entry')
P4=PARENT/'prefix-full-rows2to4-v2'
P7=PARENT/'prefix-full-rows5to7-v1'


def sha_bytes(chunks):
    h=hashlib.sha256()
    for name,raw in chunks:
        h.update(name.encode()); h.update(b'\0'); h.update(raw)
    return h.hexdigest()


def main():
    own_files=sorted((OWN/'cores').glob('*.json'))
    assert len(own_files)==1140
    c4=c7=0; own4=[]; par4=[]; own7=[]; par7=[]
    for op in own_files:
        p4=P4/'cores'/op.name; p7=P7/'cores'/op.name
        assert p4.is_file() and p7.is_file(), op.name
        a=json.loads(op.read_text()); b=json.loads(p4.read_text()); c=json.loads(p7.read_text())
        assert a['status']==b['status']==c['status']=='complete'
        assert a['summary']['T']==b['summary']['T']==c['summary']['T']
        assert a['summary']['N']==b['summary']['N']==c['summary']['N']
        assert a['row4_records']==b['records'], op.name
        assert a['row7_records']==c['records'], op.name
        # Compare all per-core counters, with the independent output split at row4.
        ac=a['summary']['counts']; bc=b['summary']['counts']; cc=c['summary']['counts']
        for key in ('b_box','s2_pass','s2_reject','c2_reject','row2_states','s3_pass','s3_reject','c3_reject','row3_states','s4_pass','s4_reject','c4_reject','row4_states'):
            assert ac.get(key,0)==bc.get(key,0),(op.name,key)
        assert len(a['row4_records'])==cc['stage4_records']
        for key in ('c5_reject','row5_states','e_box','s5_pass','s5_reject','s6_pass','s6_reject','c6_reject','row6_states','cross7_pass','g_box','s7_reject','s8_reject','row7_states'):
            assert ac.get(key,0)==cc.get(key,0),(op.name,key)
        c4+=len(a['row4_records']); c7+=len(a['row7_records'])
        own4.append((op.name,json.dumps(a['row4_records'],separators=(',',':')).encode()))
        par4.append((op.name,json.dumps(b['records'],separators=(',',':')).encode()))
        own7.append((op.name,json.dumps(a['row7_records'],separators=(',',':')).encode()))
        par7.append((op.name,json.dumps(c['records'],separators=(',',':')).encode()))
    assert c4==60134 and c7==67
    assert sha_bytes(own4)==sha_bytes(par4)
    assert sha_bytes(own7)==sha_bytes(par7)
    own=json.loads((OWN/'summary.json').read_text()); p4s=json.loads((P4/'summary.json').read_text()); p7s=json.loads((P7/'summary.json').read_text())
    assert own['status']==p4s['status']==p7s['status']=='complete'
    assert own['row4_states']==p4s['row4_states']==p7s['counts']['stage4_records']==60134
    assert own['row7_states']==p7s['row7_states']==67
    out={'status':'passed_all_1140_core_records_and_multiplicities','cores':1140,'row4_records':c4,'row7_records':c7,
         'row4_payload_sha256':sha_bytes(own4),'row7_payload_sha256':sha_bytes(own7),
         'independent_summary_sha256':hashlib.sha256((OWN/'summary.json').read_bytes()).hexdigest(),
         'parent_stage4_summary_sha256':hashlib.sha256((P4/'summary.json').read_bytes()).hexdigest(),
         'parent_stage57_summary_sha256':hashlib.sha256((P7/'summary.json').read_bytes()).hexdigest()}
    (OWN/'comparison.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps(out,sort_keys=True))


if __name__=='__main__': main()
