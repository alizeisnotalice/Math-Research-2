"""Bind complete original and checked degree-32 ledgers and freeze part files."""
from argparse import ArgumentParser
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path

from census32 import OLD, sha


def same_receipts(original,checked,candidates):
    with original.open() as a, checked.open() as b:
        start=json.loads(next(a))
        confirm=json.loads(next(b))
        identity=('start','stop','census_sha256','sources_sha256')
        assert all(start[k]==confirm[k] for k in identity)
        assert start['mode']=='generate' and confirm['mode']=='verify'
        assert confirm['witness_sha256']==sha(original)
        counts=Counter()
        unresolved=[]
        for n in range(start['start'],start['stop']):
            previous,independent=next(a),next(b)
            assert previous==independent
            row=json.loads(previous)
            assert row['index']==n and row['shape']==candidates[n]['shape']
            counts[row['status']]+=1
            if row['status']=='unit':
                counts['mod'+str(row['prime'])]+=1
                counts[row['method']]+=1
            else:
                assert row['status']=='unresolved'
                unresolved.append(row)
        end,other=json.loads(next(a)),json.loads(next(b))
        assert end['status']==other['status']=='complete' and end['counts']==other['counts']==dict(counts)
        assert end['start']==other['start']==start['start'] and end['stop']==other['stop']==start['stop']
        assert not a.read(1) and not b.read(1)
    return start,end,other,counts,unresolved


def archive(original,target):
    files=[]
    unpacked=sha256()
    with original.open('rb') as f:
        i=0
        while raw:=f.read(6_000_000):
            name=f'{original.stem}-{i:02d}.part.gz'
            encoded=gzip.compress(raw,compresslevel=6,mtime=0)
            assert len(encoded)<10*1024**2
            with (target/name).open('xb') as part:
                part.write(encoded)
            data=gzip.decompress((target/name).read_bytes())
            assert data==raw
            unpacked.update(data)
            files.append({'file':name,'raw_bytes':len(raw),'compressed_bytes':len(encoded),
                          'raw_sha256':sha256(raw).hexdigest(),'compressed_sha256':sha256(encoded).hexdigest()})
            i+=1
    assert unpacked.hexdigest()==sha(original)
    return files


def check_real(row,real):
    if row['shape']!=real['shape']:
        return False
    assert row['exact_system']['variable_indices']==real['variable_indices']
    terms=dict(row['exact_system']['early_zero_equations'])[real['row']]
    assert terms==real['exact_coefficients']
    assert real['strict_lower_bound']=='53/200'
    return True


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--census',type=Path,required=True)
    parser.add_argument('--census-audit',type=Path,required=True)
    parser.add_argument('--real',type=Path,required=True)
    parser.add_argument('--cubic',type=Path,required=True)
    parser.add_argument('--restart',type=Path,action='append',required=True)
    parser.add_argument('--units',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    original=json.loads(a.census.read_bytes())
    audit=json.loads(a.census_audit.read_bytes())
    real=json.loads(a.real.read_bytes())
    cubic=json.loads(a.cubic.read_bytes())
    restarts=[dict(json.loads(p.read_bytes()),sha256=sha(p)) for p in a.restart]
    assert all(q['generator_exit_code']==143 for q in restarts)
    assert real['status']==audit['status']==cubic['status']=='verified'
    assert audit['census_sha256']==cubic['census_sha256']==sha(a.census)
    assert [r['row24_cubic_coefficient'] for r in cubic['rows']]==[1,-1]
    labs=original['ledger']['candidates']
    def starting_index(p):
        with p.open() as rows:
            return json.loads(next(rows))['start']
    unitfiles=sorted(a.units.glob('raw-*.jsonl'),key=starting_index)
    assert unitfiles
    a.out.mkdir(exist_ok=False)
    accepted=Counter()
    slices=[]
    unresolved=[]
    resolved=[]
    source={}
    end=0
    checkpoints=set()
    for raw in unitfiles:
        checked=raw.with_name(raw.name.replace('raw-','check-',1))
        header,terminal,replay,counts,unknowns=same_receipts(raw,checked,labs)
        assert header['start']==end and header['census_sha256']==sha(a.census)
        end=header['stop']
        if 'checkpoint' in header:
            raw_sha=sha(raw)
            restart=next(q for q in restarts if q['reframed_sha256']==raw_sha)
            assert header['checkpoint']=={k:v for k,v in restart.items() if k not in ('reframed_sha256','sha256')}
            assert terminal['seconds'] is None
            checkpoints.add(raw_sha)
        accepted.update(counts)
        for name,digest in header['sources_sha256'].items():
            assert source.setdefault(name,digest)==digest
        for row in unknowns:
            if check_real(row,real):
                resolved.append({'shape':row['shape'],'index':row['index'],'strictly_positive_residual':real['strict_lower_bound']})
            else:
                unresolved.append(row)
        slices.append({'start':header['start'],'stop':header['stop'],'counts':terminal['counts'],
                       'generation_seconds':terminal['seconds'],'identity_replay_seconds':replay['seconds'],
                       'raw_sha256':sha(raw),'identity_sha256':sha(checked),'parts':archive(raw,a.out)})
    assert len(checkpoints)==len(restarts)
    assert end==len(labs)==accepted['unit']+accepted['unresolved']
    assert accepted['unresolved']==len(resolved)+len(unresolved)
    for supplement in (audit,real):
        for name,digest in supplement['sources_sha256'].items():
            assert source.setdefault(name,digest)==digest
    here=Path(__file__).resolve().parent
    previous=here.parent/'round066-gap7-next'
    source['verify_cubic32.py']=cubic['source_sha256']
    for restart in restarts:
        assert source.setdefault('seal_partial32.py',restart['source_sha256'])==restart['source_sha256']
    for name,digest in source.items():
        location=next(path/name for path in (here,previous,OLD) if (path/name).is_file())
        assert sha(location)==digest
    count=original['counts']
    assert count['all']==count['early_rejections']+count['late_rejections']+accepted['unit']+len(resolved)+len(unresolved)
    result={'status':'verified' if not unresolved else 'incomplete','degree':32,
            'counts':{'shapes':count,'modular':dict(accepted),'real_certificate':len(resolved),'unresolved_final':len(unresolved)},
            'census':{'sha256':sha(a.census),'parts':archive(a.census,a.out),'seconds':original['seconds']},
            'census_audit':{'sha256':sha(a.census_audit),'seconds':audit['seconds']},
            'real_check':{'sha256':sha(a.real),'resolved':resolved},
            'cubic_check':{'sha256':sha(a.cubic)},
            'restarts':restarts,'unresolved':unresolved,
            'slices':slices,'python':header['python'],'sympy':header['sympy'],'sources_sha256':source}
    with (a.out/'verification.json').open('x') as f:
        f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)


if __name__=='__main__':
    main()
